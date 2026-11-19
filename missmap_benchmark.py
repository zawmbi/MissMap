"""
MissMap LLM benchmark (Aim 3): run several models k times per species, then
score them against authority reference sets.

Run (calls are logged and cached in the JSONL; re-running never re-bills a
finished call):

  python missmap_benchmark.py run --species-file sample.csv \
      --models anthropic:claude-opus-5-5 anthropic:claude-sonnet-5-5 \
               anthropic:claude-haiku-5-5 groq:llama-3.3-70b-versatile \
               groq:llama-3.1-8b-instant anthropic:claude-sonnet-5-5+web \
      --runs 5 --log missmap_results/benchmark_calls.jsonl

Score:

  python missmap_benchmark.py score --log missmap_results/benchmark_calls.jsonl \
      --wcvp wcvp_names.csv --wfo classification.csv [--gbif] [--gold gold.csv] \
      --out missmap_results/benchmark

Outputs (in --out): per_model.csv, per_species.csv, pairwise_mcnemar.csv,
cochran_q.json. Every reference set is scored separately (per authority,
consensus = named by >=2 authorities, union, and the gold set if given),
because there is no single true synonym list.
"""

import argparse
import csv
import json
import os
import sys
from collections import defaultdict
from itertools import combinations

from missmap_authorities import load_wcvp, load_wfo, load_dwc, GBIFBackbone, authority_synonyms
from missmap_llm import make_provider, query_synonyms, synonym_names, LLMLog
from missmap_names import parse_name, normalise_key
from missmap_stats import (wilson_ci, mcnemar_exact, cochrans_q, holm, mean_pairwise_jaccard,
                           precision_recall, bootstrap_ci)


def read_species(path=None, inline=None):
    if inline:
        return [s.strip() for s in inline.split(",") if s.strip()]
    with open(path, encoding="utf-8") as f:
        first = f.readline()
        f.seek(0)
        if "," in first or "species" in first.lower():
            rdr = csv.DictReader(f)
            col = next((c for c in rdr.fieldnames if c.lower() in ("species", "taxon_name", "name")),
                       rdr.fieldnames[0])
            return [r[col].strip() for r in rdr if r.get(col, "").strip()]
        return [l.strip() for l in f if l.strip()]


def cmd_run(args):
    species = read_species(args.species_file, args.species)
    log = LLMLog(args.log)
    providers = [make_provider(m, temperature=args.temperature, effort=args.effort) for m in args.models]
    total = len(species) * len(providers) * args.runs
    done = 0
    for prov in providers:
        errors = 0
        for sp in species:
            for k in range(args.runs):
                rec = query_synonyms(prov, sp, known=None, run_index=k, log=log)
                errors += bool(rec.get("error"))
                done += 1
                if done % 25 == 0 or done == total:
                    print(f"[{done}/{total}] {prov.name}:{prov.model} ({prov.mode}) {sp} run {k}")
        print(f"[INFO] {prov.name}:{prov.model} ({prov.mode}): {errors} failed call(s).")
    print(f"[INFO] Log: {args.log}")


def _key(name):
    p = parse_name(name)
    return p.key() if p.parsed else normalise_key(name)


def build_references(species_list, checklists, gold=None):
    """
    {species: {ref_name: set(keys)}}. References: each checklist, 'consensus'
    (in >=2 checklists, or the only checklist if just one), 'union', and 'gold'.
    """
    refs = {}
    for sp in species_list:
        per = {}
        votes = defaultdict(set)
        for cl in checklists:
            names = authority_synonyms(sp, [cl])
            keys = {_key(n) for n in names}
            per[cl.name] = keys
            for k in keys:
                votes[k].add(cl.name)
        need = 2 if len(checklists) >= 2 else 1
        per["consensus"] = {k for k, v in votes.items() if len(v) >= need}
        per["union"] = set(votes)
        if gold is not None:
            per["gold"] = {_key(n) for n in gold.get(sp, [])}
        refs[sp] = per
    return refs


def known_anywhere(key, checklists, gbif, raw_name):
    if any(cl.by_key.get(key) for cl in checklists):
        return True
    if gbif is not None:
        return bool(gbif.name_exists(raw_name))
    return False


def cmd_score(args):
    checklists = []
    if args.wcvp:
        checklists.append(load_wcvp(args.wcvp))
    if args.wfo:
        checklists.append(load_wfo(args.wfo))
    for item in args.dwc or []:
        n, p = item.split("=", 1)
        checklists.append(load_dwc(p, n))
    if not checklists:
        sys.exit("Need at least one checklist (--wcvp/--wfo/--dwc) to score against.")
    gbif = GBIFBackbone() if args.gbif else None
    gold = None
    if args.gold:
        gold = defaultdict(list)
        with open(args.gold, encoding="utf-8") as f:
            for r in csv.DictReader(f):
                gold[r["species"].strip()].append(r["name"].strip())

    recs = [r for r in LLMLog(args.log).records() if not args.prompt_version
            or r["prompt_version"] == args.prompt_version]
    # last record wins per (model, mode, species, run)
    by = {}
    for r in recs:
        by[(r["provider"], r["model"], r["mode"], r["species"], r["run_index"])] = r
    models = sorted({(p, m, mo) for (p, m, mo, _s, _k) in by})
    species = sorted({s for (_p, _m, _mo, s, _k) in by})
    refs = build_references(species, checklists, gold)
    ref_names = sorted({rn for v in refs.values() for rn in v})
    os.makedirs(args.out, exist_ok=True)

    per_species_rows, per_model_rows = [], []
    recovered = {}   # (model, ref) -> {(species, key): 0/1}  (any run)
    for (p, m, mo) in models:
        label = f"{p}:{m}" + ("+web" if mo == "web" else "")
        halluc_n = halluc_tot = calls = fails = 0
        stab = []
        for ref in ref_names:
            recovered[(label, ref)] = {}
        agg = defaultdict(lambda: {"tp": 0, "pred": 0, "ref": 0})
        for sp in species:
            runs = [by.get((p, m, mo, sp, k)) for k in range(args.runs)]
            runs = [r for r in runs if r is not None]
            calls += len(runs)
            fails += sum(1 for r in runs if r.get("error"))
            ok = [r for r in runs if not r.get("error")]
            run_sets = [{_key(n) for n in synonym_names(r)} - {_key(sp)} for r in ok]
            raw = {}
            for r in ok:
                for n in synonym_names(r):
                    raw.setdefault(_key(n), n)
            if len(run_sets) >= 2:
                stab.append(mean_pairwise_jaccard(run_sets))
            any_run = set().union(*run_sets) if run_sets else set()
            for k in any_run:
                halluc_tot += 1
                if not known_anywhere(k, checklists, gbif, raw[k]):
                    halluc_n += 1
            for ref in ref_names:
                R = refs[sp][ref]
                pr = precision_recall(any_run, R)
                agg[ref]["tp"] += pr["tp"]
                agg[ref]["pred"] += pr["n_pred"]
                agg[ref]["ref"] += pr["n_ref"]
                per_species_rows.append({"model": label, "species": sp, "reference": ref,
                                         "runs_ok": len(ok), **pr})
                for k in R:
                    recovered[(label, ref)][(sp, k)] = int(k in any_run)
        h = wilson_ci(halluc_n, halluc_tot)
        stab_med = bootstrap_ci(stab, n_boot=2000) if stab else (float("nan"),) * 3
        for ref in ref_names:
            a = agg[ref]
            prec = wilson_ci(a["tp"], a["pred"])
            rec = wilson_ci(a["tp"], a["ref"])
            f1 = (2 * prec[0] * rec[0] / (prec[0] + rec[0])) if a["pred"] and a["ref"] and (prec[0] + rec[0]) else float("nan")
            per_model_rows.append({
                "model": label, "reference": ref, "calls": calls, "failed_calls": fails,
                "precision": prec[0], "precision_lo": prec[1], "precision_hi": prec[2],
                "recall": rec[0], "recall_lo": rec[1], "recall_hi": rec[2], "f1": f1,
                "hallucination_rate": h[0], "halluc_lo": h[1], "halluc_hi": h[2],
                "run_jaccard_median": stab_med[0], "run_jaccard_lo": stab_med[1],
                "run_jaccard_hi": stab_med[2],
            })

    # paired tests per reference: Cochran's Q across models, then McNemar + Holm
    labels = sorted({lab for (lab, _r) in recovered})
    pair_rows, q_out = [], {}
    for ref in ref_names:
        items = sorted(set().union(*[recovered[(lab, ref)].keys() for lab in labels]))
        if len(labels) >= 2 and items:
            mat = [[recovered[(lab, ref)].get(it, 0) for lab in labels] for it in items]
            q_out[ref] = {"models": labels, "n_items": len(items), **cochrans_q(mat)}
            tests = []
            for a, b in combinations(labels, 2):
                va = [recovered[(a, ref)].get(it, 0) for it in items]
                vb = [recovered[(b, ref)].get(it, 0) for it in items]
                tests.append({"reference": ref, "model_a": a, "model_b": b, **mcnemar_exact(va, vb)})
            adj = holm([t["p"] for t in tests])
            for t, pa in zip(tests, adj):
                t["p_holm"] = pa
            pair_rows += tests

    def write(name, rows):
        if not rows:
            return
        with open(os.path.join(args.out, name), "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)

    write("per_model.csv", per_model_rows)
    write("per_species.csv", per_species_rows)
    write("pairwise_mcnemar.csv", pair_rows)
    with open(os.path.join(args.out, "cochran_q.json"), "w") as f:
        json.dump(q_out, f, indent=1)
    print(f"[INFO] Scored {len(models)} model configuration(s) on {len(species)} species "
          f"against {len(ref_names)} reference set(s). Output: {args.out}")
    for r in per_model_rows:
        if r["reference"] == "consensus":
            print(f"  {r['model']:<45} P={r['precision']:.2f} R={r['recall']:.2f} "
                  f"halluc={r['hallucination_rate']:.2f} stability={r['run_jaccard_median']:.2f}")


def main(argv=None):
    ap = argparse.ArgumentParser(prog="missmap_benchmark")
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run", help="Query models and log responses")
    r.add_argument("--species-file")
    r.add_argument("--species", help="Comma-separated species (alternative to --species-file)")
    r.add_argument("--models", nargs="+", required=True)
    r.add_argument("--runs", type=int, default=5)
    r.add_argument("--temperature", type=float, default=0.0,
                   help="Used only by OpenAI-compatible providers (Claude models reject it).")
    r.add_argument("--effort", choices=["low", "medium", "high", "xhigh", "max"],
                   help="Claude effort level (recorded in the log).")
    r.add_argument("--log", default="missmap_results/benchmark_calls.jsonl")
    s = sub.add_parser("score", help="Score logged responses against authorities")
    s.add_argument("--log", default="missmap_results/benchmark_calls.jsonl")
    s.add_argument("--wcvp")
    s.add_argument("--wfo")
    s.add_argument("--dwc", action="append")
    s.add_argument("--gbif", action="store_true", help="Also count GBIF exact matches as 'name exists'.")
    s.add_argument("--gold", help="CSV with columns species,name (hand-adjudicated synonyms).")
    s.add_argument("--runs", type=int, default=5)
    s.add_argument("--prompt-version", help="Only score calls made with this prompt version.")
    s.add_argument("--out", default="missmap_results/benchmark")
    args = ap.parse_args(argv)
    if args.cmd == "run":
        if not (args.species_file or args.species):
            ap.error("run needs --species-file or --species")
        cmd_run(args)
    else:
        cmd_score(args)


if __name__ == "__main__":
    main()
