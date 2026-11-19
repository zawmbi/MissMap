"""
Stratified random sample of accepted species from a checklist (WCVP by default).

Strata = synonym-count bin x (optionally) family size class, with a per-family
cap so very large families don't dominate. Output includes a design weight
(stratum population / stratum sample) so proportions can be weighted back to
the whole flora.

  python missmap_sample.py --wcvp wcvp_names.csv --n 400 --seed 2026 \
      --family-cap 8 --out sample.csv

A second stratifier, NCBI "fame" (record-count quartile), needs NCBI queries;
run MissMap on a larger first-stage sample and stratify on its
'Total NCBI (no filters)' column if you want it.
"""

import argparse
import csv
import random
from collections import defaultdict

from missmap_authorities import load_wcvp, load_dwc

SYN_BINS = [(0, 0, "0"), (1, 2, "1-2"), (3, 9, "3-9"), (10, 10 ** 9, "10+")]


def syn_bin(n):
    for lo, hi, lab in SYN_BINS:
        if lo <= n <= hi:
            return lab
    return "10+"


def stratified_sample(checklist, n, seed=2026, family_cap=None, exclude_hybrids=True):
    rng = random.Random(seed)
    pop = defaultdict(list)
    for u in checklist.by_id.values():
        if not u.is_accepted or u.taxonRank != "species":
            continue
        if exclude_hybrids and "×" in u.scientificName:
            continue
        nsyn = len([s for s in checklist.synonymy(u.taxonID) if s.taxonRank in ("species", "")])
        pop[syn_bin(nsyn)].append((u, nsyn))
    strata = [lab for _lo, _hi, lab in SYN_BINS if pop.get(lab)]
    if not strata:
        return []
    per = max(1, n // len(strata))   # equal allocation, then weight back
    out = []
    for lab in strata:
        items = pop[lab][:]
        rng.shuffle(items)
        fam_count = defaultdict(int)
        chosen = []
        for u, nsyn in items:
            if family_cap and fam_count[u.family] >= family_cap:
                continue
            chosen.append((u, nsyn))
            fam_count[u.family] += 1
            if len(chosen) >= per:
                break
        w = len(pop[lab]) / len(chosen) if chosen else 0
        for u, nsyn in chosen:
            out.append({"species": u.scientificName, "authorship": u.scientificNameAuthorship,
                        "family": u.family, "taxonID": u.taxonID, "n_synonyms": nsyn,
                        "stratum": lab, "stratum_population": len(pop[lab]),
                        "stratum_sample": len(chosen), "design_weight": round(w, 4)})
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(prog="missmap_sample")
    ap.add_argument("--wcvp")
    ap.add_argument("--dwc", help="NAME=PATH of another Darwin Core checklist")
    ap.add_argument("--n", type=int, default=400)
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--family-cap", type=int, default=8)
    ap.add_argument("--out", default="sample.csv")
    a = ap.parse_args(argv)
    if a.wcvp:
        cl = load_wcvp(a.wcvp)
    elif a.dwc:
        name, path = a.dwc.split("=", 1)
        cl = load_dwc(path, name)
    else:
        ap.error("need --wcvp or --dwc")
    rows = stratified_sample(cl, a.n, a.seed, a.family_cap)
    with open(a.out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"Wrote {len(rows)} species to {a.out} (seed {a.seed}; checklist {cl.version}).")


if __name__ == "__main__":
    main()
