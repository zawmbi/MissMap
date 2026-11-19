"""
Tests for the reconciliation layer: name parsing/matching, authority loaders,
TaxID disambiguation, hidden-record classification, the LLM harness, the
statistics, the benchmark scorer, and the sampler. Everything external (NCBI,
GBIF, LLM APIs) is faked.

Fixture scenario (mirrors the proposal's example; IDs are made up):
  WCVP: Lampranthus elegans is a SYNONYM of L. emarginatoides
  WFO:  both are ACCEPTED (separate species)
  NCBI: L. emarginatoides = txid 1001 (no records), L. elegans = txid 1002 (6 records)
"""

import csv
import importlib.util
import json
import os

import pytest

import missmap_names as mn
import missmap_authorities as ma
import missmap_ncbi as mncbi
import missmap_classify as mc
import missmap_llm as ml
import missmap_stats as ms
import missmap_benchmark as mb
import missmap_sample as msamp


def _load_mm():
    here = os.path.dirname(os.path.abspath(__file__))
    spec = importlib.util.spec_from_file_location("missmap_rec_under_test", os.path.join(here, "MissMap.py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


mm = _load_mm()


# ======================================================================
# Fixtures
# ======================================================================
WCVP_ROWS = [
    # plant_name_id|taxon_name|taxon_authors|taxon_rank|taxon_status|accepted_plant_name_id|family|ipni_id|powo_id
    ("1", "Lampranthus emarginatoides", "(Haw.) N.E.Br.", "Species", "Accepted", "1", "Aizoaceae", "ipni-1", "p1"),
    ("2", "Lampranthus elegans", "(Jacq.) Schwantes", "Species", "Synonym", "1", "Aizoaceae", "ipni-2", "p2"),
    ("3", "Mesembryanthemum emarginatoides", "Haw.", "Species", "Synonym", "1", "Aizoaceae", "ipni-3", "p3"),
    ("4", "Lampranthus aureus", "(L.) N.E.Br.", "Species", "Accepted", "4", "Aizoaceae", "ipni-4", "p4"),
    ("5", "Lampranthus emarginatoides var. minor", "L.Bolus", "Variety", "Synonym", "1", "Aizoaceae", "", ""),
    ("6", "Solanum lycopersicum", "L.", "Species", "Accepted", "6", "Solanaceae", "ipni-6", "p6"),
    ("7", "Lycopersicon esculentum", "Mill.", "Species", "Synonym", "6", "Solanaceae", "ipni-7", "p7"),
    ("8", "Solanum pimpinellifolium", "L.", "Species", "Accepted", "8", "Solanaceae", "ipni-8", "p8"),
]

WFO_ROWS = [
    # taxonID scientificName scientificNameAuthorship taxonRank taxonomicStatus acceptedNameUsageID family scientificNameID
    ("wfo-1", "Lampranthus emarginatoides (Haw.) N.E.Br.", "(Haw.) N.E.Br.", "species", "Accepted", "", "Aizoaceae", "urn:lsid:ipni.org:names:ipni-1"),
    ("wfo-2", "Lampranthus elegans (Jacq.) Schwantes", "(Jacq.) Schwantes", "species", "Accepted", "", "Aizoaceae", "urn:lsid:ipni.org:names:ipni-2"),
    ("wfo-3", "Mesembryanthemum emarginatoides Haw.", "Haw.", "species", "Synonym", "wfo-1", "Aizoaceae", ""),
    ("wfo-6", "Solanum lycopersicum L.", "L.", "species", "Accepted", "", "Solanaceae", ""),
    ("wfo-7", "Lycopersicon esculentum Mill.", "Mill.", "species", "Synonym", "wfo-6", "Solanaceae", ""),
    ("wfo-8", "Solanum pimpinellifolium L.", "L.", "species", "Synonym", "wfo-6", "Solanaceae", ""),
]


@pytest.fixture
def wcvp_path(tmp_path):
    p = tmp_path / "wcvp_names.csv"
    with open(p, "w", encoding="utf-8") as f:
        f.write("plant_name_id|taxon_name|taxon_authors|taxon_rank|taxon_status|"
                "accepted_plant_name_id|family|ipni_id|powo_id\n")
        for r in WCVP_ROWS:
            f.write("|".join(r) + "\n")
    return str(p)


@pytest.fixture
def wfo_path(tmp_path):
    p = tmp_path / "classification.csv"
    with open(p, "w", encoding="utf-8") as f:
        f.write("taxonID\tscientificName\tscientificNameAuthorship\ttaxonRank\ttaxonomicStatus\t"
                "acceptedNameUsageID\tfamily\tscientificNameID\n")
        for r in WFO_ROWS:
            f.write("\t".join(r) + "\n")
    return str(p)


@pytest.fixture
def checklists(wcvp_path, wfo_path):
    return [ma.load_wcvp(wcvp_path), ma.load_wfo(wfo_path)]


def taxon(taxid, name, rank="species", lineage=("Eukaryota", "Viridiplantae", "Streptophyta"),
          synonyms=()):
    return mncbi.NcbiTaxon(taxid=str(taxid), scientific_name=name, rank=rank,
                           lineage=list(lineage), synonyms=list(synonyms))


# ======================================================================
# Names
# ======================================================================
def test_parse_binomial_with_author():
    p = mn.parse_name("Solanum lycopersicum L.")
    assert p.parsed and p.canonical == "Solanum lycopersicum" and p.authorship == "L."


def test_parse_infraspecific_and_hybrid():
    p = mn.parse_name("Zea mays subsp. mexicana (Schrad.) Iltis")
    assert p.infra_rank == "subsp." and p.infra_epithet == "mexicana" and p.is_infraspecific
    assert p.canonical == "Zea mays subsp. mexicana"
    assert mn.parse_name("Zea mays ssp. mexicana").key() == p.key()
    h = mn.parse_name("Mentha × piperita L.")
    assert h.hybrid and h.canonical == "Mentha × piperita"


def test_common_name_is_not_scientific():
    assert not mn.looks_like_scientific_name("tomato")
    assert not mn.looks_like_scientific_name("Solanum sp.")
    assert mn.looks_like_scientific_name("Lycopersicon esculentum Mill.")


def test_match_kinds():
    assert mn.match_names("Lampranthus elegans", "Lampranthus elegans (Jacq.) Schwantes").kind == "exact"
    assert mn.match_names("Rosa canina", "Rosa caninus").kind == "orthographic_variant"
    assert mn.match_names("Rosa canina", "Rosa caniina").kind == "fuzzy_candidate"
    assert not mn.match_names("Rosa canina", "Rosa caniina").accepted
    # different genus never matches by string similarity
    assert mn.match_names("Lycopersicon esculentum", "Solanum esculentum").kind == "no_match"


def test_authorship_relation():
    assert mn.authorship_relation("(L.) Mill.", "(L.)Mill.") == "same"
    assert mn.authorship_relation("(Haw.) N.E.Br.", "N.E.Br.") == "basionym_differs"
    assert mn.authorship_relation("L.", "Mill.") == "different"
    assert mn.authorship_relation("", "Mill.") == "unknown"


# ======================================================================
# Authorities
# ======================================================================
def test_load_wcvp_and_resolve(checklists):
    wcvp, wfo = checklists
    assert len(wcvp) == len(WCVP_ROWS)
    acc = wcvp.accepted_for("Lampranthus elegans")
    assert [a.scientificName for a in acc] == ["Lampranthus emarginatoides"]
    assert wcvp.by_id["2"].externalIDs["ipni_id"] == "ipni-2"
    # WFO splits authorship off scientificName and keeps the IPNI id
    assert wfo.by_id["wfo-2"].scientificName == "Lampranthus elegans"
    assert wfo.by_id["wfo-2"].externalIDs["ipni"] == "ipni-2"
    assert wfo.status_of("Lampranthus elegans") == "accepted"


def test_compare_synonymy_lampranthus(checklists):
    cmp_ = ma.compare_synonymy("Lampranthus emarginatoides", checklists, ncbi_synonyms=[])
    rows = {r["name"]: r for r in cmp_["names"]}
    eleg = rows["Lampranthus elegans"]
    assert eleg["WCVP"] == "synonym" and eleg["WFO"] == "accepted_other" and eleg["NCBI"] == "not_synonym"
    assert eleg["class"] == "lump_split"
    assert cmp_["jaccard"]["WCVP|WFO"] == pytest.approx(1 / 2)


def test_tomato_split_and_jaccard(checklists):
    cmp_ = ma.compare_synonymy("Solanum lycopersicum", checklists,
                               ncbi_synonyms=["Lycopersicon esculentum"])
    rows = {r["name"]: r for r in cmp_["names"]}
    # WFO lumps S. pimpinellifolium into tomato; WCVP keeps it separate
    assert rows["Solanum pimpinellifolium"]["class"] == "lump_split"
    assert rows["Lycopersicon esculentum"]["class"] == "agreed"


def test_verify_crosswalk_with_fake_wikidata():
    fake = lambda name: [{"item": "Q1", "ncbi": ["4081"], "ipni": ["ipni-6"]}]
    assert ma.verify_crosswalk("Solanum lycopersicum", ncbi_taxid="4081", lookup=fake)["status"] == "consistent"
    assert ma.verify_crosswalk("Solanum lycopersicum", ncbi_taxid="999", lookup=fake)["status"] == "conflict"
    two = lambda name: [{"item": "Q1"}, {"item": "Q2"}]
    assert ma.verify_crosswalk("Prunella", lookup=two)["status"] == "unverified"


# ======================================================================
# NCBI taxonomy: XML, disambiguation, taxdump, cache
# ======================================================================
TAXON_XML = """<?xml version="1.0" ?>
<TaxaSet>
<Taxon><TaxId>4081</TaxId><ScientificName>Solanum lycopersicum</ScientificName>
<OtherNames><GenbankCommonName>tomato</GenbankCommonName><Synonym>Lycopersicon esculentum</Synonym>
<Authority>Solanum lycopersicum L.</Authority><EquivalentName>Lycopersicum esculentum</EquivalentName></OtherNames>
<Rank>species</Rank><Division>Plants and Fungi</Division>
<Lineage>cellular organisms; Eukaryota; Viridiplantae; Streptophyta; Solanaceae; Solanum</Lineage></Taxon>
</TaxaSet>"""


def test_parse_taxon_xml_separates_common_names():
    t = mncbi.parse_taxon_xml(TAXON_XML)[0]
    assert t.taxid == "4081" and t.rank == "species"
    assert "Lycopersicon esculentum" in t.synonyms and "Lycopersicum esculentum" in t.synonyms
    assert "tomato" not in t.synonyms and t.common_names == ["tomato"]
    assert t.in_group("Viridiplantae")


def test_choose_taxon_homonym_resolved_by_group():
    plant = taxon(1, "Prunella vulgaris")
    bird = taxon(2, "Prunella vulgaris", lineage=("Eukaryota", "Metazoa", "Aves"))
    chosen, flags = mncbi.choose_taxon("Prunella vulgaris", [bird, plant], expected_group="Viridiplantae")
    assert chosen.taxid == "1" and "homonym_other_group" in flags
    chosen, flags = mncbi.choose_taxon("Prunella vulgaris", [bird, plant])
    assert chosen is None and "ambiguous" in flags


def test_choose_taxon_prefers_exact_scientific_name_and_rank():
    sp = taxon(10, "Zea mays")
    ssp = taxon(11, "Zea mays subsp. mays", rank="subspecies")
    other = taxon(12, "Zea perennis", synonyms=["Zea mays"])
    chosen, flags = mncbi.choose_taxon("Zea mays", [other, ssp, sp])
    assert chosen.taxid == "10" and flags == []
    chosen, flags = mncbi.choose_taxon("Euchlaena perennis", [other])
    assert chosen.taxid == "12" and "matched_via_synonym" in flags


def test_taxdump_reader(tmp_path):
    names = [("1", "root", "scientific name"), ("2759", "Eukaryota", "scientific name"),
             ("33090", "Viridiplantae", "scientific name"), ("4081", "Solanum lycopersicum", "scientific name"),
             ("4081", "Lycopersicon esculentum", "synonym"), ("4081", "tomato", "genbank common name"),
             ("99", "Solanum lycopersicum var. cerasiforme", "scientific name")]
    nodes = [("1", "1", "no rank"), ("2759", "1", "superkingdom"), ("33090", "2759", "kingdom"),
             ("4081", "33090", "species"), ("99", "4081", "varietas")]
    with open(tmp_path / "names.dmp", "w") as f:
        for t, n, c in names:
            f.write(f"{t}\t|\t{n}\t|\t\t|\t{c}\t|\n")
    with open(tmp_path / "nodes.dmp", "w") as f:
        for t, p, r in nodes:
            f.write(f"{t}\t|\t{p}\t|\t{r}\t|\t\t|\n")
    td = mncbi.NcbiTaxdump(str(tmp_path))
    c = td.candidates("Lycopersicon esculentum")
    assert [x.taxid for x in c] == ["4081"]
    assert c[0].lineage == ["Eukaryota", "Viridiplantae"] and c[0].common_names == ["tomato"]
    assert td.children("4081") == ["99"]


def test_entrez_cache_roundtrip_and_expiry(tmp_path, monkeypatch):
    cache = mncbi.EntrezCache(str(tmp_path / "c.sqlite"), max_age_hours=1)
    cache.set("nuccore_count", "q", 5)
    assert cache.get("nuccore_count", "q") == 5
    now = mncbi.time.time()
    monkeypatch.setattr(mncbi.time, "time", lambda: now + 7200)
    assert cache.get("nuccore_count", "q") is None


# ======================================================================
# Classification + queries
# ======================================================================
def test_classify_buckets(checklists):
    sp, primary = "Lampranthus emarginatoides", "1001"
    c = mc.classify_name("Lampranthus elegans", sp, primary, [taxon(1002, "Lampranthus elegans")], checklists)
    assert c["bucket"] == "other_taxid_confirmed" and c["hidden_taxids"] == ["1002"]
    c = mc.classify_name("Lampranthus aureus", sp, primary, [taxon(1003, "Lampranthus aureus")], checklists)
    assert c["bucket"] == "other_taxid_wrong_species"
    c = mc.classify_name("Mesembryanthemum emarginatoides", sp, primary,
                         [taxon(1001, "Lampranthus emarginatoides")], checklists)
    assert c["bucket"] == "same_taxid"
    c = mc.classify_name("Lampranthus inventus", sp, primary, [], checklists, name_exists=lambda n: False)
    assert c["bucket"] == "not_found"
    c = mc.classify_name("Lampranthus inventus", sp, primary, [], checklists, name_exists=lambda n: True)
    assert c["bucket"] == "valid_not_in_ncbi"
    assert mc.classify_name("ice plant", sp, primary, [], checklists)["bucket"] == "not_a_scientific_name"
    # a child taxon (variety) of the primary TaxID is already covered by [Organism:exp]
    c = mc.classify_name("Lampranthus emarginatoides var. minor", sp, primary,
                         [taxon(1009, "Lampranthus emarginatoides var. minor", rank="varietas")],
                         checklists, descendant_taxids={"1009"})
    assert c["bucket"] == "same_taxid"


def test_hidden_and_overcount_queries():
    q = mc.hidden_records_query(["1002", "1003"], "1001", '"rbcL"[Gene]')
    assert q == '((txid1002[Organism:exp] OR txid1003[Organism:exp]) NOT txid1001[Organism:exp]) AND ("rbcL"[Gene])'
    assert mc.hidden_records_query([], "1001") == ""
    assert mc.overcount_query("Solanum lycopersicum", "4081") == \
        '"Solanum lycopersicum"[All Fields] NOT txid4081[Organism:exp]'


# ======================================================================
# MissMap end-to-end (Lampranthus scenario), NCBI and LLM faked
# ======================================================================
NCBI_TAXA = {
    "lampranthus emarginatoides": [taxon(1001, "Lampranthus emarginatoides",
                                         synonyms=["Mesembryanthemum emarginatoides"])],
    "mesembryanthemum emarginatoides": [taxon(1001, "Lampranthus emarginatoides")],
    "lampranthus elegans": [taxon(1002, "Lampranthus elegans")],
    "lampranthus aureus": [taxon(1003, "Lampranthus aureus")],
}


def fake_count(term, context=""):
    if "NOT txid1001" in term:                      # hidden records (under L. elegans)
        if '"rbcL"' in term:
            return 2
        if "[All Fields]" in term or '"matK"' in term:
            return 0
        return 6
    if "[All Fields]" in term:
        return 3
    if "txid1001" in term:
        return 0                                     # nothing filed under the accepted name
    return 0


@pytest.fixture
def mm_program(monkeypatch, tmp_path, checklists):
    monkeypatch.setattr(mm, "RESULTS_DIR", str(tmp_path))
    monkeypatch.setattr(mm.time, "sleep", lambda *a, **k: None)
    p = mm.MissMap(checklists=checklists, group="Viridiplantae")
    p.delay = 0
    p.species_list = ["Lampranthus emarginatoides"]
    p.dt_list = ["rbcL", "matK"]
    monkeypatch.setattr(p, "_ncbi_candidates",
                        lambda name, filter_group=True: NCBI_TAXA.get(name.strip().lower(), []))
    monkeypatch.setattr(p, "_descendant_taxa", lambda taxid: [])
    monkeypatch.setattr(p, "lookup_synonyms_with_ai",
                        lambda sp, manual: ["Lampranthus elegans", "Lampranthus aureus",
                                            "Lampranthus fabricatus", "ice plant"])
    monkeypatch.setattr(p, "_esearch_count", fake_count)
    return p


def test_end_to_end_hidden_records(mm_program, tmp_path):
    p = mm_program
    p.data_counts()
    row = p.rows[0]
    assert row["TaxID"] == "1001"
    assert row["rbcL"] == "0"
    assert row["Hidden TaxIDs (authority synonyms)"] == "1002"
    assert row["Hidden Records (authority synonyms)"] == "6"
    assert row["Hidden Records (all)"] == "6"
    assert row["Markers Gained via Synonyms"] == "rbcL(+2)"
    assert row["AI: other_taxid_confirmed"] == "Lampranthus elegans"
    assert row["AI: other_taxid_wrong_species"] == "Lampranthus aureus"
    assert row["AI: not_found"] == "Lampranthus fabricatus"
    assert row["AI: not_a_scientific_name"] == "ice plant"
    assert "Lampranthus elegans" in row["AI-validated Synonyms"]
    assert "Lampranthus aureus" in row["AI Rejected Synonyms"]
    assert row["WFO Accepted Name"] == "Lampranthus emarginatoides"
    assert row["Hidden Share of Reconciled Total (%)"] == "100"
    assert p._matrix_records[0]["counts_reconciled"]["rbcL"] == 2

    p.display_and_save_matrix()
    files = os.listdir(tmp_path)
    for prefix in ("missmap_presence_reconciled_", "name_classification_",
                   "synonymy_comparison_", "crosswalk_dwc_", "run_manifest_"):
        assert any(f.startswith(prefix) for f in files), prefix
    man = json.load(open(tmp_path / next(f for f in files if f.startswith("run_manifest_"))))
    assert set(man["authorities"]) == {"WCVP", "WFO"}
    assert man["expected_lineage_group"] == "Viridiplantae"


def test_homonym_species_gets_no_taxid(monkeypatch, tmp_path):
    monkeypatch.setattr(mm, "RESULTS_DIR", str(tmp_path))
    p = mm.MissMap()
    bird = taxon(2, "Prunella modularis", lineage=("Eukaryota", "Metazoa"))
    plant = taxon(1, "Prunella modularis")
    monkeypatch.setattr(p, "_ncbi_candidates", lambda name, filter_group=True: [bird, plant])
    syns, taxid = p.lookup_synonyms_manually("Prunella modularis")
    assert taxid is None and "ambiguous" in p._taxon_info["Prunella modularis"]["flags"]
    p.group = "Metazoa"
    syns, taxid = p.lookup_synonyms_manually("Prunella modularis")
    assert taxid == "2"


def test_esearch_count_uses_cache(monkeypatch, tmp_path):
    p = mm.MissMap(cache=mncbi.EntrezCache(str(tmp_path / "c.sqlite")))
    p.delay = 0
    calls = {"n": 0}

    class H:
        def __enter__(self): return self
        def __exit__(self, *a): return False

    def read(h):
        calls["n"] += 1
        return {"Count": "9"}
    monkeypatch.setattr(mm.Entrez, "esearch", lambda **k: H())
    monkeypatch.setattr(mm.Entrez, "read", read)
    assert p._esearch_count("t") == 9 and p._esearch_count("t") == 9
    assert calls["n"] == 1


def test_cli_flags_parse(tmp_path):
    cfg = mm.Configuration(str(tmp_path / "none.txt"))
    cfg.terminal_inputs(["-s", "Zea mays", "--wcvp", "w.csv", "--wfo", "f.csv", "--dwc", "MDD=m.tsv",
                         "--group", "Viridiplantae", "--llm", "anthropic:claude-sonnet-5-5",
                         "--no-cache", "--non-interactive"])
    assert cfg.wcvp_path == "w.csv" and cfg.dwc_paths == {"MDD": "m.tsv"}
    assert cfg.group == "Viridiplantae" and cfg.llm_spec == "anthropic:claude-sonnet-5-5"
    assert cfg.use_cache is False


# ======================================================================
# LLM harness
# ======================================================================
def test_extract_and_normalise_json():
    txt = 'Here you go:\n```json\n{"accepted_name": "X y", "synonyms": ["A b"], "common_names": []}\n```'
    out = ml.normalise_result(ml.extract_json(txt))
    assert out["synonyms"][0] == {"name": "A b", "authorship": "", "type": "unknown", "source": ""}
    with pytest.raises(ValueError):
        ml.extract_json("no json here")


class FakeProvider:
    name, model, mode = "fake", "fake-1", "closed"
    settings = {}

    def __init__(self, replies):
        self.replies = list(replies)
        self.calls = 0

    def complete(self, species, known=None):
        self.calls += 1
        r = self.replies.pop(0)
        if isinstance(r, Exception):
            raise r
        return r, {"input_tokens": 1, "output_tokens": 1}


def test_query_synonyms_logs_and_caches(tmp_path):
    log = ml.LLMLog(str(tmp_path / "calls.jsonl"))
    good = json.dumps({"accepted_name": "S l", "synonyms": [{"name": "Lycopersicon esculentum",
                       "authorship": "Mill.", "type": "homotypic", "source": "WCVP"}], "common_names": ["tomato"]})
    prov = FakeProvider([RuntimeError("HTTP 500"), good])
    r1 = ml.query_synonyms(prov, "Solanum lycopersicum", run_index=0, log=log)
    assert r1["error"].startswith("RuntimeError") and ml.synonym_names(r1) == []
    r2 = ml.query_synonyms(prov, "Solanum lycopersicum", run_index=0, log=log)
    assert ml.synonym_names(r2) == ["Lycopersicon esculentum"]
    # third call is served from the log; the provider is not called again
    ml.query_synonyms(prov, "Solanum lycopersicum", run_index=0, log=ml.LLMLog(str(tmp_path / "calls.jsonl")))
    assert prov.calls == 2
    assert len(log.records()) == 2


def test_make_provider_specs(monkeypatch):
    monkeypatch.setenv("GROQ_API_KEY", "x")
    g = ml.make_provider("groq:llama-3.3-70b-versatile")
    assert isinstance(g, ml.OpenAICompatProvider) and g.base_url.startswith("https://api.groq.com")
    a = ml.make_provider("anthropic:claude-sonnet-5-5+web", effort="high")
    assert isinstance(a, ml.AnthropicProvider) and a.web_search and a.mode == "web"
    with pytest.raises(ValueError):
        ml.make_provider("claude-sonnet-5-5")


class _Obj:
    def __init__(self, **kw):
        self.__dict__.update(kw)


class FakeAnthropicClient:
    def __init__(self, responses):
        self.responses = list(responses)
        self.sent = []
        self.messages = self

    def create(self, **kw):
        self.sent.append(kw)
        return self.responses.pop(0)


def _resp(text, stop="end_turn"):
    return _Obj(content=[_Obj(type="text", text=text)], stop_reason=stop,
                usage=_Obj(input_tokens=10, output_tokens=5, server_tool_use=None), stop_details=None)


def test_anthropic_provider_request_shape():
    body = json.dumps({"accepted_name": "A b", "synonyms": [], "common_names": []})
    client = FakeAnthropicClient([_resp(body)])
    prov = ml.AnthropicProvider("claude-haiku-5-5", effort="low", client=client)
    text, usage = prov.complete("A b")
    sent = client.sent[0]
    assert sent["model"] == "claude-haiku-5-5"
    assert sent["output_config"]["format"]["type"] == "json_schema"
    assert sent["output_config"]["effort"] == "low"
    assert "temperature" not in sent and "tools" not in sent
    assert json.loads(text)["accepted_name"] == "A b"


def test_anthropic_provider_web_mode_and_refusal():
    body = '{"accepted_name": "A b", "synonyms": [], "common_names": []}'
    client = FakeAnthropicClient([_resp("searching", "pause_turn"), _resp(body)])
    prov = ml.AnthropicProvider("claude-sonnet-5-5", web_search=True, client=client)
    text, _ = prov.complete("A b")
    assert client.sent[0]["tools"][0]["type"] == "web_search_20260209"
    assert "format" not in client.sent[0].get("output_config", {})
    assert len(client.sent) == 2 and json.loads(text)["accepted_name"] == "A b"
    refuse = _resp("", "refusal")
    refuse.stop_details = _Obj(category="bio")
    prov2 = ml.AnthropicProvider("claude-sonnet-5-5", client=FakeAnthropicClient([refuse]))
    with pytest.raises(RuntimeError, match="refusal"):
        prov2.complete("A b")


# ======================================================================
# Statistics
# ======================================================================
def test_wilson_ci():
    p, lo, hi = ms.wilson_ci(0, 10)
    assert p == 0 and lo == 0 and hi == pytest.approx(0.2775, abs=1e-3)
    p, lo, hi = ms.wilson_ci(50, 100)
    assert lo == pytest.approx(0.4038, abs=1e-3) and hi == pytest.approx(0.5962, abs=1e-3)


def test_mcnemar_and_holm():
    a = [1] * 10 + [0] * 10
    b = [1] * 10 + [1] * 8 + [0] * 2
    r = ms.mcnemar_exact(a, b)
    assert (r["b01"], r["b10"]) == (8, 0) and r["p"] == pytest.approx(2 / 256)
    assert ms.holm([0.01, 0.04, 0.03]) == pytest.approx([0.03, 0.06, 0.06])


def test_cochrans_q_known_value():
    # Textbook example: 3 treatments; Q = 2.0 for this table
    mat = [[1, 1, 0], [1, 0, 0], [1, 1, 1], [0, 1, 0], [1, 0, 0], [1, 1, 0]]
    q = ms.cochrans_q(mat)
    col = [5, 4, 1]
    k, N = 3, 10
    row = [2, 1, 3, 1, 1, 2]
    expected = (k - 1) * (k * sum(c * c for c in col) - N * N) / (k * N - sum(r * r for r in row))
    assert q["Q"] == pytest.approx(expected) and q["df"] == 2 and 0 < q["p"] < 1


def test_fleiss_kappa_perfect_and_chance():
    assert ms.fleiss_kappa([[3, 0], [0, 3], [3, 0]]) == pytest.approx(1.0)
    assert ms.fleiss_kappa([[2, 1], [1, 2], [2, 1], [1, 2]]) < 0


def test_jaccard_and_bootstrap():
    assert ms.jaccard({1, 2}, {2, 3}) == pytest.approx(1 / 3)
    assert ms.mean_pairwise_jaccard([{1}, {1}, {1}]) == 1.0
    med, lo, hi = ms.bootstrap_ci([1, 2, 3, 4, 100], n_boot=500)
    assert med == 3 and lo <= med <= hi


# ======================================================================
# Benchmark scoring + sampler
# ======================================================================
def test_benchmark_score(tmp_path, wcvp_path, wfo_path):
    log = ml.LLMLog(str(tmp_path / "calls.jsonl"))
    replies = {
        "good": ["Lycopersicon esculentum"],
        "bad": ["Lycopersicon esculentum", "Solanum inventum"],
    }
    for model, names in replies.items():
        prov = FakeProvider([json.dumps({"accepted_name": "", "synonyms": names, "common_names": []})] * 2)
        prov.model = model
        for k in range(2):
            ml.query_synonyms(prov, "Solanum lycopersicum", run_index=k, log=log)
    out = tmp_path / "bench"
    mb.main(["score", "--log", str(tmp_path / "calls.jsonl"), "--wcvp", wcvp_path,
             "--wfo", wfo_path, "--runs", "2", "--out", str(out)])
    rows = list(csv.DictReader(open(out / "per_model.csv")))
    cons = {r["model"]: r for r in rows if r["reference"] == "consensus"}
    assert float(cons["fake:good"]["recall"]) == 1.0 and float(cons["fake:good"]["precision"]) == 1.0
    assert float(cons["fake:bad"]["precision"]) == 0.5
    assert float(cons["fake:bad"]["hallucination_rate"]) == 0.5
    assert float(cons["fake:good"]["run_jaccard_median"]) == 1.0
    assert (out / "pairwise_mcnemar.csv").exists() and (out / "cochran_q.json").exists()


def test_stratified_sample(wcvp_path):
    cl = ma.load_wcvp(wcvp_path)
    rows = msamp.stratified_sample(cl, n=4, seed=1)
    assert rows and all(r["design_weight"] >= 1 for r in rows)
    assert {r["stratum"] for r in rows} <= {"0", "1-2", "3-9", "10+"}
    assert rows == msamp.stratified_sample(cl, n=4, seed=1)   # reproducible


def test_trinomial_without_rank_marker():
    p = mn.parse_name("Panthera tigris altaica Temminck, 1844")
    assert p.infra_epithet == "altaica" and p.infra_rank == "" and p.authorship.startswith("Temminck")
    assert mn.parse_name("Lycopersicon esculentum cerasiforme").is_infraspecific
    # an author particle is not an epithet
    assert not mn.parse_name("Aster amellus de Wild.").is_infraspecific


def test_gbif_canonical_restores_rank_marker():
    rec = {"canonicalName": "Lycopersicon esculentum cerasiforme", "rank": "VARIETY"}
    assert ma.GBIFBackbone._canonical_with_rank(rec) == "Lycopersicon esculentum var. cerasiforme"
    assert ma.GBIFBackbone._canonical_with_rank({"canonicalName": "Solanum pomiferum", "rank": "SPECIES"}) \
        == "Solanum pomiferum"


class FakeGBIF:
    def synonym_names(self, name):
        return name, ["Lampranthus novus"]

    def name_exists(self, name):
        return False


def test_gbif_synonym_confirms_other_taxid(mm_program, monkeypatch):
    p = mm_program
    p.gbif = FakeGBIF()
    taxa = dict(NCBI_TAXA)
    taxa["lampranthus novus"] = [taxon(1005, "Lampranthus novus")]
    monkeypatch.setattr(p, "_ncbi_candidates", lambda name, filter_group=True: taxa.get(name.strip().lower(), []))
    p.data_counts()
    row = p.rows[0]
    assert "1005" in row["Hidden TaxIDs (authority synonyms)"].split(",")
    nr = [r for r in p._name_rows if r["Name"] == "Lampranthus novus"][0]
    assert nr["Bucket"] == "other_taxid_confirmed" and nr["Authority Sources"] == "GBIF"
