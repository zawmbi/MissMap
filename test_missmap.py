"""
Unit + integration tests for MissMap, with NCBI Entrez fully mocked.

These tests cover the correctness-critical behaviour added on top of the
taxid-search base:
  * NA-vs-0: a failed query is recorded as NA, never as a genuine 0;
  * retry/backoff on transient HTTP/network errors;
  * deduplicated union count vs the (overlapping) sum;
  * the presence/absence + count heatmaps, coverage summary, and manifest;
  * the improved infraspecific-name filter;
  * non-interactive guarding of input() prompts.

Run:  pytest -q
"""

import importlib.util
import os

import pytest


# ----------------------------------------------------------------------
# Load MissMap.py as a module (its __main__ guard keeps main() from running).
# ----------------------------------------------------------------------
def _load_mm():
    here = os.path.dirname(os.path.abspath(__file__))
    spec = importlib.util.spec_from_file_location("missmap_under_test",
                                                  os.path.join(here, "MissMap.py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


mm = _load_mm()


class _FakeHandle:
    """Stand-in for the Entrez.esearch context manager handle."""
    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


def _program(**kwargs):
    p = mm.MissMap(**kwargs)
    p.delay = 0  # never actually sleep in tests
    return p


@pytest.fixture(autouse=True)
def _no_sleep(monkeypatch):
    """Neutralise all sleeps so retry/backoff tests are instant."""
    monkeypatch.setattr(mm.time, "sleep", lambda *_a, **_k: None)


# ======================================================================
# _esearch_count: the single seam through which every count flows
# ======================================================================
def test_esearch_count_success(monkeypatch):
    p = _program()
    monkeypatch.setattr(mm.Entrez, "esearch", lambda **k: _FakeHandle())
    monkeypatch.setattr(mm.Entrez, "read", lambda h: {"Count": "42"})
    assert p._esearch_count("term") == 42


def test_genuine_zero_is_zero_not_na(monkeypatch):
    p = _program()
    monkeypatch.setattr(mm.Entrez, "esearch", lambda **k: _FakeHandle())
    monkeypatch.setattr(mm.Entrez, "read", lambda h: {"Count": "0"})
    result = p._esearch_count("term")
    assert result == 0
    assert result is not mm.COUNT_NA
    assert p._failed_cells == []


def test_persistent_http_error_records_na(monkeypatch):
    p = _program()

    def boom(h):
        raise mm.HTTPError("http://ncbi", 500, "server error", {}, None)

    monkeypatch.setattr(mm.Entrez, "esearch", lambda **k: _FakeHandle())
    monkeypatch.setattr(mm.Entrez, "read", boom)

    result = p._esearch_count("term", context="Sp | rbcL")
    assert result is mm.COUNT_NA          # NA, never 0
    assert p._failed_cells == ["Sp | rbcL"]


def test_retries_then_succeeds(monkeypatch):
    p = _program()
    calls = {"n": 0}

    def flaky(h):
        calls["n"] += 1
        if calls["n"] < 3:
            raise mm.HTTPError("http://ncbi", 503, "busy", {}, None)
        return {"Count": "7"}

    monkeypatch.setattr(mm.Entrez, "esearch", lambda **k: _FakeHandle())
    monkeypatch.setattr(mm.Entrez, "read", flaky)

    assert p._esearch_count("term") == 7
    assert calls["n"] == 3               # retried twice, then succeeded
    assert p._failed_cells == []


def test_non_retryable_fails_fast(monkeypatch):
    p = _program()
    calls = {"n": 0}

    def bad(h):
        calls["n"] += 1
        raise ValueError("malformed query syntax")  # not a transient marker

    monkeypatch.setattr(mm.Entrez, "esearch", lambda **k: _FakeHandle())
    monkeypatch.setattr(mm.Entrez, "read", bad)

    assert p._esearch_count("term") is mm.COUNT_NA
    assert calls["n"] == 1               # no retries for a non-transient error


def test_missing_count_field_is_na(monkeypatch):
    p = _program()
    monkeypatch.setattr(mm.Entrez, "esearch", lambda **k: _FakeHandle())
    monkeypatch.setattr(mm.Entrez, "read", lambda h: {"IdList": []})  # no "Count"
    assert p._esearch_count("term") is mm.COUNT_NA


# ======================================================================
# Formatting + infraspecific filter
# ======================================================================
def test_fmt():
    assert mm.MissMap._fmt(0) == "0"
    assert mm.MissMap._fmt(1234) == "1234"
    assert mm.MissMap._fmt(mm.COUNT_NA) == "NA"


def test_infraspecific_filter_catches_many_ranks():
    p = _program(exclude_infraspecific=True)
    names = [
        "Zea mays",
        "Zea mays var. japonica",
        "Zea mays subsp. mexicana",
        "Zea mays ssp. parviglumis",
        "Zea mays f. ceratina",
        "Zea perennis × Zea mays",
        "maize",
    ]
    out = p._filter_infraspecific(names, "Zea mays")
    assert "Zea mays" in out          # canonical always kept
    assert "maize" in out             # common name kept
    for dropped in [
        "Zea mays var. japonica",
        "Zea mays subsp. mexicana",
        "Zea mays ssp. parviglumis",
        "Zea mays f. ceratina",
        "Zea perennis × Zea mays",
    ]:
        assert dropped not in out


def test_infraspecific_filter_noop_when_disabled():
    p = _program(exclude_infraspecific=False)
    names = ["Zea mays", "Zea mays var. japonica"]
    assert p._filter_infraspecific(names, "Zea mays") == names


# ======================================================================
# data_counts: union vs sum, NA propagation
# ======================================================================
def _fake_counts(**overrides):
    """Build a term->count resolver for data_counts integration tests."""
    def resolver(term, context=""):
        if ") OR (" in term:                    # the deduplicated union query
            return overrides.get("union", 100)
        if "organelle:chloroplast" in term:
            return overrides.get("chloroplast", 30)
        if "rbcL" in term:
            return overrides.get("rbcL", 50)
        if "matK" in term:
            return overrides.get("matK", 40)
        return overrides.get("total", 200)      # no-filter totals
    return resolver


def _wire_species(p, monkeypatch, resolver):
    p.species_list = ["Test species"]
    p.dt_list = ["chloroplast", "rbcL", "matK"]
    monkeypatch.setattr(p, "lookup_synonyms_manually", lambda sp: ([], "999"))
    monkeypatch.setattr(p, "lookup_synonyms_with_ai", lambda sp, manual: [])
    monkeypatch.setattr(p, "_esearch_count", resolver)


def test_union_and_sum_are_distinct(monkeypatch):
    p = _program()
    _wire_species(p, monkeypatch, _fake_counts())
    p.data_counts()
    row = p.rows[0]
    # union is a real record count from NCBI's dedup; sum double-counts
    assert row["Records w/ >=1 Selected Marker (deduplicated)"] == "100"
    assert row["Sum of Marker Counts (overlapping, not a record count)"] == "120"
    assert row["Cells Failed (NA)"] == "0"
    # sum (120) exceeding union (100) is exactly the double-counting the fix exposes
    assert int(row["Sum of Marker Counts (overlapping, not a record count)"]) > \
        int(row["Records w/ >=1 Selected Marker (deduplicated)"])


def test_na_cell_is_not_counted_as_zero(monkeypatch):
    p = _program()
    _wire_species(p, monkeypatch, _fake_counts(matK=mm.COUNT_NA))
    p.data_counts()
    row = p.rows[0]
    assert row["matK"] == "NA"                    # not "0"
    assert row["chloroplast"] == "30"
    assert row["rbcL"] == "50"
    assert row["Cells Failed (NA)"] == "1"
    # the NA cell is excluded from the sum rather than added as 0
    assert row["Sum of Marker Counts (overlapping, not a record count)"] == "80"


# ======================================================================
# Outputs: heatmap, coverage summary, manifest, Excel
# ======================================================================
def _matrix_program(monkeypatch, tmp_path):
    monkeypatch.setattr(mm, "RESULTS_DIR", str(tmp_path))
    p = mm.MissMap()
    p.delay = 0
    p.dt_list = ["chloroplast", "rbcL", "matK"]
    p._matrix_records = [
        {"Species": "Sp one", "taxid": "1",
         "counts": {"chloroplast": 5, "rbcL": 0, "matK": mm.COUNT_NA},
         "union": 5, "total_all": 10},
        {"Species": "Sp two", "taxid": "2",
         "counts": {"chloroplast": 0, "rbcL": 3, "matK": 2},
         "union": 4, "total_all": 8},
    ]
    p.species_list = ["Sp one", "Sp two"]
    return p


def test_heatmaps_written(monkeypatch, tmp_path):
    p = _matrix_program(monkeypatch, tmp_path)
    paths = p.save_heatmap("testts")
    assert len(paths) == 2
    for pth in paths:
        assert os.path.exists(pth)
        assert os.path.getsize(pth) > 0


def test_coverage_summary_and_manifest(monkeypatch, tmp_path):
    p = _matrix_program(monkeypatch, tmp_path)
    cov = p.save_coverage_summary("testts")
    man = p.save_manifest("testts")
    assert os.path.exists(cov)
    assert os.path.exists(man)
    import json
    with open(man) as f:
        manifest = json.load(f)
    assert manifest["entrez_db"] == "nuccore"
    assert manifest["datatypes"] == ["chloroplast", "rbcL", "matK"]
    assert "filters" in manifest


def test_end_to_end_writes_all_artifacts(monkeypatch, tmp_path):
    monkeypatch.setattr(mm, "RESULTS_DIR", str(tmp_path))
    p = mm.MissMap()
    p.delay = 0
    _wire_species(p, monkeypatch, _fake_counts())
    p.data_counts()
    p.display_and_save_matrix()

    files = os.listdir(tmp_path)
    assert any(f.startswith("data_availability_matrix_") and f.endswith(".csv") for f in files)
    assert any(f.startswith("data_availability_matrix_") and f.endswith(".xlsx") for f in files)
    assert any(f.startswith("missmap_presence_") for f in files)
    assert any(f.startswith("missmap_counts_") for f in files)
    assert any(f.startswith("coverage_summary_") for f in files)
    assert any(f.startswith("run_manifest_") for f in files)


def test_excel_engine_available():
    assert mm.MissMap._pick_excel_engine() in ("xlsxwriter", "openpyxl")


# ======================================================================
# Non-interactive guarding
# ======================================================================
def test_non_interactive_species_exits():
    p = mm.MissMap(non_interactive=True)
    with pytest.raises(SystemExit):
        p.ask_species()


def test_non_interactive_datatypes_defaults_to_all():
    p = mm.MissMap(non_interactive=True)
    p.ask_data_types()
    assert p.dt_list == list(mm.MissMap.FILTERS.keys())
