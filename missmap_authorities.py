"""
Taxonomic authority adapters for MissMap.

Each authority has its OWN identifier system. They are not interchangeable:

    NCBI Taxonomy   taxid            e.g. 4081 (Solanum lycopersicum)
    WCVP / POWO     plant_name_id    + powo_id / ipni_id (urn:lsid:ipni.org:names:...)
    WFO             taxonID          e.g. wfo-0000xxxxxx (records often carry an IPNI ID)
    GBIF backbone   usageKey         integer
    Catalogue of Life  ID            short alphanumeric

So "the same taxon across platforms" has to be established by a bridge:
  1. a shared external ID where one exists (IPNI ID links WCVP <-> WFO for plants);
  2. otherwise exact canonical name (+ authorship check) — the only bridge to NCBI,
     which does not store IPNI/WFO IDs;
  3. an independent crosswalk (Wikidata stores NCBI, GBIF, IPNI, POWO and WFO IDs
     on one item) used to *verify* 1-2, never as ground truth on its own.
Every crosswalk row records which bridge was used (`matchMethod`) so it can be audited.

Loaders read the bulk downloads (versioned files) rather than per-name web APIs,
so every run can be pinned to a dated release. The GBIF and Wikidata adapters are
online and meant for baselines and verification.
"""

import csv
import hashlib
import io
import json
import os
import sys
import time
import zipfile
from collections import defaultdict
from dataclasses import dataclass, asdict

import requests

from missmap_names import parse_name, normalise_key, authorship_relation

csv.field_size_limit(min(sys.maxsize, 2 ** 31 - 1))


# ----------------------------------------------------------------------
# Common model (Darwin Core Taxon terms)
# ----------------------------------------------------------------------
@dataclass
class NameUsage:
    taxonID: str
    scientificName: str              # canonical name, no authorship
    scientificNameAuthorship: str
    taxonRank: str                   # species, subspecies, variety, form, genus ...
    taxonomicStatus: str             # accepted | synonym | unplaced | other (normalised)
    acceptedNameUsageID: str         # == taxonID for accepted names
    family: str = ""
    source: str = ""
    externalIDs: dict = None         # e.g. {"ipni": "...", "powo": "..."}

    @property
    def is_accepted(self):
        return self.taxonomicStatus == "accepted"

    def to_dwc(self):
        d = asdict(self)
        d["externalIDs"] = json.dumps(self.externalIDs or {}, sort_keys=True)
        return d


def _norm_status(s):
    s = (s or "").strip().lower()
    if s in ("accepted", "valid", "accepted name"):
        return "accepted"
    if "synonym" in s or s in ("illegitimate", "invalid", "misapplied", "orthographic"):
        return "synonym"
    if s in ("unplaced", "unchecked", "doubtful", "uncertain", "provisionally accepted", ""):
        return "unplaced"
    return s


def _norm_rank(r):
    r = (r or "").strip().lower().rstrip(".")
    return {"var": "variety", "subsp": "subspecies", "f": "form", "forma": "form"}.get(r, r)


def file_fingerprint(path):
    """Size, mtime and the first 16 hex of the SHA-256: enough to pin a release in the manifest."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    st = os.stat(path)
    return {"path": os.path.abspath(path), "bytes": st.st_size,
            "modified": time.strftime("%Y-%m-%d", time.localtime(st.st_mtime)),
            "sha256_16": h.hexdigest()[:16]}


# ----------------------------------------------------------------------
# Checklist: an in-memory index over NameUsage rows
# ----------------------------------------------------------------------
class Checklist:
    def __init__(self, name, version=None):
        self.name = name
        self.version = version or {}
        self.by_id = {}
        self.by_key = defaultdict(list)
        self.synonyms_by_accepted = defaultdict(list)

    def add(self, u):
        u.source = self.name
        self.by_id[u.taxonID] = u
        self.by_key[normalise_key(u.scientificName)].append(u)
        if u.taxonomicStatus == "synonym" and u.acceptedNameUsageID:
            self.synonyms_by_accepted[u.acceptedNameUsageID].append(u)

    def __len__(self):
        return len(self.by_id)

    def find(self, name, authorship=None):
        """All usages whose canonical name equals `name`; filtered by authorship if given."""
        p = parse_name(name)
        key = p.key() if p.parsed else normalise_key(name)
        hits = list(self.by_key.get(key, []))
        auth = authorship or (p.authorship if p.parsed else "")
        if auth and len(hits) > 1:
            same = [u for u in hits if authorship_relation(auth, u.scientificNameAuthorship)
                    in ("same", "basionym_differs")]
            if same:
                hits = same
        return hits

    def accepted_for(self, name, authorship=None):
        """
        Resolve a name to the accepted usage(s) this checklist assigns it to.
        Returns a de-duplicated list (more than one => the name is ambiguous here,
        e.g. a homonym or a pro parte synonym).
        """
        out = {}
        for u in self.find(name, authorship):
            if u.is_accepted:
                out[u.taxonID] = u
            elif u.taxonomicStatus == "synonym" and u.acceptedNameUsageID in self.by_id:
                acc = self.by_id[u.acceptedNameUsageID]
                out[acc.taxonID] = acc
        return list(out.values())

    def synonymy(self, accepted_id, include_infraspecific=True):
        """Every synonym the checklist files under this accepted taxon."""
        syns = self.synonyms_by_accepted.get(accepted_id, [])
        if not include_infraspecific:
            syns = [s for s in syns if s.taxonRank in ("species", "")]
        return syns

    def status_of(self, name):
        """'accepted' | 'synonym' | 'unplaced' | 'absent' | 'ambiguous' for a bare name."""
        hits = self.find(name)
        if not hits:
            return "absent"
        statuses = {u.taxonomicStatus for u in hits}
        if len(statuses) > 1:
            return "ambiguous"
        return statuses.pop()


# ----------------------------------------------------------------------
# Loaders
# ----------------------------------------------------------------------
def _open_text(path, member_hint=None):
    """Open a plain file or the matching member of a zip (DwC-A downloads are zips)."""
    if zipfile.is_zipfile(path):
        z = zipfile.ZipFile(path)
        names = z.namelist()
        member = next((n for n in names if member_hint and n.endswith(member_hint)), None) \
            or next((n for n in names if n.endswith((".csv", ".txt", ".tsv"))), None)
        if member is None:
            raise ValueError(f"No tabular file found inside {path}")
        return io.TextIOWrapper(z.open(member), encoding="utf-8", errors="replace")
    return open(path, encoding="utf-8", errors="replace", newline="")


def _sniff_delimiter(first_line):
    counts = {d: first_line.count(d) for d in ("|", "\t", ",")}
    return max(counts, key=counts.get)


def load_wcvp(path, genera=None):
    """
    Load the WCVP names file (wcvp_names.csv, pipe-delimited; also accepts the
    rWCVP-exported CSV). Column names follow the Kew release:
      plant_name_id, taxon_name, taxon_authors, taxon_rank, taxon_status,
      accepted_plant_name_id, family, ipni_id, powo_id
    `genera`: optional set of genus names to keep (speeds up large loads).
    """
    cl = Checklist("WCVP", version=file_fingerprint(path))
    genera = {g.lower() for g in genera} if genera else None
    with _open_text(path, "wcvp_names.csv") as f:
        first = f.readline()
        delim = _sniff_delimiter(first)
        header = next(csv.reader([first], delimiter=delim))
        reader = csv.DictReader(f, fieldnames=header, delimiter=delim)
        for r in reader:
            name = (r.get("taxon_name") or "").strip()
            if not name:
                continue
            if genera and name.split(" ")[0].lower() not in genera:
                continue
            tid = (r.get("plant_name_id") or "").strip()
            acc = (r.get("accepted_plant_name_id") or "").strip()
            status = _norm_status(r.get("taxon_status"))
            ext = {k: r[k].strip() for k in ("ipni_id", "powo_id") if (r.get(k) or "").strip()}
            cl.add(NameUsage(
                taxonID=tid, scientificName=name,
                scientificNameAuthorship=(r.get("taxon_authors") or "").strip(),
                taxonRank=_norm_rank(r.get("taxon_rank")), taxonomicStatus=status,
                acceptedNameUsageID=tid if status == "accepted" else acc,
                family=(r.get("family") or "").strip(), externalIDs=ext,
            ))
    return cl


def load_dwc(path, name, genera=None):
    """
    Load any Darwin Core Taxon table (WFO 'classification.csv', COL/ChecklistBank
    DwC-A 'Taxon.tsv', MDD/ASW exports converted to DwC). Required columns:
    taxonID, scientificName, taxonomicStatus, acceptedNameUsageID.
    `scientificName` may include authorship; it is split off using
    scientificNameAuthorship when present.
    """
    cl = Checklist(name, version=file_fingerprint(path))
    genera = {g.lower() for g in genera} if genera else None
    with _open_text(path, "classification.csv") as f:
        first = f.readline()
        delim = _sniff_delimiter(first)
        header = [h.split(":")[-1] for h in next(csv.reader([first], delimiter=delim))]  # 'dwc:taxonID'
        reader = csv.DictReader(f, fieldnames=header, delimiter=delim)
        for r in reader:
            sci = (r.get("scientificName") or "").strip()
            auth = (r.get("scientificNameAuthorship") or "").strip()
            if auth and sci.endswith(auth):
                sci = sci[: -len(auth)].strip()
            if not sci:
                continue
            if genera and sci.split(" ")[0].lower() not in genera:
                continue
            tid = (r.get("taxonID") or "").strip()
            status = _norm_status(r.get("taxonomicStatus"))
            acc = (r.get("acceptedNameUsageID") or "").strip()
            ext = {}
            sid = (r.get("scientificNameID") or "").strip()
            if "ipni" in sid.lower():
                ext["ipni"] = sid.split(":")[-1]
            elif sid:
                ext["scientificNameID"] = sid
            cl.add(NameUsage(
                taxonID=tid, scientificName=sci, scientificNameAuthorship=auth,
                taxonRank=_norm_rank(r.get("taxonRank")), taxonomicStatus=status,
                acceptedNameUsageID=tid if status == "accepted" else acc,
                family=(r.get("family") or "").strip(), externalIDs=ext,
            ))
    return cl


def load_wfo(path, genera=None):
    return load_dwc(path, "WFO", genera=genera)


# ----------------------------------------------------------------------
# Online adapters (baselines / verification; cache their answers)
# ----------------------------------------------------------------------
GBIF_API = "https://api.gbif.org/v1"
WIKIDATA_SPARQL = "https://query.wikidata.org/sparql"

# Wikidata properties holding each authority's identifier.
WIKIDATA_ID_PROPS = {
    "ncbi": "P685",    # NCBI taxonomy ID
    "gbif": "P846",    # GBIF taxon ID
    "ipni": "P961",    # IPNI plant ID
    "powo": "P5037",   # Plants of the World Online ID
    "wfo": "P7715",    # World Flora Online ID
}


def _get_json(url, params=None, headers=None, timeout=30, retries=3):
    last = None
    for attempt in range(retries + 1):
        try:
            r = requests.get(url, params=params, headers=headers, timeout=timeout)
            if r.status_code in (429, 500, 502, 503, 504):
                raise requests.HTTPError(f"HTTP {r.status_code}")
            r.raise_for_status()
            return r.json()
        except Exception as e:  # noqa: BLE001
            last = e
            if attempt < retries:
                time.sleep(1.0 * 2 ** attempt)
    raise last


class GBIFBackbone:
    """GBIF backbone via the public species API. Used as the rule-based (no-LLM) baseline."""

    def __init__(self, kingdom=None):
        self.kingdom = kingdom
        self._cache = {}

    def match(self, name):
        key = ("match", name, self.kingdom)
        if key not in self._cache:
            params = {"name": name, "strict": "true", "verbose": "false"}
            if self.kingdom:
                params["kingdom"] = self.kingdom
            self._cache[key] = _get_json(f"{GBIF_API}/species/match", params)
        return self._cache[key]

    def synonyms(self, usage_key):
        key = ("syn", usage_key)
        if key not in self._cache:
            out, offset = [], 0
            while True:
                page = _get_json(f"{GBIF_API}/species/{usage_key}/synonyms",
                                 {"limit": 100, "offset": offset})
                out.extend(page.get("results", []))
                if page.get("endOfRecords", True):
                    break
                offset += 100
            self._cache[key] = out
        return self._cache[key]

    _RANK_MARKER = {"SUBSPECIES": "subsp.", "VARIETY": "var.", "FORM": "f.",
                    "SUBVARIETY": "subvar.", "SUBFORM": "subf."}

    @classmethod
    def _canonical_with_rank(cls, rec):
        """GBIF canonicalName drops the rank marker ('A b c'); put it back ('A b var. c')."""
        cn = (rec.get("canonicalName") or "").strip()
        marker = cls._RANK_MARKER.get((rec.get("rank") or "").upper())
        parts = cn.split(" ")
        if marker and len(parts) == 3:
            return f"{parts[0]} {parts[1]} {marker} {parts[2]}"
        return cn

    def synonym_names(self, name):
        """(accepted canonical name, [synonym canonical names]) or (None, [])."""
        m = self.match(name)
        if m.get("matchType") not in ("EXACT",):
            return None, []
        acc_key = m.get("acceptedUsageKey") or m.get("usageKey")
        acc_name = m.get("species") if m.get("acceptedUsageKey") else m.get("canonicalName")
        names = [self._canonical_with_rank(s) for s in self.synonyms(acc_key) if s.get("canonicalName")]
        return acc_name, sorted(set(n for n in names if n))

    def name_exists(self, name):
        """True/False/None(unknown) — used as a hallucination check."""
        try:
            return self.match(name).get("matchType") == "EXACT"
        except Exception:  # noqa: BLE001
            return None


def wikidata_ids_for_name(name, user_agent="MissMap/1.0 (taxonomy research tool)"):
    """
    Look up a taxon name on Wikidata and return the authority IDs recorded on
    each matching item: [{"item": "Q...", "ncbi": "4081", "gbif": "...", ...}].
    More than one item => homonym or duplicate items; treat as unverified.
    """
    selects = " ".join(f"?{k}" for k in WIKIDATA_ID_PROPS)
    optionals = " ".join(f"OPTIONAL {{ ?item wdt:{p} ?{k} . }}" for k, p in WIKIDATA_ID_PROPS.items())
    safe = name.replace('"', "")
    q = f'SELECT ?item {selects} WHERE {{ ?item wdt:P225 "{safe}" . {optionals} }}'
    data = _get_json(WIKIDATA_SPARQL, {"query": q, "format": "json"},
                     headers={"User-Agent": user_agent, "Accept": "application/sparql-results+json"})
    rows = {}
    for b in data.get("results", {}).get("bindings", []):
        item = b["item"]["value"].rsplit("/", 1)[-1]
        d = rows.setdefault(item, {"item": item})
        for k in WIKIDATA_ID_PROPS:
            if k in b:
                d.setdefault(k, set()).add(b[k]["value"])
    return [{k: (sorted(v) if isinstance(v, set) else v) for k, v in d.items()} for d in rows.values()]


def verify_crosswalk(name, ncbi_taxid=None, ipni_id=None, wfo_id=None, lookup=wikidata_ids_for_name):
    """
    Independent check of a crosswalk row against Wikidata. Returns
    {'status': 'consistent'|'conflict'|'unverified', 'detail': ...}.
    'conflict' means Wikidata links this name to a different ID for at least
    one authority — those rows need a human look.
    """
    try:
        items = lookup(name)
    except Exception as e:  # noqa: BLE001
        return {"status": "unverified", "detail": f"lookup failed: {e}"}
    if not items:
        return {"status": "unverified", "detail": "name not on Wikidata"}
    if len(items) > 1:
        return {"status": "unverified", "detail": f"{len(items)} Wikidata items share this name"}
    it = items[0]
    conflicts, checked = [], 0
    for key, mine in (("ncbi", ncbi_taxid), ("ipni", ipni_id), ("wfo", wfo_id)):
        if mine and it.get(key):
            checked += 1
            if str(mine) not in it[key]:
                conflicts.append(f"{key}: ours={mine} wikidata={','.join(it[key])}")
    if conflicts:
        return {"status": "conflict", "detail": "; ".join(conflicts), "item": it["item"]}
    if checked == 0:
        return {"status": "unverified", "detail": "no overlapping IDs to compare", "item": it["item"]}
    return {"status": "consistent", "detail": f"{checked} ID(s) agree", "item": it["item"]}


# ----------------------------------------------------------------------
# Cross-authority comparison for one species (Aim 1 unit of analysis)
# ----------------------------------------------------------------------
def compare_synonymy(species, checklists, ncbi_synonyms=None, include_infraspecific=False):
    """
    For one species, gather each source's synonym set and classify every name.

    Returns dict with:
      accepted:   {source: [accepted names this source resolves `species` to]}
      synonyms:   {source: sorted list of synonym canonical names}
      jaccard:    {"A|B": float} pairwise over synonym sets (species-rank by default)
      names:      [{name, <source>: status_for_species, ..., class}]
    Status per source for each name: 'synonym' (of this species), 'accepted_other'
    (the source treats it as a different accepted taxon), 'synonym_other',
    'unplaced', 'absent', or 'self'.
    """
    sets, accepted = {}, {}
    acc_ids = {}
    for cl in checklists:
        accs = cl.accepted_for(species)
        accepted[cl.name] = sorted({a.scientificName for a in accs})
        acc_ids[cl.name] = {a.taxonID for a in accs}
        names = set()
        for a in accs:
            for s in cl.synonymy(a.taxonID, include_infraspecific):
                names.add(s.scientificName)
        sets[cl.name] = names
    if ncbi_synonyms is not None:
        keep = []
        for n in ncbi_synonyms:
            p = parse_name(n)
            if p.parsed and (include_infraspecific or not p.is_infraspecific):
                keep.append(p.canonical)
        sets["NCBI"] = set(keep)
        accepted["NCBI"] = [species]

    universe = sorted(set().union(*sets.values()) - {species}) if sets else []
    rows = []
    for n in universe:
        row = {"name": n}
        for cl in checklists:
            if n in sets[cl.name]:
                row[cl.name] = "synonym"
                continue
            st = cl.status_of(n)
            if st == "accepted":
                row[cl.name] = "accepted_other"
            elif st == "synonym":
                other = {a.taxonID for a in cl.accepted_for(n)}
                row[cl.name] = "synonym" if other & acc_ids[cl.name] else "synonym_other"
            else:
                row[cl.name] = st
        if "NCBI" in sets:
            row["NCBI"] = "synonym" if n in sets["NCBI"] else "not_synonym"
        statuses = [v for k, v in row.items() if k != "name"]
        n_syn = sum(v == "synonym" for v in statuses)
        if n_syn == len(statuses):
            row["class"] = "agreed"
        elif any(v in ("accepted_other", "synonym_other") for v in statuses):
            row["class"] = "lump_split"            # one source merges what another separates
        elif any(v == "absent" for v in statuses):
            row["class"] = "name_missing"          # a source simply does not know the name
        else:
            row["class"] = "status_differs"        # e.g. unplaced in one, synonym in another
        rows.append(row)

    jac = {}
    keys = sorted(sets)
    for i, a in enumerate(keys):
        for b in keys[i + 1:]:
            u = sets[a] | sets[b]
            jac[f"{a}|{b}"] = (len(sets[a] & sets[b]) / len(u)) if u else 1.0
    accepted_agree = len({tuple(v) for v in accepted.values() if v}) <= 1
    return {"species": species, "accepted": accepted, "accepted_agree": accepted_agree,
            "synonyms": {k: sorted(v) for k, v in sets.items()}, "jaccard": jac, "names": rows}


def authority_synonyms(species, checklists, include_infraspecific=False):
    """Union of synonym names for `species` across the given checklists, with sources."""
    out = defaultdict(set)
    for cl in checklists:
        for a in cl.accepted_for(species):
            for s in cl.synonymy(a.taxonID, include_infraspecific):
                out[s.scientificName].add(cl.name)
    return {k: sorted(v) for k, v in out.items()}
