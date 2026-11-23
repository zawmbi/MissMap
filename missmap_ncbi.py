"""
NCBI Taxonomy helpers for MissMap: TaxID disambiguation, offline taxdump
lookups, and a small on-disk cache for Entrez answers.

TaxID disambiguation
--------------------
A name search in NCBI Taxonomy can return several TaxIDs (cross-kingdom genus
homonyms such as Prunella / Morus / Oenanthe, or a name that is a synonym of
one taxon and the scientific name of another). Taking the first hit silently
picks one. `choose_taxon` instead:

  1. keeps only candidates whose lineage contains the expected group
     (e.g. 'Viridiplantae') when one is given;
  2. prefers a candidate whose *scientific name* equals the query over one that
     merely lists it as a synonym;
  3. prefers rank 'species' for a binomial (or the matching infraspecific rank);
  4. reports anything still ambiguous instead of guessing.
"""

import json
import os
import sqlite3
import time
from dataclasses import dataclass, field
from xml.etree import ElementTree as ET

from missmap_names import parse_name, normalise_key

# Name classes in NCBI OtherNames that are *scientific* alternative names.
SYNONYM_CLASSES = ("Synonym", "EquivalentName", "Includes")
COMMON_CLASSES = ("GenbankCommonName", "CommonName", "GenbankSynonym")


@dataclass
class NcbiTaxon:
    taxid: str
    scientific_name: str
    rank: str = ""
    lineage: list = field(default_factory=list)
    division: str = ""
    synonyms: list = field(default_factory=list)        # Synonym + EquivalentName + Includes
    authorities: list = field(default_factory=list)     # 'Genus epithet Author' strings
    misspellings: list = field(default_factory=list)
    common_names: list = field(default_factory=list)

    def in_group(self, group):
        g = group.lower()
        return any(x.lower() == g for x in self.lineage) or self.division.lower() == g


def parse_taxon_xml(xml):
    """Parse an efetch db=taxonomy XML payload into NcbiTaxon objects."""
    root = ET.fromstring(xml)
    out = []
    for t in root.findall("./Taxon"):
        tx = NcbiTaxon(
            taxid=(t.findtext("TaxId") or "").strip(),
            scientific_name=(t.findtext("ScientificName") or "").strip(),
            rank=(t.findtext("Rank") or "").strip(),
            division=(t.findtext("Division") or "").strip(),
        )
        lin = (t.findtext("Lineage") or "").strip()
        tx.lineage = [x.strip() for x in lin.split(";") if x.strip()]
        on = t.find("OtherNames")
        if on is not None:
            for cls in SYNONYM_CLASSES:
                tx.synonyms += [e.text.strip() for e in on.findall(cls) if e.text]
            tx.authorities = [e.text.strip() for e in on.findall("Authority") if e.text]
            tx.misspellings = [e.text.strip() for e in on.findall("Misspelling") if e.text]
            for cls in COMMON_CLASSES:
                tx.common_names += [e.text.strip() for e in on.findall(cls) if e.text]
            # <Name><ClassCDE>synonym</ClassCDE><DispName>..</DispName></Name> form
            for nm in on.findall("Name"):
                cls = (nm.findtext("ClassCDE") or "").strip().lower()
                disp = (nm.findtext("DispName") or "").strip()
                if disp and cls in ("synonym", "equivalent name", "includes") and disp not in tx.synonyms:
                    tx.synonyms.append(disp)
        out.append(tx)
    return out


_RANK_FOR_INFRA = {"subsp.": "subspecies", "var.": "varietas", "f.": "forma"}


def choose_taxon(query, candidates, expected_group=None):
    """
    Pick the TaxID for `query` from candidate NcbiTaxon objects.
    Returns (taxon_or_None, flags) where flags is a list of strings such as
    'homonym_other_group', 'ambiguous', 'matched_via_synonym', 'rank_mismatch'.
    """
    flags = []
    if not candidates:
        return None, ["not_found"]
    pool = list(candidates)
    if expected_group:
        in_group = [c for c in pool if c.in_group(expected_group)]
        if len(in_group) < len(pool):
            flags.append("homonym_other_group")
        pool = in_group
        if not pool:
            return None, flags + ["not_in_expected_group"]

    p = parse_name(query)
    qkey = p.key() if p.parsed else normalise_key(query)
    exact = [c for c in pool if normalise_key(parse_name(c.scientific_name).canonical
                                              or c.scientific_name) == qkey]
    if exact:
        pool = exact
    else:
        flags.append("matched_via_synonym")

    want_rank = _RANK_FOR_INFRA.get(p.infra_rank, "species") if p.parsed else None
    if want_rank:
        ranked = [c for c in pool if c.rank == want_rank]
        if ranked:
            pool = ranked
        else:
            flags.append("rank_mismatch")
    if len(pool) > 1:
        flags.append("ambiguous")
        return None, flags
    return pool[0], flags


# ----------------------------------------------------------------------
# Offline taxdump (names.dmp + nodes.dmp from ftp.ncbi.nih.gov/pub/taxonomy/)
# ----------------------------------------------------------------------
class NcbiTaxdump:
    """
    Minimal reader for the NCBI taxdump. Loads names.dmp and nodes.dmp into
    memory (~2-3 GB RAM for the full dump; pass `keep_genera` to load only
    names in those genera, which is enough for a study sample).
    """

    def __init__(self, directory, keep_genera=None):
        self.directory = directory
        self.parent, self.rank = {}, {}
        self.sci_name = {}
        self.names_by_taxid = {}
        self.taxids_by_key = {}
        keep = {g.lower() for g in keep_genera} if keep_genera else None
        names_path = os.path.join(directory, "names.dmp")
        nodes_path = os.path.join(directory, "nodes.dmp")
        self.version = {"directory": os.path.abspath(directory),
                        "names_modified": time.strftime("%Y-%m-%d", time.localtime(os.path.getmtime(names_path)))}
        wanted = set()
        with open(names_path, encoding="utf-8", errors="replace") as f:
            for line in f:
                parts = [x.strip() for x in line.rstrip("\t|\n").split("\t|\t")]
                if len(parts) < 4:
                    continue
                tid, name, _uniq, cls = parts[:4]
                if keep is not None and name.split(" ")[0].lower() not in keep:
                    continue
                wanted.add(tid)
                self.names_by_taxid.setdefault(tid, []).append((name, cls))
                if cls == "scientific name":
                    self.sci_name[tid] = name
                if cls in ("scientific name", "synonym", "equivalent name", "includes",
                           "authority", "misspelling", "genbank synonym"):
                    k = parse_name(name).key() or normalise_key(name)
                    self.taxids_by_key.setdefault(k, set()).add(tid)
        with open(nodes_path, encoding="utf-8", errors="replace") as f:
            for line in f:
                parts = [x.strip() for x in line.rstrip("\t|\n").split("\t|\t")]
                if len(parts) < 3:
                    continue
                tid, par, rank = parts[:3]
                self.parent[tid] = par
                self.rank[tid] = rank
        # lineage names need scientific names of ancestors even when filtered
        if keep is not None:
            self._fill_ancestor_names(names_path, wanted)

    def _fill_ancestor_names(self, names_path, wanted):
        need = set()
        for tid in wanted:
            cur, seen = tid, 0
            while cur in self.parent and self.parent[cur] != cur and seen < 100:
                cur = self.parent[cur]
                need.add(cur)
                seen += 1
        need -= set(self.sci_name)
        if not need:
            return
        with open(names_path, encoding="utf-8", errors="replace") as f:
            for line in f:
                parts = [x.strip() for x in line.rstrip("\t|\n").split("\t|\t")]
                if len(parts) >= 4 and parts[3] == "scientific name" and parts[0] in need:
                    self.sci_name[parts[0]] = parts[1]

    def lineage(self, taxid):
        out, cur, seen = [], self.parent.get(taxid), 0
        while cur and cur != "1" and seen < 100:
            out.append(self.sci_name.get(cur, cur))
            nxt = self.parent.get(cur)
            if nxt == cur:
                break
            cur, seen = nxt, seen + 1
        return list(reversed(out))

    def taxon(self, taxid):
        if taxid not in self.sci_name and taxid not in self.names_by_taxid:
            return None
        tx = NcbiTaxon(taxid=taxid, scientific_name=self.sci_name.get(taxid, ""),
                       rank=self.rank.get(taxid, ""), lineage=self.lineage(taxid))
        for name, cls in self.names_by_taxid.get(taxid, []):
            if cls in ("synonym", "equivalent name", "includes"):
                tx.synonyms.append(name)
            elif cls == "authority":
                tx.authorities.append(name)
            elif cls == "misspelling":
                tx.misspellings.append(name)
            elif cls in ("genbank common name", "common name"):
                tx.common_names.append(name)
        return tx

    def candidates(self, name):
        k = parse_name(name).key() or normalise_key(name)
        return [t for t in (self.taxon(x) for x in sorted(self.taxids_by_key.get(k, ()))) if t]

    def children(self, taxid):
        return [t for t, p in self.parent.items() if p == taxid and t != taxid]


# ----------------------------------------------------------------------
# On-disk cache for Entrez answers
# ----------------------------------------------------------------------
class EntrezCache:
    """
    SQLite cache keyed by (kind, query). Answers older than `max_age_hours`
    are ignored, so same-day reruns are identical and free while a new day
    re-queries (NCBI counts change daily). Set max_age_hours=None to never expire
    (useful for freezing a dataset for the thesis).
    """

    def __init__(self, path, max_age_hours=24.0):
        self.path = path
        self.max_age = max_age_hours
        d = os.path.dirname(os.path.abspath(path))
        os.makedirs(d, exist_ok=True)
        self.db = sqlite3.connect(path)
        self.db.execute("CREATE TABLE IF NOT EXISTS cache (kind TEXT, query TEXT, value TEXT, "
                        "ts REAL, PRIMARY KEY (kind, query))")
        self.db.commit()
        self.hits = self.misses = 0

    def get(self, kind, query):
        row = self.db.execute("SELECT value, ts FROM cache WHERE kind=? AND query=?",
                              (kind, query)).fetchone()
        if row is None or (self.max_age is not None and time.time() - row[1] > self.max_age * 3600):
            self.misses += 1
            return None
        self.hits += 1
        return json.loads(row[0])

    def set(self, kind, query, value):
        self.db.execute("INSERT OR REPLACE INTO cache VALUES (?, ?, ?, ?)",
                        (kind, query, json.dumps(value), time.time()))
        self.db.commit()

    def close(self):
        self.db.close()
