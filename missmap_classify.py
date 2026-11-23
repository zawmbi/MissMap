"""
Classify candidate synonym names (from LLMs, authorities, anyone) by what they
mean for NCBI retrieval, and build the Entrez queries that measure hidden and
over-counted records.

Buckets (one per name):
  same_taxid                    NCBI already files the name under the species' TaxID
                                (searching it finds nothing new)
  other_taxid_confirmed         NCBI files it under a DIFFERENT TaxID and at least one
                                authority says it is a synonym of the species
                                -> records under that TaxID are HIDDEN from a TaxID search
  other_taxid_wrong_species     NCBI files it elsewhere and an authority says it is a
                                different accepted species -> adding it would contaminate
  other_taxid_unconfirmed       NCBI files it elsewhere; no loaded authority knows it
  valid_not_in_ncbi             not in NCBI Taxonomy, but an authority/GBIF knows the name
  not_found                     nowhere we can check (possible hallucination)
  not_a_scientific_name         common name, cultivar, malformed string
"""

from missmap_names import parse_name, looks_like_scientific_name

BUCKETS = (
    "same_taxid", "other_taxid_confirmed", "other_taxid_wrong_species",
    "other_taxid_unconfirmed", "valid_not_in_ncbi", "not_found", "not_a_scientific_name",
)


def authority_verdict(name, species, checklists):
    """
    Per checklist: 'synonym_of_species' | 'species_itself' | 'different_taxon' |
    'unplaced' | 'absent'.
    """
    out = {}
    for cl in checklists:
        target = {a.taxonID for a in cl.accepted_for(species)}
        accs = cl.accepted_for(name)
        if not cl.find(name):
            out[cl.name] = "absent"
        elif not accs:
            out[cl.name] = "unplaced"
        elif {a.taxonID for a in accs} & target:
            same = any(u.is_accepted and u.taxonID in target for u in cl.find(name))
            out[cl.name] = "species_itself" if same else "synonym_of_species"
        else:
            out[cl.name] = "different_taxon"
    return out


def classify_name(name, species, primary_taxid, ncbi_candidates, checklists=(),
                  name_exists=None, descendant_taxids=()):
    """
    name_exists: optional callable(name) -> True/False/None (e.g. GBIFBackbone.name_exists)
    ncbi_candidates: list of NcbiTaxon the name resolves to in NCBI (in-group only)
    descendant_taxids: taxids under the primary taxon (already covered by [Organism:exp])
    """
    info = {"name": name, "ncbi_taxids": [c.taxid for c in ncbi_candidates], "authority": {}}
    if not looks_like_scientific_name(name):
        info["bucket"] = "not_a_scientific_name"
        return info
    verdicts = authority_verdict(name, species, checklists) if checklists else {}
    info["authority"] = verdicts
    covered = {str(primary_taxid)} | {str(t) for t in descendant_taxids} if primary_taxid else set()
    ids = {str(c.taxid) for c in ncbi_candidates}
    known_somewhere = any(v != "absent" for v in verdicts.values())

    if ids and ids <= covered:
        info["bucket"] = "same_taxid"
    elif ids:
        if any(v in ("synonym_of_species", "species_itself") for v in verdicts.values()):
            info["bucket"] = "other_taxid_confirmed"
        elif any(v == "different_taxon" for v in verdicts.values()):
            info["bucket"] = "other_taxid_wrong_species"
        else:
            info["bucket"] = "other_taxid_unconfirmed"
        info["hidden_taxids"] = sorted(ids - covered)
    else:
        if not known_somewhere and name_exists is not None:
            ex = name_exists(name)
            info["gbif_exists"] = ex
            known_somewhere = bool(ex)
        info["bucket"] = "valid_not_in_ncbi" if known_somewhere else "not_found"
    return info


def taxid_union(taxids, exp=True):
    tag = "Organism:exp" if exp else "Organism:noexp"
    ids = sorted({str(t) for t in taxids})
    if not ids:
        return ""
    if len(ids) == 1:
        return f"txid{ids[0]}[{tag}]"
    return "(" + " OR ".join(f"txid{t}[{tag}]" for t in ids) + ")"


def hidden_records_query(hidden_taxids, primary_taxid, filter_term=None):
    """
    Records under the hidden TaxIDs that a search on the primary TaxID does
    NOT already return. NCBI evaluates the NOT server-side, so overlap (e.g.
    a record cross-filed under both) is never double counted.
    """
    if not hidden_taxids:
        return ""
    q = f"{taxid_union(hidden_taxids)} NOT txid{primary_taxid}[Organism:exp]"
    if filter_term:
        q = f"({q}) AND ({filter_term})"
    return q


def overcount_query(species, primary_taxid):
    """Free-text hits for the name that are NOT records of the taxon (precision loss)."""
    return f'"{species}"[All Fields] NOT txid{primary_taxid}[Organism:exp]'


def all_fields_query(species):
    return f'"{species}"[All Fields]'


def is_infraspecific_name(name):
    return parse_name(name).is_infraspecific
