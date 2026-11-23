"""
Conservative scientific-name parsing and matching for MissMap.

Why this exists
---------------
Matching names across NCBI, WCVP, WFO, GBIF and LLM output is where false
"synonyms" sneak in. The rules here are deliberately conservative:

  * exact canonical match (genus + epithet[s] + rank marker) is the only
    automatic match;
  * authorship is compared separately and reported, never used to *force*
    a match;
  * an edit-distance match (<= 2 on the epithet, same genus) is returned as
    a *candidate* for review, never silently accepted;
  * names that differ only in a gender ending (-us/-a/-um) or -ii/-i are
    reported as orthographic-variant candidates.

No third-party dependencies, so it is easy to test and to port.
"""

import re
import unicodedata
from dataclasses import dataclass, field

# Rank markers, normalised to a canonical spelling.
_RANK_ALIASES = {
    "subsp": "subsp.", "ssp": "subsp.", "subspecies": "subsp.",
    "var": "var.", "variety": "var.",
    "subvar": "subvar.",
    "f": "f.", "fo": "f.", "forma": "f.", "form": "f.",
    "subf": "subf.",
    "nothosubsp": "nothosubsp.", "nothovar": "nothovar.",
    "cv": "cv.",
}
_HYBRID_MARKERS = {"×", "x"}

# Words that show a string is not a scientific name at all.
_NON_NAME_HINTS = re.compile(
    r"\b(cultivar|group|complex|sp\.|spp\.|cf\.|aff\.|sensu|s\.l\.|s\.str\.|clone|strain|isolate)\b",
    re.IGNORECASE,
)


def strip_accents(s):
    return "".join(c for c in unicodedata.normalize("NFKD", s) if not unicodedata.combining(c))


@dataclass
class ParsedName:
    raw: str
    genus: str = ""
    epithet: str = ""
    infra_rank: str = ""
    infra_epithet: str = ""
    authorship: str = ""
    hybrid: bool = False
    parsed: bool = False          # False => could not read as a binomial/trinomial
    warnings: list = field(default_factory=list)

    @property
    def canonical(self):
        """Genus + epithet(s) with rank marker, no authorship. Lower-case epithets."""
        if not self.parsed:
            return ""
        parts = [self.genus]
        if self.hybrid:
            parts.append("×")
        parts.append(self.epithet)
        if self.infra_epithet:
            if self.infra_rank:
                parts.append(self.infra_rank)
            parts.append(self.infra_epithet)
        return " ".join(parts)

    @property
    def species_canonical(self):
        """Binomial only (drops infraspecific part)."""
        if not self.parsed:
            return ""
        return f"{self.genus} {'× ' if self.hybrid else ''}{self.epithet}"

    @property
    def is_infraspecific(self):
        return bool(self.infra_epithet)

    def key(self):
        """Lookup key: canonical name, lower-cased, accents stripped."""
        return normalise_key(self.canonical)


def normalise_key(name):
    """Case-, accent- and whitespace-insensitive key for a canonical name."""
    s = strip_accents(name or "").lower().replace("×", "x ")
    s = re.sub(r"\s+", " ", s).strip()
    # treat the hybrid sign as a separate token so 'x' and '×' agree
    s = re.sub(r"\bx\s+", "x ", s)
    return s


def parse_name(raw):
    """
    Parse a scientific name string into its parts.

    Handles: 'Genus epithet', 'Genus epithet Author', 'Genus epithet subsp. infra Author',
    'Genus × epithet', '×Genus epithet' (hybrid genus is kept in genus), and
    bracketed basionym authors '(L.) Mill.'.
    """
    p = ParsedName(raw=raw or "")
    s = (raw or "").strip()
    s = re.sub(r"\s+", " ", s)
    if not s:
        p.warnings.append("empty")
        return p
    if _NON_NAME_HINTS.search(s):
        p.warnings.append("contains non-name qualifier")

    tokens = s.split(" ")
    i = 0
    genus = tokens[0]
    if genus in _HYBRID_MARKERS and len(tokens) > 1:   # '× Genus epithet'
        genus = "×" + tokens[1]
        i = 1
    if not re.match(r"^×?[A-Z][a-z\-]+$", strip_accents(genus)):
        p.warnings.append("genus not capitalised latin word")
        return p
    p.genus = genus
    i += 1
    if i >= len(tokens):
        p.warnings.append("genus only")
        return p

    if tokens[i] in _HYBRID_MARKERS:
        p.hybrid = True
        i += 1
        if i >= len(tokens):
            return p
    epi = tokens[i]
    if not re.match(r"^[a-z][a-z\-]+$", strip_accents(epi)):
        p.warnings.append("epithet not lower-case latin word")
        return p
    p.epithet = strip_accents(epi)
    i += 1

    # optional author before an infraspecific rank, then rank + infra epithet
    rest = tokens[i:]
    rank_idx = None
    for j, t in enumerate(rest):
        if t.rstrip(".").lower() in _RANK_ALIASES and j + 1 < len(rest) \
                and re.match(r"^[a-z][a-z\-]+$", strip_accents(rest[j + 1])):
            rank_idx = j
            break
    if rank_idx is not None:
        species_author = " ".join(rest[:rank_idx])
        p.infra_rank = _RANK_ALIASES[rest[rank_idx].rstrip(".").lower()]
        p.infra_epithet = strip_accents(rest[rank_idx + 1])
        infra_author = " ".join(rest[rank_idx + 2:])
        # The authorship that belongs to the full name is the terminal one.
        p.authorship = infra_author or ""
        if species_author and not infra_author:
            p.warnings.append("author only on species part")
    elif rest and _is_bare_infra_epithet(rest[0]):
        # Trinomial without a rank marker: standard in zoology ('Panthera tigris
        # altaica') and in GBIF canonical names ('Lycopersicon esculentum cerasiforme').
        p.infra_epithet = strip_accents(rest[0])
        p.authorship = " ".join(rest[1:])
    else:
        p.authorship = " ".join(rest)
    p.parsed = True
    return p


# Lower-case words that start author strings, not epithets.
_AUTHOR_PARTICLES = {"de", "da", "del", "della", "der", "den", "di", "du", "van", "von", "la", "le",
                     "ex", "et", "in", "f", "fil", "ms", "nom", "auct", "sensu", "non", "emend", "hort"}


def _is_bare_infra_epithet(tok):
    t = strip_accents(tok)
    return (re.match(r"^[a-z][a-z\-]{2,}$", t) is not None
            and t not in _AUTHOR_PARTICLES and t.rstrip(".") in t)


def normalise_authorship(a):
    """
    Loose author-string normalisation for comparison:
    strip accents, spaces, dots, 'ex' prefixes; keep bracketed basionym
    authors separately so '(L.) Mill.' and 'Mill.' compare as 'combination
    author matches, basionym author differs'.
    """
    a = strip_accents(a or "")
    m = re.match(r"^\((.*?)\)\s*(.*)$", a.strip())
    basionym, combo = (m.group(1), m.group(2)) if m else ("", a)

    def clean(x):
        x = re.sub(r"\b\d{4}\b", "", x)               # zoological years
        x = re.sub(r".*\bex\b\s*", "", x)               # 'A ex B' -> 'B'
        x = x.replace("&", ",").replace(" et ", ",")
        x = re.sub(r"[\s\.]", "", x).lower()
        return x.strip(",")
    return clean(basionym), clean(combo)


def authorship_relation(a, b):
    """'same' | 'basionym_differs' | 'different' | 'unknown' (either missing)."""
    if not (a or "").strip() or not (b or "").strip():
        return "unknown"
    ba, ca = normalise_authorship(a)
    bb, cb = normalise_authorship(b)
    if ca == cb and ba == bb:
        return "same"
    if ca == cb:
        return "basionym_differs"
    # 'L.' vs 'Linnaeus' style abbreviations: prefix match on the combo author
    if ca and cb and (ca.startswith(cb) or cb.startswith(ca)):
        return "same"
    return "different"


def levenshtein(a, b, max_dist=None):
    """Plain edit distance; returns max_dist+1 early when exceeded."""
    if a == b:
        return 0
    if abs(len(a) - len(b)) > (max_dist if max_dist is not None else 10 ** 9):
        return (max_dist or 0) + 1
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i]
        best = i
        for j, cb in enumerate(b, 1):
            v = min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb))
            cur.append(v)
            best = min(best, v)
        if max_dist is not None and best > max_dist:
            return max_dist + 1
        prev = cur
    return prev[-1]


_GENDER_ENDINGS = [("us", "a"), ("us", "um"), ("a", "um"), ("is", "e"), ("er", "ra"),
                   ("er", "rum"), ("ii", "i"), ("iae", "ae")]


def _gender_variant(e1, e2):
    for x, y in _GENDER_ENDINGS:
        for s1, s2 in ((x, y), (y, x)):
            if e1.endswith(s1) and e2.endswith(s2) and e1[: -len(s1)] == e2[: -len(s2)]:
                return True
    return False


@dataclass
class MatchResult:
    kind: str            # exact | orthographic_variant | fuzzy_candidate | no_match | unparsed
    distance: int = 0
    authorship: str = "unknown"

    @property
    def accepted(self):
        """Only exact canonical matches are accepted automatically."""
        return self.kind == "exact"


def match_names(a, b, max_fuzzy=2):
    """
    Compare two names conservatively. Different genera never match (a genus
    transfer is a *synonymy* fact that must come from an authority, not from
    string similarity).
    """
    pa, pb = (a if isinstance(a, ParsedName) else parse_name(a),
              b if isinstance(b, ParsedName) else parse_name(b))
    if not (pa.parsed and pb.parsed):
        return MatchResult("unparsed")
    auth = authorship_relation(pa.authorship, pb.authorship)
    if pa.key() == pb.key():
        return MatchResult("exact", 0, auth)
    if normalise_key(pa.genus) != normalise_key(pb.genus) or pa.hybrid != pb.hybrid:
        return MatchResult("no_match")
    if (pa.infra_rank, bool(pa.infra_epithet)) != (pb.infra_rank, bool(pb.infra_epithet)):
        return MatchResult("no_match")
    pairs = [(pa.epithet, pb.epithet)]
    if pa.infra_epithet:
        pairs.append((pa.infra_epithet, pb.infra_epithet))
    if all(x == y or _gender_variant(x, y) for x, y in pairs):
        return MatchResult("orthographic_variant", 0, auth)
    dist = sum(levenshtein(x, y, max_fuzzy) for x, y in pairs)
    if dist <= max_fuzzy:
        return MatchResult("fuzzy_candidate", dist, auth)
    return MatchResult("no_match", dist, auth)


def looks_like_scientific_name(s):
    """True if the string parses as a binomial/trinomial without red flags."""
    p = parse_name(s)
    return p.parsed and not any("non-name" in w for w in p.warnings)
