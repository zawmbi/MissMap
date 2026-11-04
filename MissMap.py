# MissMap.py
# Integrated AI-based synonym enrichment in a single file (Groq version + optional subspecies exclusion)

import os                   # file handling
import sys                  # error stuff
import time                 # to delay API calls
import argparse             # accepting terminal input
import pandas as pd         # dataframe for the output matrix
from Bio import Entrez      # accessing NCBI Entrez API
from urllib.error import HTTPError, URLError  # handles API request errors
from xml.etree import ElementTree as ET  # parse XML from NCBI
import requests             # Groq API calls
import json                 # parse JSON from LLM
import re                   # infraspecific-rank detection

# ======================================================================
#  MissMap: A pipeline for visualizing sequence data availability in
#  plant clades with AI-driven synonym enrichment
# ======================================================================

MISSMAP_INTRO = """
================================================================================#
  MissMap: A pipeline for visualizing sequence data availability in plant clades
  with AI-driven synonym enrichment
================================================================================#

MissMap utilizes the Entrez API from NCBI GenBank, thus recommends that you input 
your NCBI registered email address and API key in order to increase your API request 
rate from 3 requests per second to 10 requests per second. If you do not have an NCBI 
account/API, you can register on NCBI for free, or simply skip this step if you are not
performing a large number of requests. 

Manual entry with email and API key: 
$ python MissMap.py -e yourncbiemail@address.com -k yourapikeylettersandnumbers -s "species1,species2" -dt "chloroplast,chloroplast,rbcL"

Manual entry without email and API key using config file (recommended): 
$ python MissMap.py -s "species1" -dt "rbcL,matK,ITS" -c config.txt

Storing email & API key and Groq settings in your environment (config.txt or .env):
# email yourncbiemail@address.com
# GROQ_API_KEY
# GROQ_MODEL

You can also optionally set:
# species Arabidopsis thaliana,Zea mays
# datatypes chloroplast,mitochondrion,rbcL

Optional flags:
--verbose         Print per-datatype progress and count details.
--ai-unverified   Allow AI synonyms not validated in NCBI taxonomy.
--debug-queries   Print exact NCBI query strings used for counts.
--ai-validated-only  Use only NCBI-validated AI synonyms in searches (strict).
--no-subspecies   Exclude infraspecific/hybrid names from synonym search.
--non-interactive Never prompt; fail if required inputs are missing (HPC/batch/CI).

Outputs (in missmap_results/): matrix CSV/XLSX, presence/absence + count
heatmaps (PNG), a coverage summary CSV, and a JSON run manifest.

AI Disclaimer:
AI output can vary by model and input context. Synonym lists may differ across runs,
models, or prompt updates. Treat AI synonyms as suggestions to be verified.

================================================================================#
"""

RESULTS_DIR = "missmap_results"

# ----------------------------------------------------------------------
# Counting reliability constants
# ----------------------------------------------------------------------
# Sentinel for "the query failed and we do NOT know the count". This is kept
# strictly distinct from an integer 0 (a genuine "no records") so a transient
# network blip is never silently recorded as a biological absence.
COUNT_NA = None
NA_LABEL = "NA"

# HTTP status codes worth retrying (rate limiting + transient server errors).
RETRYABLE_HTTP = {429, 500, 502, 503, 504}
# Total attempts = 1 initial + MAX_RETRIES retries.
MAX_RETRIES = 4
# Exponential backoff: wait = BACKOFF_BASE * 2**attempt (seconds).
BACKOFF_BASE = 1.0


def _is_retryable_error(err):
    """True if an Entrez/urllib error is a transient condition worth retrying."""
    if isinstance(err, HTTPError):
        return getattr(err, "code", None) in RETRYABLE_HTTP
    if isinstance(err, URLError):
        return True  # DNS / connection reset / timeout at the network layer
    # Biopython raises RuntimeError carrying NCBI's own error text for some
    # rate-limit / server conditions; sniff the message as a last resort.
    msg = str(err).lower()
    transient_markers = (
        "429", "too many requests", "500", "502", "503", "504",
        "timeout", "timed out", "temporarily", "try again", "connection reset",
    )
    return any(t in msg for t in transient_markers)

# ======================================================================
# === AI-BASED SYNONYM ENRICHMENT HELPERS (Groq) =======================
# ======================================================================

def load_dotenv(dotenv_path=".env"):
    """Lightweight .env loader (KEY=VALUE). Only sets vars not already set."""
    try:
        with open(dotenv_path, "r") as f:
            for line in f:
                raw = line.strip()
                if not raw or raw.startswith("#"):
                    continue
                if "=" not in raw:
                    continue
                key, val = raw.split("=", 1)
                key = key.strip()
                val = val.strip().strip('"').strip("'")
                if key and key not in os.environ:
                    os.environ[key] = val
    except FileNotFoundError:
        # .env is optional
        pass

def read_config(fl):
    """Read KEY VALUE config lines into a dict (supports commented lines)."""
    data = {}
    try:
        with open(fl, "r") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                # allow commented config like '# key val'
                if line.startswith("#"):
                    line = line.lstrip("#").strip()
                    if not line:
                        continue
                parts = line.split()
                if len(parts) >= 2:
                    key, val = parts[0], parts[1]
                    data[key] = val
    except FileNotFoundError:
        # Missing config file is not fatal; handled later.
        pass
    return data


def get_llm_client(configfile):
    """
    Initialize and return a simple Groq config dict using config file or env vars.
    We keep the name 'get_llm_client' so the rest of the code doesn't change.
    """
    cfg = read_config(configfile)

    api_key = cfg.get("GROQ_API_KEY") or os.getenv("GROQ_API_KEY")
    model = cfg.get("GROQ_MODEL") or os.getenv("GROQ_MODEL") or ""

    if not api_key:
        raise KeyError(
            "Missing GROQ_API_KEY. Set it in config.txt or as an environment variable."
        )
    if not model:
        raise KeyError(
            "Missing GROQ_MODEL. Set it in config.txt, .env, or as an environment variable."
        )

    return {
        "api_key": api_key,
        "model": model,
    }


# Template for enriching synonyms via AI (unchanged)
AI_MAX_SYNONYMS = 200
SYN_TEMPLATE = (
    "You are a botanical taxonomy assistant. Given a species name and a list of known synonyms, "
    "return a comprehensive list of additional genuine taxonomic synonyms (not subspecies or varieties), "
    "plus well-established common names or spelling variants. Be exhaustive but limit to at most "
    f"{AI_MAX_SYNONYMS} total items.\n"
    "Do NOT include cultivars, processing terms, traits, or non-synonym species.\n"
    "Input JSON: {{\"species\": \"{species}\", \"known_synonyms\": {known}}}\n"
    "Output JSON schema:\n"
    "{{\"synonyms\": [\"syn1\", \"syn2\", ...]}}\n"
    "If none, return {{\"synonyms\": []}}.\n"
    "Respond with JSON ONLY, no extra text."
)

GROQ_URL = "https://api.groq.com/openai/v1/chat/completions"


def run_template_query(llm, species, known_list):
    """
    Call the Groq chat API using SYN_TEMPLATE.
    llm is a dict: {'api_key': ..., 'model': ...}
    """
    system_prompt = SYN_TEMPLATE.format(species=species, known=known_list)

    messages = [
        {"role": "system", "content": system_prompt},
        {"role": "user", "content": ""},  # mirrors the previous empty HumanMessage
    ]

    payload = {
        "model": llm["model"],
        "messages": messages,
        "temperature": 0,
        "max_tokens": 1500,
    }

    headers = {
        "Authorization": f"Bearer {llm['api_key']}",
        "Content-Type": "application/json",
    }

    try:
        resp = requests.post(GROQ_URL, headers=headers, json=payload, timeout=30)
    except Exception as e:
        print(f"[WARN] Groq request failed for '{species}': {e}")
        return {"synonyms": []}

    if resp.status_code != 200:
        print(f"[WARN] Groq API error for '{species}': {resp.status_code}")
        try:
            print("[DEBUG]", resp.text[:500])
        except Exception:
            pass
        return {"synonyms": []}

    try:
        data = resp.json()
        content = data["choices"][0]["message"]["content"]
    except Exception as e:
        print(f"[WARN] Unexpected Groq response format for '{species}': {e}")
        return {"synonyms": []}

    # Content should be raw JSON because of SYN_TEMPLATE.
    try:
        parsed = json.loads(content)
    except Exception as e:
        print(f"[WARN] Failed to parse AI synonym JSON for '{species}': {e}")
        print(f"[DEBUG] Raw AI response: {content}")
        return {"synonyms": []}

    return parsed  # expect {'synonyms': [...]}


def lookup_synonyms_with_ai(species, manual_list, llm):
    """Return AI-enriched synonyms by merging manual list and LLM output."""
    if not llm:
        return []
    known = manual_list or []
    try:
        result = run_template_query(llm, species, known)
        extra = result.get("synonyms", [])
        return [s.strip() for s in extra if isinstance(s, str) and s.strip()]
    except Exception as e:
        print(f"[WARN] AI synonym lookup failed for '{species}': {e}")
        return []

# ======================================================================
# === Configuration class ==============================================
# ======================================================================
class Configuration:
    def __init__(self, config_file):
        self.config_file       = config_file   # used by Terminal inputs, Entrez, and LLM creds
        self.email             = None
        self.api_key           = None
        self.species           = None         # list of species
        self.data_types        = None         # list of data types (strings)
        self.no_ai             = False
        self.exclude_infraspecific = False    # new: option to exclude subspecies/varieties
        self.ai_unverified     = False        # allow AI synonyms not validated by NCBI taxonomy
        self.verbose           = False        # verbose progress output
        self.debug_queries     = False        # print exact NCBI query strings
        self.ai_validated_only = False        # only use NCBI-validated AI synonyms in searches
        # Non-interactive: never block on input() (for HPC/batch/CI). Auto-on
        # when stdin is not a TTY; can also be forced with --non-interactive.
        self.non_interactive   = not sys.stdin.isatty()

    def load_config(self):
        """Load email, api_key, species, and datatypes from config file."""
        try:
            with open(self.config_file, "r") as r:
                for line in r:
                    line = line.strip()
                    if not line:
                        continue
                    if line.startswith("#"):
                        line = line.lstrip("#").strip()
                        if not line:
                            continue
                    parts = line.split(maxsplit=1)
                    if len(parts) < 2:
                        continue
                    key, val = parts[0].lower(), parts[1]
                    if key in ("email", "entrez_email"):
                        self.email = val
                    elif key in ("apikey", "api_key", "entrez_api_key"):
                        self.api_key = val
                    elif key in ("species",):
                        self.species = [s.strip() for s in val.split(",") if s.strip()]
                    elif key in ("datatypes", "data_types"):
                        self.data_types = [d.strip() for d in val.split(",") if d.strip()]
        except FileNotFoundError:
            # Config is optional; just skip if missing.
            pass

    def terminal_inputs(self):
        parser = argparse.ArgumentParser(
            prog="MissMap",
            description="GenBank data availability visualizer."
        )
        parser.add_argument("-e", "--email",  help="NCBI email", required=False)
        parser.add_argument("-k", "--apikey", help="NCBI API key", required=False)
        parser.add_argument("-s", "--species",
                            help="Comma-separated list of species", required=False)
        parser.add_argument("-dt", "--datatypes",
                            help="Comma-separated list of data types", required=False)
        parser.add_argument("-c", "--configfile", type=str, required=False,
                            help="Config file for Groq/OpenAI and Entrez")
        parser.add_argument("--no-ai", action="store_true",
                            help="Disable AI-based synonym enrichment.")
        parser.add_argument("--no-subspecies", action="store_true",
                            help="Exclude infraspecific names (e.g., 'var.' or 'subsp.') from synonym search.")
        parser.add_argument("--ai-unverified", action="store_true",
                            help="Allow AI synonyms not validated in NCBI taxonomy (name-based only).")
        parser.add_argument("--verbose", action="store_true",
                            help="Print per-datatype progress and count details.")
        parser.add_argument("--debug-queries", action="store_true",
                            help="Print exact NCBI query strings used for counts.")
        parser.add_argument("--ai-validated-only", action="store_true",
                            help="Use only NCBI-validated AI synonyms in searches (strict mode).")
        parser.add_argument("--non-interactive", action="store_true",
                            help="Never prompt; fail if required inputs are missing (HPC/batch/CI).")
        args = parser.parse_args()

        # Optional override of config file
        if args.configfile:
            self.config_file       = args.configfile
            # reload config with new file path
            self.load_config()

        # Non-interactive mode: CLI flag forces it; otherwise inherit the
        # TTY-based default set in __init__.
        if args.non_interactive:
            self.non_interactive = True

        # Email
        if args.email:
            self.email = args.email
        elif not self.email and os.getenv("ENTREZ_EMAIL"):
            self.email = os.getenv("ENTREZ_EMAIL")
        elif not self.email and not self.non_interactive:
            e = input("NCBI Email Address (press Enter to skip): ").strip()
            if e:
                self.email = e

        # API key
        if args.apikey:
            self.api_key = args.apikey
        elif not self.api_key and os.getenv("ENTREZ_API_KEY"):
            self.api_key = os.getenv("ENTREZ_API_KEY")
        elif not self.api_key and not self.non_interactive:
            k = input("API Key (press Enter to skip): ").strip()
            if k:
                self.api_key = k

        # Species from CLI (overrides config)
        if args.species:
            self.species = [s.strip() for s in args.species.split(",") if s.strip()]

        # Data types from CLI (overrides config)
        if args.datatypes:
            self.data_types = [d.strip() for d in args.datatypes.split(",") if d.strip()]

        # Flags
        self.no_ai = args.no_ai
        self.exclude_infraspecific = args.no_subspecies
        self.ai_unverified = args.ai_unverified
        self.verbose = args.verbose
        self.debug_queries = args.debug_queries
        self.ai_validated_only = args.ai_validated_only

        # Configure Entrez
        if self.email:
            Entrez.email = self.email
        if self.api_key:
            Entrez.api_key = self.api_key

# ======================================================================
# === MissMap class ====================================================
# ======================================================================
class MissMap:
    FILTERS = {
        "chloroplast":   "organelle:chloroplast",
        "mitochondrion": "mitochondrion[Filter]",
        "nuclear":       "biomol_genomic[PROP] NOT srcdb_refseq[PROP]",
        "transcriptome": "biomol_mrna[PROP] NOT srcdb_refseq[PROP]",
        "RefSeq":        "srcdb_refseq[PROP]",
        "WGS":           "wgs[All Fields]",
        "rbcL":          "\"rbcL\"[Gene]",
        "matK":          "\"matK\"[Gene]",
        "ITS":           "internal transcribed spacer[All Fields]",
        "trnL":          "\"trnL\"[Gene]",
        "psbA-trnH":     "\"psbA-trnH\"[All Fields]",
        "atpF-atpH":     "\"atpF\"[Gene] AND \"atpH\"[Gene]",
        "psbK-psbI":     "\"psbK-psbI\"[All Fields]",
        "rpoB":          "\"rpoB\"[Gene]",
        "rpoC1":         "\"rpoC1\"[Gene]",
        "ycf1":          "\"ycf1\"[Gene]",
        "ndhF":          "\"ndhF\"[Gene]",
    }

    # Datatypes that count large genomic resources (whole chromosomes, WGS
    # contigs, every mRNA). Their magnitudes dwarf barcode markers, so the
    # visualisation groups and scales them separately rather than adding them.
    GENOMIC_CATEGORIES = {
        "chloroplast", "mitochondrion", "nuclear", "transcriptome", "RefSeq", "WGS",
    }

    def __init__(self, exclude_infraspecific=False, ai_unverified=False, verbose=False,
                 debug_queries=False, ai_validated_only=False, non_interactive=False):
        os.makedirs(RESULTS_DIR, exist_ok=True)
        self.species_list = []
        self.dt_list      = []
        self.rows         = []
        self.llm          = None  # will be set by init_ai_client
        self.exclude_infraspecific = exclude_infraspecific
        self.ai_unverified = ai_unverified
        self.verbose = verbose
        self.debug_queries = debug_queries
        self.ai_validated_only = ai_validated_only
        self.non_interactive = non_interactive
        self._tax_name_cache = {}  # cache taxonomy lookups by name

        # Reliability / reproducibility bookkeeping.
        self._failed_cells = []     # list of (species, datatype) that returned NA
        self._matrix_records = []   # numeric per-species rows for the heatmap
        self.run_started = time.strftime("%Y-%m-%d %H:%M:%S %z")

        # Be a bit more conservative if no API key is set.
        self.delay = 0.12 if getattr(Entrez, "api_key", None) else 0.4

    def init_ai_client(self, config_file):
        """Initialize Groq LLM config for synonym enrichment (reusing config file path)."""
        try:
            self.llm = get_llm_client(config_file)
            print("[INFO] Groq client config initialized for synonym enrichment.")
        except Exception as e:
            print(f"[WARN] Could not initialize Groq client: {e}")
            print("[WARN] Continuing without AI-based synonym enrichment.")
            self.llm = None

    def lookup_synonyms_with_ai(self, species, manual_list):
        """Return extra synonyms using the LLM."""
        return lookup_synonyms_with_ai(species, manual_list, self.llm)

    def ask_species(self):
        if self.non_interactive:
            print("[ERROR] No species provided and running non-interactively. "
                  "Pass -s/--species or set 'species' in the config. Exiting.")
            sys.exit(2)
        inp = input("Enter species (comma-separated): ")
        for sp in inp.split(","):
            sp = sp.strip()
            if sp:
                self.species_list.append(sp)
        if not self.species_list:
            print("No species provided. Exiting.")
            sys.exit(1)

    def ask_data_types(self):
        if self.non_interactive:
            # Sensible default for batch use: run the full predefined set.
            print("[INFO] No datatypes provided; running all predefined datatypes (non-interactive).")
            self.dt_list = list(self.FILTERS.keys())
            return
        print("Predefined data types:", ", ".join(self.FILTERS.keys()))
        inp = input("Enter data types or press Enter for all: ")
        if not inp.strip():
            self.dt_list = list(self.FILTERS.keys())
        else:
            self.dt_list = [dt.strip() for dt in inp.split(",") if dt.strip()]

    def validate_data_types(self):
        """Validate dt_list against FILTERS to avoid KeyErrors."""
        if not self.dt_list:
            self.dt_list = list(self.FILTERS.keys())

        valid = []
        for dt in self.dt_list:
            if dt not in self.FILTERS:
                print(f"[WARN] Unknown data type '{dt}' – ignoring.")
            else:
                valid.append(dt)
        self.dt_list = valid

        if not self.dt_list:
            print("No valid data types remain after validation. Exiting.")
            sys.exit(1)

    def handle_errors(self, e):
        print(f"[ERROR] {e}")

    @staticmethod
    def _fmt(value):
        """Render a count for output: integers as-is, the NA sentinel as 'NA'."""
        return NA_LABEL if value is COUNT_NA else str(value)

    def _esearch_count(self, term, context=""):
        """
        Run an nuccore esearch and return an integer record count, or COUNT_NA
        (None) if the query fails after retries.

        This is the single seam through which every count flows. It:
          * distinguishes a genuine 0 from an error (returns None, never 0, on
            failure) so a network blip is never recorded as a real absence;
          * retries transient HTTP 429/5xx and network errors with exponential
            backoff, honouring NCBI's rate limits;
          * sleeps self.delay after a *successful* call to pace requests.
        """
        last_err = None
        for attempt in range(MAX_RETRIES + 1):
            try:
                with Entrez.esearch(db="nuccore", term=term, retmax=0) as q:
                    rec = Entrez.read(q)
                # A well-formed esearch response always carries "Count".
                # Its absence means we got something we can't trust -> treat as
                # an error rather than defaulting to 0.
                if "Count" not in rec:
                    raise ValueError("Entrez response missing 'Count' field")
                time.sleep(self.delay)  # pace only after a real success
                return int(rec["Count"])
            except Exception as e:  # noqa: BLE001 - we re-raise/inspect below
                last_err = e
                if attempt < MAX_RETRIES and _is_retryable_error(e):
                    wait = BACKOFF_BASE * (2 ** attempt)
                    if self.verbose:
                        print(
                            f"[RETRY] {context or term}: {type(e).__name__} "
                            f"({e}); attempt {attempt + 1}/{MAX_RETRIES}, "
                            f"waiting {wait:.1f}s"
                        )
                    time.sleep(wait)
                    continue
                break

        # Exhausted retries (or non-transient failure): record NA, not 0.
        self.handle_errors(f"{context or term}: query failed, recording NA ({last_err})")
        self._failed_cells.append(context or term)
        return COUNT_NA

    def lookup_synonyms_manually(self, taxon):
        syns = []
        taxid = None
        try:
            # taxonomy search by scientific name
            with Entrez.esearch(
                db="taxonomy",
                term=f"{taxon.lower()}[Scientific Name]"
            ) as q:
                res = Entrez.read(q)
            ids = res.get("IdList", [])

            # fallback: all names
            if not ids:
                with Entrez.esearch(
                    db="taxonomy",
                    term=f"{taxon.lower()}[All Names]"
                ) as q:
                    res = Entrez.read(q)
                ids = res.get("IdList", [])

            if not ids:
                return syns, taxid  # no synonyms, no taxid

            # Take the first taxonomy ID (could refine later if needed)
            taxid = ids[0]

            with Entrez.efetch(db="taxonomy", id=taxid, retmode="xml") as q:
                xml = q.read()

            root = ET.fromstring(xml)
            for e in root.findall(".//OtherNames/Synonym"):
                if e.text:
                    syns.append(e.text.strip())
            for e in root.findall(".//OtherNames/GenbankCommonName"):
                if e.text:
                    txt = e.text.strip()
                    if txt and txt not in syns:
                        syns.append(txt)
        except Exception as e:
            self.handle_errors(e)
        return syns, taxid


    # Infraspecific / hybrid rank connectors. Matched as whole tokens so we
    # catch 'var.', 'subsp.', 'ssp.', 'f.', 'forma', 'nothosubsp.', 'subvar.'
    # and the hybrid marker '×' regardless of position -- not just the two
    # space-padded substrings the original code looked for.
    _INFRASPECIFIC_RE = re.compile(
        r"(?:^|\s)(?:var|subsp|ssp|subvar|nothosubsp|nothovar|forma|subf|f|cv)\.?(?:\s|$)"
        r"|(?:^|\s)[×x](?:\s|$)"
        r"|×",  # multiplication sign anywhere (hybrid)
        re.IGNORECASE,
    )

    def _is_infraspecific(self, name):
        """Heuristically decide if a name denotes an infraspecific/hybrid taxon."""
        return bool(self._INFRASPECIFIC_RE.search(name))

    def _filter_infraspecific(self, names, species):
        """
        Optionally drop infraspecific/hybrid names from the list.
        Keeps the original species name regardless.
        """
        if not self.exclude_infraspecific:
            return names

        filtered = []
        for n in names:
            # Always keep the canonical species name
            if n.lower() == species.lower():
                filtered.append(n)
                continue
            if self._is_infraspecific(n):
                continue
            filtered.append(n)
        return filtered

    def _taxonomy_name_search(self, name):
        """Return taxonomy ID list for a given name using NCBI taxonomy search."""
        key = name.strip().lower()
        if not key:
            return []
        if key in self._tax_name_cache:
            return self._tax_name_cache[key]

        ids = []
        try:
            with Entrez.esearch(
                db="taxonomy",
                term=f"\"{name}\"[All Names]"
            ) as q:
                res = Entrez.read(q)
            ids = res.get("IdList", []) or []
        except Exception as e:
            self.handle_errors(e)
            ids = []

        self._tax_name_cache[key] = ids
        time.sleep(self.delay)
        return ids

    def _validate_ai_synonyms(self, species, taxid, ai_synonyms):
        """
        Keep only AI synonyms that resolve to the SAME NCBI TaxID.
        This ensures synonym accuracy using NCBI as the source of truth.
        """
        if not ai_synonyms:
            return []
        if self.ai_unverified:
            # Permissive mode: accept AI synonyms even if not present in NCBI taxonomy.
            return [n.strip() for n in ai_synonyms if isinstance(n, str) and n.strip()]
        if not taxid:
            print(
                f"[WARN] No TaxID for '{species}'; skipping AI synonyms to preserve accuracy."
            )
            return []

        validated = []
        for name in ai_synonyms:
            if not isinstance(name, str):
                continue
            n = name.strip()
            if not n:
                continue
            if n.lower() == species.lower():
                continue
            ids = self._taxonomy_name_search(n)
            if taxid in ids:
                validated.append(n)
        return validated

    def data_counts(self):
        """
        Query NCBI for each species × datatype, using TaxID when available, and collect counts.
        Also compare name-based search counts using NCBI taxonomy synonyms vs NCBI+AI synonyms.
        """
        for species in self.species_list:
            row = {"Species": species}

            # 1) Get manual synonyms + TaxID from taxonomy
            manual, taxid = self.lookup_synonyms_manually(species)

            # 2) AI enrichment (uses manual list only), then validate against NCBI TaxID (labeling only)
            ai_extra_raw = self.lookup_synonyms_with_ai(species, manual)
            ai_proposed = [s.strip() for s in (ai_extra_raw or []) if isinstance(s, str) and s.strip()]
            ai_validated = self._validate_ai_synonyms(species, taxid, ai_proposed)
            ai_rejected = [s for s in ai_proposed if s not in ai_validated]

            # include the original species name in the sets
            base_names = [species] + (manual or [])
            if self.ai_validated_only:
                expanded_names = base_names + (ai_validated or [])
            else:
                expanded_names = base_names + (ai_proposed or [])

            # deduplicate and strip
            base_names = sorted(set(n.strip() for n in base_names if n.strip()))
            expanded_names = sorted(set(n.strip() for n in expanded_names if n.strip()))

            # optionally drop subspecies/varieties (only affects display + AI, not TaxID)
            base_names = self._filter_infraspecific(base_names, species)
            expanded_names = self._filter_infraspecific(expanded_names, species)

            # save synonyms (without duplicating the original species name)
            syn_only_base = [n for n in base_names if n.lower() != species.lower()]
            syn_only_ai = [n for n in ai_validated if n.lower() != species.lower()]
            syn_only_all = [n for n in expanded_names if n.lower() != species.lower()]
            row["NCBI Synonyms"] = ",".join(syn_only_base)
            row["AI Proposed Synonyms"] = ",".join(ai_proposed)
            row["AI-validated Synonyms"] = ",".join(syn_only_ai)
            row["AI Rejected Synonyms"] = ",".join(ai_rejected)
            row["AI Validation Source"] = "NCBI Taxonomy"
            row["synonyms"] = ",".join(syn_only_all)

            # 3) Decide whether to use TaxID or fall back to name-based search
            use_taxid = taxid is not None

            # Name-based queries for comparison (NCBI synonyms vs NCBI+AI synonyms)
            def build_organism_query(names, fallback_species):
                if names:
                    q = " OR ".join([f'"{name}"[Organism]' for name in names])
                else:
                    q = f'"{fallback_species}"[Organism]'
                return f"({q})"

            organism_query_base = build_organism_query(base_names, species)
            organism_query_expanded = build_organism_query(expanded_names, species)
            # Save query strings for manual NCBI checks
            row["NCBI Query (no filters, primary)"] = (
                f"txid{taxid}[Organism:exp]" if taxid else organism_query_expanded
            )
            row["NCBI Query (no filters, NCBI synonyms)"] = organism_query_base
            row["NCBI Query (no filters, NCBI+AI)"] = organism_query_expanded

            # Name-based fallback query if TaxID is missing (used for main counts)
            organism_query = organism_query_expanded

            # numeric_counts holds the authoritative (TaxID-based when available)
            # per-datatype count for this species: int, or COUNT_NA on failure.
            numeric_counts = {}
            name_based_total_base = 0
            name_based_total_expanded = 0

            # 4) For each datatype, search nuccore
            for dt in self.dt_list:
                if self.verbose:
                    print(f"[INFO] {species} | {dt} | running search...")
                if use_taxid:
                    # TaxID-based search (preferred, more robust)
                    term = f"txid{taxid}[Organism:exp] AND {self.FILTERS[dt]}"
                else:
                    # Fallback: name-based search using species + synonyms
                    term = f"{organism_query} AND {self.FILTERS[dt]}"
                if self.debug_queries:
                    print(f"[QUERY] {species} | {dt} | primary: {term}")

                count = self._esearch_count(term, context=f"{species} | {dt} | primary")
                numeric_counts[dt] = count
                row[dt] = self._fmt(count)
                if self.verbose:
                    print(f"[INFO] {species} | {dt} | TaxID/primary count: {self._fmt(count)}")

                # Name-based comparison counts (NCBI synonyms vs NCBI+AI).
                # NA cells are skipped from the running totals rather than
                # silently treated as zero.
                base_count = self._esearch_count(
                    f"{organism_query_base} AND {self.FILTERS[dt]}",
                    context=f"{species} | {dt} | name-based NCBI",
                )
                expanded_count = self._esearch_count(
                    f"{organism_query_expanded} AND {self.FILTERS[dt]}",
                    context=f"{species} | {dt} | name-based NCBI+AI",
                )
                if isinstance(base_count, int):
                    name_based_total_base += base_count
                if isinstance(expanded_count, int):
                    name_based_total_expanded += expanded_count
                if self.verbose:
                    print(
                        f"[INFO] {species} | {dt} | Name-based NCBI: {self._fmt(base_count)} "
                        f"| NCBI+AI: {self._fmt(expanded_count)}"
                    )
                if self.debug_queries:
                    print(
                        f"[QUERY] {species} | {dt} | name-based NCBI: {organism_query_base} AND {self.FILTERS[dt]}"
                    )
                    print(
                        f"[QUERY] {species} | {dt} | name-based NCBI+AI: {organism_query_expanded} AND {self.FILTERS[dt]}"
                    )

            # 5) Aggregate across datatypes.
            #    The SUM double-counts: one record (e.g. a chloroplast genome)
            #    matches chloroplast, WGS, rbcL, matK, ycf1, ndhF at once. It is
            #    reported only as an overlapping tally, never as a record count.
            valid_counts = [c for c in numeric_counts.values() if c is not COUNT_NA]
            na_cells = sum(1 for c in numeric_counts.values() if c is COUNT_NA)
            sum_overlapping = sum(valid_counts)

            # 5b) DEDUPLICATED UNION: records matching AT LEAST ONE selected
            #     filter. OR-combining inside a single query lets NCBI dedupe,
            #     giving a true record count (always <= Total NCBI no-filters).
            filter_union = " OR ".join(f"({self.FILTERS[dt]})" for dt in self.dt_list)
            if use_taxid:
                union_term = f"txid{taxid}[Organism:exp] AND ({filter_union})"
            else:
                union_term = f"{organism_query} AND ({filter_union})"
            if self.debug_queries:
                print(f"[QUERY] {species} | union (dedup): {union_term}")
            union_count = self._esearch_count(union_term, context=f"{species} | union")

            row["Records w/ >=1 Selected Marker (deduplicated)"] = self._fmt(union_count)
            row["Sum of Marker Counts (overlapping, not a record count)"] = str(sum_overlapping)
            row["Cells Failed (NA)"] = str(na_cells)
            row["Name-based Total (NCBI Synonyms)"] = str(name_based_total_base)
            row["Name-based Total (NCBI+AI)"] = str(name_based_total_expanded)
            row["Name-based Total Gain (AI - NCBI)"] = str(
                name_based_total_expanded - name_based_total_base
            )

            # 6) Total NCBI (no filters)
            if use_taxid:
                base_term = f"txid{taxid}[Organism:exp]"
            else:
                base_term = organism_query
            total_all = self._esearch_count(base_term, context=f"{species} | no-filter primary")
            row["Total NCBI (no filters)"] = self._fmt(total_all)
            if self.debug_queries:
                print(f"[QUERY] {species} | no-filter primary: {base_term}")

            # Sanity invariant: a deduplicated subset cannot exceed the whole.
            if (isinstance(union_count, int) and isinstance(total_all, int)
                    and union_count > total_all):
                print(
                    f"[WARN] {species}: union count ({union_count}) exceeds unfiltered "
                    f"total ({total_all}); check filter/scoping consistency."
                )

            # Name-based unfiltered comparison (NCBI synonyms vs NCBI+AI)
            total_all_base = self._esearch_count(
                organism_query_base, context=f"{species} | no-filter name-based NCBI"
            )
            total_all_expanded = self._esearch_count(
                organism_query_expanded, context=f"{species} | no-filter name-based NCBI+AI"
            )
            row["Name-based NCBI (no filters)"] = self._fmt(total_all_base)
            row["Name-based NCBI+AI (no filters)"] = self._fmt(total_all_expanded)
            if isinstance(total_all_base, int) and isinstance(total_all_expanded, int):
                row["Name-based Gain (no filters)"] = str(total_all_expanded - total_all_base)
            else:
                row["Name-based Gain (no filters)"] = NA_LABEL
            if self.debug_queries:
                print(f"[QUERY] {species} | no-filter name-based NCBI: {organism_query_base}")
                print(f"[QUERY] {species} | no-filter name-based NCBI+AI: {organism_query_expanded}")

            # Record-keeping for the heatmap / coverage summary.
            self._matrix_records.append({
                "Species": species,
                "taxid": taxid,
                "counts": dict(numeric_counts),
                "union": union_count,
                "total_all": total_all,
            })

            if na_cells:
                print(f"[WARN] {species}: {na_cells} datatype cell(s) recorded as NA (query failure).")

            union_is_zero = (union_count == 0)
            if union_is_zero or (union_count is COUNT_NA and not valid_counts):
                if use_taxid:
                    print(
                        f'[WARN] No sequence data found for "{species}" (TaxID {taxid}) – '
                        f"check spelling, taxon resolution, or filters."
                    )
                else:
                    print(
                        f'[WARN] No sequence data found for "{species}" (name-based search) – '
                        f"check spelling or filters."
                    )

            self.rows.append(row)


    # ------------------------------------------------------------------
    # Ordered columns: genomic categories first, then barcode markers,
    # preserving the user's dt_list order within each group.
    # ------------------------------------------------------------------
    def _ordered_datatypes(self):
        genomic = [dt for dt in self.dt_list if dt in self.GENOMIC_CATEGORIES]
        markers = [dt for dt in self.dt_list if dt not in self.GENOMIC_CATEGORIES]
        return genomic, markers

    def save_heatmap(self, ts):
        """
        Render the species x datatype availability map(s) -- the "map" in
        MissMap. Produces two PNGs:
          * presence/absence/NA grid (neutralises the huge granularity gap
            between genome-scale categories and single-amplicon barcodes);
          * log10(count+1) magnitude grid for relative depth.
        Genomic categories and barcode markers are grouped with a divider so
        they are never read as an additive whole.
        """
        if not self._matrix_records:
            print("[INFO] No data to plot; skipping heatmap.")
            return []
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            from matplotlib.colors import ListedColormap, BoundaryNorm
            from matplotlib.patches import Patch
            import numpy as np
        except Exception as e:  # pragma: no cover - environment dependent
            print(f"[WARN] Plotting libraries unavailable; skipping heatmap: {e}")
            return []

        species = [r["Species"] for r in self._matrix_records]
        genomic, markers = self._ordered_datatypes()
        cols = genomic + markers
        if not cols or not species:
            print("[INFO] Nothing to plot; skipping heatmap.")
            return []
        divider = len(genomic)  # boundary between groups (0 if no genomic cols)

        n_rows, n_cols = len(species), len(cols)
        counts = np.full((n_rows, n_cols), np.nan)
        for i, rec in enumerate(self._matrix_records):
            for j, dt in enumerate(cols):
                c = rec["counts"].get(dt, COUNT_NA)
                counts[i, j] = np.nan if c is COUNT_NA else float(c)

        figw = max(6.0, 0.7 * n_cols + 3.0)
        figh = max(3.0, 0.55 * n_rows + 2.0)
        annotate = (n_rows * n_cols) <= 400  # avoid unreadable clutter
        saved = []

        def _decorate(ax):
            ax.set_xticks(range(n_cols))
            ax.set_xticklabels(cols, rotation=45, ha="right", fontsize=8)
            ax.set_yticks(range(n_rows))
            ax.set_yticklabels(species, fontstyle="italic", fontsize=8)
            ax.set_xticks([x - 0.5 for x in range(1, n_cols + 1)], minor=True)
            ax.set_yticks([y - 0.5 for y in range(1, n_rows + 1)], minor=True)
            ax.grid(which="minor", color="white", linewidth=1.0)
            ax.tick_params(which="minor", length=0)
            if 0 < divider < n_cols:
                ax.axvline(divider - 0.5, color="black", linewidth=2.0)

        # ---- Figure 1: presence / absence / NA -------------------------------
        # states: 0 absent, 1 present, 2 NA(query failed)
        state = np.where(np.isnan(counts), 2.0, (counts > 0).astype(float))
        cmap = ListedColormap(["#eef1f5", "#2c7fb8", "#e34a33"])
        norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5], cmap.N)
        fig, ax = plt.subplots(figsize=(figw, figh))
        ax.imshow(state, cmap=cmap, norm=norm, aspect="auto")
        _decorate(ax)
        if annotate:
            for i in range(n_rows):
                for j in range(n_cols):
                    c = counts[i, j]
                    label = "NA" if np.isnan(c) else f"{int(c):,}"
                    txt_color = "white" if (not np.isnan(c) and c > 0) else "#333333"
                    ax.text(j, i, label, ha="center", va="center",
                            fontsize=6.5, color=txt_color)
        ax.set_title("MissMap: sequence data presence / absence", fontsize=11)
        legend = [
            Patch(facecolor="#2c7fb8", label="present (>=1 record)"),
            Patch(facecolor="#eef1f5", label="absent (0 records)"),
            Patch(facecolor="#e34a33", label="NA (query failed)"),
        ]
        ax.legend(handles=legend, bbox_to_anchor=(1.01, 1.0),
                  loc="upper left", fontsize=7, frameon=False)
        if 0 < divider < n_cols:
            ax.text(0, 1.02, "genomic <-", transform=ax.transAxes, fontsize=7, color="#555")
            ax.text(1.0, 1.02, "-> markers", transform=ax.transAxes, fontsize=7,
                    color="#555", ha="right")
        path1 = os.path.join(RESULTS_DIR, f"missmap_presence_{ts}.png")
        fig.savefig(path1, dpi=150, bbox_inches="tight")
        plt.close(fig)
        saved.append(path1)
        print(f"Saved presence/absence map to {path1}")

        # ---- Figure 2: log10(count+1) magnitude ------------------------------
        log_counts = np.log10(counts + 1.0)  # NaN stays NaN
        fig, ax = plt.subplots(figsize=(figw, figh))
        cmap2 = plt.get_cmap("viridis").copy()
        cmap2.set_bad(color="#e34a33")  # NA cells show red
        im = ax.imshow(log_counts, cmap=cmap2, aspect="auto")
        _decorate(ax)
        if annotate:
            vmax = np.nanmax(log_counts) if np.isfinite(np.nanmax(log_counts)) else 1.0
            for i in range(n_rows):
                for j in range(n_cols):
                    c = counts[i, j]
                    label = "NA" if np.isnan(c) else f"{int(c):,}"
                    shade = 0.0 if np.isnan(log_counts[i, j]) else log_counts[i, j] / max(vmax, 1e-9)
                    txt_color = "white" if shade > 0.55 else "black"
                    ax.text(j, i, label, ha="center", va="center",
                            fontsize=6.5, color=txt_color)
        cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
        cbar.set_label("log10(count + 1)", fontsize=8)
        ax.set_title("MissMap: sequence data magnitude", fontsize=11)
        path2 = os.path.join(RESULTS_DIR, f"missmap_counts_{ts}.png")
        fig.savefig(path2, dpi=150, bbox_inches="tight")
        plt.close(fig)
        saved.append(path2)
        print(f"Saved count magnitude map to {path2}")
        return saved

    # Standard land-plant barcode markers, for the gap report.
    CORE_BARCODES = ["rbcL", "matK", "ITS", "trnL", "psbA-trnH"]

    def save_coverage_summary(self, ts):
        """
        Per-marker coverage across the species set (actionable for sampling
        decisions) plus a gap list for the standard barcodes. Writes a CSV and
        prints a concise table. Returns the CSV path.
        """
        if not self._matrix_records:
            return None
        rows = []
        for dt in self.dt_list:
            present = absent = na = 0
            for rec in self._matrix_records:
                c = rec["counts"].get(dt, COUNT_NA)
                if c is COUNT_NA:
                    na += 1
                elif c > 0:
                    present += 1
                else:
                    absent += 1
            known = present + absent
            pct = round(100.0 * present / known, 1) if known else NA_LABEL
            rows.append({
                "Datatype": dt,
                "Group": "genomic" if dt in self.GENOMIC_CATEGORIES else "marker",
                "Species with data": present,
                "Species without data": absent,
                "Species NA": na,
                "Coverage % (of known)": pct,
            })
        cov = pd.DataFrame(rows)
        path = os.path.join(RESULTS_DIR, f"coverage_summary_{ts}.csv")
        cov.to_csv(path, index=False)

        print("\n#==== Marker coverage across species ====#\n")
        print(cov.to_string(index=False))
        print("\n#========================================#\n")

        # Gap list: which species lack each core barcode.
        core = [b for b in self.CORE_BARCODES if b in self.dt_list]
        if core:
            print("Core-barcode gaps (species missing each marker):")
            for b in core:
                missing = [
                    rec["Species"] for rec in self._matrix_records
                    if rec["counts"].get(b, COUNT_NA) == 0
                ]
                if missing:
                    print(f"  {b}: {', '.join(missing)}")
                else:
                    print(f"  {b}: (none missing)")
            print("")
        print(f"Saved coverage summary to {path}")
        return path

    def save_manifest(self, ts):
        """
        Write a run manifest so results are reproducible and citable: run date,
        Entrez db/context, options, AI model, and any failed cells. Counts on
        NCBI change daily, so an undated number is not citable.
        """
        ai_model = self.llm.get("model") if isinstance(self.llm, dict) else None
        manifest = {
            "tool": "MissMap",
            "run_started": self.run_started,
            "run_saved": time.strftime("%Y-%m-%d %H:%M:%S %z"),
            "entrez_db": "nuccore",
            "entrez_email_set": bool(getattr(Entrez, "email", None)),
            "entrez_api_key_set": bool(getattr(Entrez, "api_key", None)),
            "request_delay_seconds": self.delay,
            "species": self.species_list,
            "datatypes": self.dt_list,
            "filters": {dt: self.FILTERS[dt] for dt in self.dt_list},
            "options": {
                "exclude_infraspecific": self.exclude_infraspecific,
                "ai_enabled": self.llm is not None,
                "ai_unverified": self.ai_unverified,
                "ai_validated_only": self.ai_validated_only,
                "non_interactive": self.non_interactive,
            },
            "ai_model": ai_model,
            "failed_cells": self._failed_cells,
            "note": "Counts reflect NCBI nuccore on the run date and change daily.",
        }
        path = os.path.join(RESULTS_DIR, f"run_manifest_{ts}.json")
        with open(path, "w") as f:
            json.dump(manifest, f, indent=2)
        print(f"Saved run manifest to {path}")
        return path

    @staticmethod
    def _pick_excel_engine():
        """Return an installed pandas Excel engine, or None if unavailable."""
        for eng in ("xlsxwriter", "openpyxl"):
            try:
                __import__(eng)
                return eng
            except ImportError:
                continue
        return None

    def _save_excel(self, df, ts):
        """Save the matrix as .xlsx with auto-fit columns, engine-agnostic."""
        engine = self._pick_excel_engine()
        if engine is None:
            print("[WARN] No Excel engine (xlsxwriter/openpyxl) installed; skipping .xlsx.")
            return None
        out_xlsx = os.path.join(RESULTS_DIR, f"data_availability_matrix_{ts}.xlsx")
        try:
            with pd.ExcelWriter(out_xlsx, engine=engine) as writer:
                df.to_excel(writer, index=False, sheet_name="MissMap")
                worksheet = writer.sheets["MissMap"]
                for i, col in enumerate(df.columns):
                    series = df[col].astype(str)
                    width = min(max([len(str(col))] + [len(v) for v in series]) + 2, 120)
                    if engine == "xlsxwriter":
                        worksheet.set_column(i, i, width)
                    else:  # openpyxl
                        from openpyxl.utils import get_column_letter
                        worksheet.column_dimensions[get_column_letter(i + 1)].width = width
            print(f"Saved results to {out_xlsx}")
            return out_xlsx
        except Exception as e:
            print(f"[WARN] Could not write Excel file: {e}")
            return None

    def display_and_save_matrix(self):
        df = pd.DataFrame(self.rows)
        print("\n#==== GenBank data availability matrix ====#\n")
        if not df.empty:
            print(df.to_string(index=False))
        else:
            print("[INFO] No rows to display.")
        print("\n#========================================#\n")

        if self._failed_cells:
            print(f"[WARN] {len(self._failed_cells)} query cell(s) failed and were "
                  f"recorded as NA (not 0). See run manifest for the list.\n")

        ts = time.strftime("%Y%m%d-%H%M%S")
        out_csv = os.path.join(RESULTS_DIR, f"data_availability_matrix_{ts}.csv")
        df.to_csv(out_csv, index=False)
        print(f"Saved results to {out_csv}")

        self._save_excel(df, ts)
        self.save_coverage_summary(ts)
        self.save_heatmap(ts)
        self.save_manifest(ts)

# ======================================================================
# === Main execution ===================================================
# ======================================================================
if __name__ == "__main__":
    print(MISSMAP_INTRO)

    # load .env if present
    load_dotenv()

    # parse config & credentials
    config = Configuration("config.txt")
    config.load_config()
    config.terminal_inputs()

    # build program & init AI client
    program = MissMap(
        exclude_infraspecific=config.exclude_infraspecific,
        ai_unverified=config.ai_unverified,
        verbose=config.verbose,
        debug_queries=config.debug_queries,
        ai_validated_only=config.ai_validated_only,
        non_interactive=config.non_interactive
    )
    if not config.no_ai:
        program.init_ai_client(config.config_file)
    else:
        print("[INFO] AI synonym enrichment disabled (--no-ai).")

    # species input
    if config.species:
        program.species_list = config.species
    else:
        program.ask_species()

    # data types input
    if config.data_types:
        program.dt_list = config.data_types
    else:
        program.ask_data_types()

    # validate datatypes to avoid KeyErrors
    program.validate_data_types()

    # run counts & print
    program.data_counts()
    program.display_and_save_matrix()
