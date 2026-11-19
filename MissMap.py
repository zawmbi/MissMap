# MissMap.py
# NCBI sequence-data availability per species, with cross-authority synonym
# reconciliation (WCVP / WFO / GBIF), hidden-record accounting, and optional
# LLM synonym enrichment (any provider; see missmap_llm.py).

import os                   # file handling
import sys                  # error stuff
import time                 # to delay API calls
import argparse             # accepting terminal input
import pandas as pd         # dataframe for the output matrix
from Bio import Entrez      # accessing NCBI Entrez API
from urllib.error import HTTPError, URLError  # handles API request errors
from xml.etree import ElementTree as ET  # parse XML from NCBI
import json                 # manifest + logs
import re                   # infraspecific-rank detection
from dataclasses import asdict

from missmap_llm import make_provider, query_synonyms, synonym_names, LLMLog, PROMPT_VERSION
from missmap_ncbi import parse_taxon_xml, choose_taxon, NcbiTaxdump, NcbiTaxon, EntrezCache
from missmap_authorities import (load_wcvp, load_wfo, load_dwc, GBIFBackbone,
                                 authority_synonyms, compare_synonymy)
from missmap_classify import (classify_name, hidden_records_query, overcount_query,
                              all_fields_query, BUCKETS)

# ======================================================================
#  MissMap: how much NCBI sequence data does a name hide?
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

Storing email, API keys and the LLM choice (config.txt or .env):
# email yourncbiemail@address.com
# LLM_MODEL anthropic:claude-sonnet-5-5     (or groq:llama-3.3-70b-versatile, ...)
# GROQ_API_KEY / ANTHROPIC_API_KEY / OPENROUTER_API_KEY as needed

You can also optionally set:
# species Arabidopsis thaliana,Zea mays
# datatypes chloroplast,mitochondrion,rbcL

Authorities and reconciliation (strongly recommended):
--wcvp PATH        WCVP names file (wcvp_names.csv or the release zip)
--wfo PATH         WFO classification.csv (DwC) or its zip
--dwc NAME=PATH    any other Darwin Core checklist (repeatable)
--gbif             use the GBIF backbone API as an extra authority/name check
--group NAME       expected lineage group for TaxID disambiguation (e.g. Viridiplantae)
--taxdump DIR      use a local NCBI taxdump (names.dmp/nodes.dmp) instead of Entrez lookups

Other flags:
--llm SPEC        provider:model for synonym enrichment (overrides config)
--verbose         Print per-datatype progress and count details.
--ai-unverified   Search with every AI name, even ones no authority confirms.
--debug-queries   Print exact NCBI query strings used for counts.
--ai-validated-only  Search only with AI names in the same TaxID or authority-confirmed.
--no-subspecies   Exclude infraspecific/hybrid names from synonym search.
--non-interactive Never prompt; fail if required inputs are missing (HPC/batch/CI).
--no-cache        Do not reuse same-day Entrez answers (missmap_results/entrez_cache.sqlite).

Outputs (in missmap_results/): matrix CSV/XLSX, presence/absence + count
heatmaps (PNG, before and after synonym reconciliation), a coverage summary CSV,
a per-name classification CSV, a Darwin Core crosswalk CSV, the LLM call log
(JSONL), and a JSON run manifest.

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
# === Config + AI-based synonym enrichment helpers ====================
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


def get_llm_client(configfile, spec=None):
    """
    Build the LLM provider used for synonym enrichment.

    Resolution order for the model:
      1. `spec` argument (from --llm), e.g. "anthropic:claude-sonnet-5-5"
      2. LLM_MODEL in config.txt / .env / environment (same spec format)
      3. legacy GROQ_MODEL (+ GROQ_API_KEY) -> "groq:<GROQ_MODEL>"
    """
    cfg = read_config(configfile)
    for k in ("GROQ_API_KEY", "OPENAI_API_KEY", "OPENROUTER_API_KEY",
              "TOGETHER_API_KEY", "ANTHROPIC_API_KEY"):
        if cfg.get(k) and not os.getenv(k):
            os.environ[k] = cfg[k]
    spec = spec or cfg.get("LLM_MODEL") or os.getenv("LLM_MODEL")
    if not spec:
        groq_model = cfg.get("GROQ_MODEL") or os.getenv("GROQ_MODEL")
        if not groq_model:
            raise KeyError("No LLM configured. Set LLM_MODEL (e.g. anthropic:claude-sonnet-5-5) "
                           "or GROQ_MODEL + GROQ_API_KEY, or pass --llm.")
        spec = f"groq:{groq_model}"
    return make_provider(spec)


def lookup_synonyms_with_ai(species, manual_list, llm, log=None, run_index=0):
    """
    Ask the LLM for synonyms. Returns (scientific_synonyms, record) where
    record is the logged call (raw text, model ID, errors) for the manifest.
    A failed call returns ([], record-with-error) -- it is never silently
    read as "the model knows no synonyms".
    """
    if not llm:
        return [], None
    rec = query_synonyms(llm, species, known=manual_list or None, run_index=run_index, log=log)
    if rec.get("error"):
        print(f"[WARN] AI synonym lookup failed for '{species}': {rec['error']}")
    return synonym_names(rec), rec

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
        # Authorities / reconciliation
        self.wcvp_path         = None
        self.wfo_path          = None
        self.dwc_paths         = {}           # {name: path}
        self.use_gbif          = False
        self.group             = None         # expected lineage group, e.g. Viridiplantae
        self.taxdump_dir       = None
        self.llm_spec          = None
        self.use_cache         = True
        self.cache_hours       = 24.0

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
                    elif key == "wcvp":
                        self.wcvp_path = val.strip()
                    elif key == "wfo":
                        self.wfo_path = val.strip()
                    elif key == "group":
                        self.group = val.strip()
                    elif key == "taxdump":
                        self.taxdump_dir = val.strip()
        except FileNotFoundError:
            # Config is optional; just skip if missing.
            pass

    def terminal_inputs(self, argv=None):
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
        parser.add_argument("--wcvp", help="Path to the WCVP names file (or release zip).")
        parser.add_argument("--wfo", help="Path to WFO classification.csv (or zip).")
        parser.add_argument("--dwc", action="append", default=[], metavar="NAME=PATH",
                            help="Extra Darwin Core checklist, e.g. MDD=mdd_taxon.tsv (repeatable).")
        parser.add_argument("--gbif", action="store_true",
                            help="Use the GBIF backbone API as an extra name check.")
        parser.add_argument("--group",
                            help="Expected lineage group for TaxID disambiguation (e.g. Viridiplantae, Mammalia).")
        parser.add_argument("--taxdump", help="Directory with NCBI names.dmp + nodes.dmp.")
        parser.add_argument("--llm", help="LLM spec provider:model, e.g. anthropic:claude-sonnet-5-5.")
        parser.add_argument("--no-cache", action="store_true",
                            help="Do not reuse cached Entrez answers.")
        parser.add_argument("--cache-hours", type=float, default=24.0,
                            help="Reuse cached Entrez answers younger than this (default 24; 0 = never expire).")
        args = parser.parse_args(argv)

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
        self.wcvp_path = args.wcvp or self.wcvp_path
        self.wfo_path = args.wfo or self.wfo_path
        for item in args.dwc:
            if "=" not in item:
                parser.error(f"--dwc expects NAME=PATH, got '{item}'")
            n, pth = item.split("=", 1)
            self.dwc_paths[n.strip()] = pth.strip()
        self.use_gbif = args.gbif
        self.group = args.group or self.group
        self.taxdump_dir = args.taxdump or self.taxdump_dir
        self.llm_spec = args.llm
        self.use_cache = not args.no_cache
        self.cache_hours = None if args.cache_hours == 0 else args.cache_hours

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
                 debug_queries=False, ai_validated_only=False, non_interactive=False,
                 checklists=None, gbif=None, group=None, taxdump=None, cache=None):
        os.makedirs(RESULTS_DIR, exist_ok=True)
        self.species_list = []
        self.dt_list      = []
        self.rows         = []
        self.llm          = None  # provider object, set by init_ai_client
        self.llm_log      = None
        self.exclude_infraspecific = exclude_infraspecific
        self.ai_unverified = ai_unverified
        self.verbose = verbose
        self.debug_queries = debug_queries
        self.ai_validated_only = ai_validated_only
        self.non_interactive = non_interactive
        self.checklists = list(checklists or [])   # Checklist objects (WCVP, WFO, ...)
        self.gbif = gbif                            # GBIFBackbone or None
        self.group = group                          # e.g. "Viridiplantae"
        self.taxdump = taxdump                      # NcbiTaxdump or None
        self.cache = cache                          # EntrezCache or None
        self._tax_name_cache = {}  # name -> [NcbiTaxon]
        self._taxon_info = {}      # species -> {"flags": [...], "taxon": NcbiTaxon}

        # Reliability / reproducibility bookkeeping.
        self._failed_cells = []     # list of (species, datatype) that returned NA
        self._matrix_records = []   # numeric per-species rows for the heatmap
        self._name_rows = []        # one row per (species, candidate name) classification
        self._comparisons = []      # per-species cross-authority synonymy comparison
        self._llm_records = []      # logged LLM calls used in this run
        self.run_started = time.strftime("%Y-%m-%d %H:%M:%S %z")

        # Be a bit more conservative if no API key is set.
        self.delay = 0.12 if getattr(Entrez, "api_key", None) else 0.4

    def init_ai_client(self, config_file, spec=None):
        """Initialize the LLM provider for synonym enrichment (any provider; see missmap_llm)."""
        try:
            self.llm = get_llm_client(config_file, spec)
            self.llm_log = LLMLog(os.path.join(RESULTS_DIR, "llm_calls.jsonl"))
            print(f"[INFO] LLM for synonym enrichment: {self.llm.name}:{self.llm.model} "
                  f"(prompt {PROMPT_VERSION}).")
        except Exception as e:
            print(f"[WARN] Could not initialize LLM client: {e}")
            print("[WARN] Continuing without AI-based synonym enrichment.")
            self.llm = None

    def lookup_synonyms_with_ai(self, species, manual_list):
        """Return scientific synonyms proposed by the LLM (call is logged)."""
        names, rec = lookup_synonyms_with_ai(species, manual_list, self.llm, log=self.llm_log)
        if rec is not None:
            self._llm_records.append({k: rec.get(k) for k in
                                      ("provider", "model", "mode", "species", "run_index",
                                       "prompt_version", "timestamp", "error")})
        return names

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
        if self.cache is not None:
            hit = self.cache.get("nuccore_count", term)
            if hit is not None:
                return int(hit)
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
                count = int(rec["Count"])
                if self.cache is not None:
                    self.cache.set("nuccore_count", term, count)  # only successes are cached
                return count
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

    # ------------------------------------------------------------------
    # NCBI Taxonomy lookups
    # ------------------------------------------------------------------
    def _entrez_retry(self, func, context, **kwargs):
        """Call an Entrez function with the same retry/backoff policy as counts."""
        last_err = None
        for attempt in range(MAX_RETRIES + 1):
            try:
                with func(**kwargs) as h:
                    data = h.read() if func is Entrez.efetch else Entrez.read(h)
                time.sleep(self.delay)
                return data
            except Exception as e:  # noqa: BLE001
                last_err = e
                if attempt < MAX_RETRIES and _is_retryable_error(e):
                    time.sleep(BACKOFF_BASE * (2 ** attempt))
                    continue
                break
        self.handle_errors(f"{context}: {last_err}")
        return None

    def _ncbi_candidates(self, name, filter_group=True):
        """
        Every NCBI taxon the name resolves to (as scientific name, synonym,
        equivalent name, ...). Filtered to the expected lineage group if set.
        """
        key = name.strip().lower()
        if not key:
            return []
        if key not in self._tax_name_cache:
            cands = None
            if self.taxdump is not None:
                cands = self.taxdump.candidates(name)
            else:
                cached = self.cache.get("taxonomy_candidates", key) if self.cache else None
                if cached is not None:
                    cands = [NcbiTaxon(**c) for c in cached]
                else:
                    res = self._entrez_retry(Entrez.esearch, f"taxonomy search '{name}'",
                                             db="taxonomy", term=f"\"{name}\"[All Names]", retmax=50)
                    if res is None:
                        cands = None   # lookup failed: do not cache as 'not found'
                    else:
                        ids = res.get("IdList", []) or []
                        cands = []
                        if ids:
                            xml = self._entrez_retry(Entrez.efetch, f"taxonomy fetch '{name}'",
                                                     db="taxonomy", id=",".join(ids), retmode="xml")
                            cands = parse_taxon_xml(xml) if xml else None
                        if cands is not None and self.cache is not None:
                            self.cache.set("taxonomy_candidates", key, [asdict(c) for c in cands])
            if cands is None:
                return []
            self._tax_name_cache[key] = cands
        cands = self._tax_name_cache[key]
        if filter_group and self.group:
            cands = [c for c in cands if c.in_group(self.group)]
        return cands

    def lookup_synonyms_manually(self, taxon):
        """
        Resolve the species to ONE NCBI TaxID (rank- and lineage-checked) and
        return (scientific synonyms, taxid). Ambiguous names return taxid None
        with flags explaining why, instead of silently taking the first hit.
        """
        cands = self._ncbi_candidates(taxon, filter_group=False)
        chosen, flags = choose_taxon(taxon, cands, expected_group=self.group)
        info = {"flags": flags, "taxon": chosen,
                "candidates": [f"{c.taxid}:{c.scientific_name}:{c.rank}" for c in cands]}
        self._taxon_info[taxon] = info
        if chosen is None:
            if "ambiguous" in flags or "not_in_expected_group" in flags:
                print(f"[WARN] '{taxon}': TaxID not resolved ({', '.join(flags)}). "
                      f"Candidates: {'; '.join(info['candidates']) or 'none'}")
            return [], None
        if flags and self.verbose:
            print(f"[INFO] '{taxon}': TaxID {chosen.taxid} ({', '.join(flags)})")
        return list(chosen.synonyms), chosen.taxid

    def _descendant_taxa(self, taxid):
        """Taxa below `taxid` (subspecies, varieties, ...) -- already in [Organism:exp]."""
        if not taxid:
            return []
        if self.taxdump is not None:
            out, stack = [], [str(taxid)]
            while stack:
                for ch in self.taxdump.children(stack.pop()):
                    t = self.taxdump.taxon(ch)
                    if t:
                        out.append(t)
                        stack.append(ch)
            return out
        key = str(taxid)
        cached = self.cache.get("taxonomy_subtree", key) if self.cache else None
        if cached is not None:
            return [NcbiTaxon(**c) for c in cached]
        res = self._entrez_retry(Entrez.esearch, f"subtree txid{taxid}",
                                 db="taxonomy", term=f"txid{taxid}[Subtree]", retmax=500)
        if res is None:
            return []
        ids = [i for i in (res.get("IdList", []) or []) if str(i) != key]
        out = []
        if ids:
            xml = self._entrez_retry(Entrez.efetch, f"subtree fetch txid{taxid}",
                                     db="taxonomy", id=",".join(ids), retmode="xml")
            out = parse_taxon_xml(xml) if xml else []
        if self.cache is not None:
            self.cache.set("taxonomy_subtree", key, [asdict(c) for c in out])
        return out

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
        """Return the in-group NCBI TaxIDs a name resolves to."""
        return [c.taxid for c in self._ncbi_candidates(name)]

    def _classify(self, species, taxid, names, descendants, source):
        """Classify candidate names; returns {name: info}. See missmap_classify.BUCKETS."""
        out = {}
        name_exists = self.gbif.name_exists if self.gbif is not None else None
        desc_ids = {d.taxid for d in descendants}
        for n in names:
            n = (n or "").strip()
            if not n or n.lower() == species.lower():
                continue
            info = classify_name(n, species, taxid, self._ncbi_candidates(n), self.checklists,
                                 name_exists=name_exists, descendant_taxids=desc_ids)
            info["source"] = source
            out[n] = info
        return out

    def _validate_ai_synonyms(self, species, taxid, ai_synonyms, descendants=()):
        """
        AI names usable for searching. A name qualifies if NCBI already files
        it under this TaxID, or if NCBI files it under ANOTHER TaxID and an
        authority confirms it as a synonym (that second group is where hidden
        records are). With --ai-unverified every scientific name qualifies.
        """
        if not ai_synonyms:
            return []
        if self.ai_unverified:
            return [n.strip() for n in ai_synonyms if isinstance(n, str) and n.strip()]
        if not taxid:
            print(f"[WARN] No TaxID for '{species}'; skipping AI synonyms to preserve accuracy.")
            return []
        cls = self._classify(species, taxid, ai_synonyms, descendants, "AI")
        return [n for n, i in cls.items() if i["bucket"] in ("same_taxid", "other_taxid_confirmed")]

    # ------------------------------------------------------------------
    # Main loop
    # ------------------------------------------------------------------
    def data_counts(self):
        """
        Per species:
          1. resolve one TaxID (rank/lineage checked) and NCBI's synonyms;
          2. collect each authority's synonymy and compare (Aim 1);
          3. ask the LLM, classify every AI and authority name (Aim 3);
          4. count records per datatype on the TaxID (primary), plus hidden
             records under other TaxIDs that authorities call synonyms, plus
             free-text over-counts (Aim 2).
        """
        for species in self.species_list:
            row = {"Species": species}

            # 1) TaxID + NCBI synonyms
            manual, taxid = self.lookup_synonyms_manually(species)
            tinfo = self._taxon_info.get(species, {})
            taxon = tinfo.get("taxon")
            row["TaxID"] = taxid if taxid else NA_LABEL
            row["TaxID Flags"] = ",".join(tinfo.get("flags", []))
            row["NCBI Common Names"] = ",".join(taxon.common_names) if taxon else ""

            # 2) Authorities
            auth_syns = {}
            if self.checklists:
                cmp_ = compare_synonymy(species, self.checklists, ncbi_synonyms=manual,
                                        include_infraspecific=not self.exclude_infraspecific)
                self._comparisons.append(cmp_)
                for cl in self.checklists:
                    row[f"{cl.name} Accepted Name"] = "; ".join(cmp_["accepted"].get(cl.name, []))
                    row[f"{cl.name} Synonyms"] = ",".join(cmp_["synonyms"].get(cl.name, []))
                row["Accepted Name Agrees Across Authorities"] = str(cmp_["accepted_agree"])
                auth_syns = authority_synonyms(species, self.checklists,
                                               include_infraspecific=not self.exclude_infraspecific)
            if self.gbif is not None:
                try:
                    _acc, gsyns = self.gbif.synonym_names(species)
                    for n in gsyns:
                        if self.exclude_infraspecific and self._is_infraspecific(n):
                            continue
                        auth_syns.setdefault(n, []).append("GBIF")
                except Exception as e:  # noqa: BLE001
                    print(f"[WARN] GBIF lookup failed for '{species}': {e}")
            ncbi_syn_keys = {s.lower() for s in manual}
            row["Authority Synonyms Missing From NCBI"] = ",".join(
                sorted(n for n in auth_syns if n.lower() not in ncbi_syn_keys))

            # 3) AI enrichment + classification of every candidate name
            ai_extra_raw = self.lookup_synonyms_with_ai(species, manual)
            ai_proposed = [s.strip() for s in (ai_extra_raw or []) if isinstance(s, str) and s.strip()]
            descendants = self._descendant_taxa(taxid) if (taxid and (ai_proposed or auth_syns)) else []
            ai_cls = self._classify(species, taxid, ai_proposed, descendants, "AI") if taxid else {}
            auth_cls = self._classify(species, taxid, list(auth_syns), descendants, "authority") if taxid else {}
            # A name that any authority source (incl. GBIF, which is not a loaded
            # checklist) lists as a synonym of this species is confirmed by it.
            for cls_ in (ai_cls, auth_cls):
                for n, i in cls_.items():
                    srcs = auth_syns.get(n)
                    if not srcs:
                        continue
                    i["confirmed_by"] = srcs
                    if i["bucket"] == "other_taxid_unconfirmed":
                        i["bucket"] = "other_taxid_confirmed"
                    elif i["bucket"] == "not_found":
                        i["bucket"] = "valid_not_in_ncbi"
            for n, i in list(ai_cls.items()) + list(auth_cls.items()):
                self._name_rows.append({
                    "Species": species, "TaxID": taxid, "Name": n, "Source": i["source"],
                    "Authority Sources": ",".join(auth_syns.get(n, [])),
                    "Bucket": i["bucket"], "NCBI TaxIDs": ",".join(i.get("ncbi_taxids", [])),
                    "Hidden TaxIDs": ",".join(i.get("hidden_taxids", [])),
                    "Authority Verdicts": json.dumps(i.get("authority", {}), sort_keys=True),
                    "GBIF Exact Match": i.get("gbif_exists", ""),
                })

            if self.ai_unverified:
                ai_validated = [n for n in ai_proposed if n.lower() != species.lower()]
            elif not taxid:
                if ai_proposed:
                    print(f"[WARN] No TaxID for '{species}'; skipping AI synonyms to preserve accuracy.")
                ai_validated = []
            else:
                ai_validated = [n for n, i in ai_cls.items()
                                if i["bucket"] in ("same_taxid", "other_taxid_confirmed")]
            ai_rejected = [s for s in ai_proposed if s not in ai_validated]
            for b in BUCKETS:
                row[f"AI: {b}"] = ",".join(n for n, i in ai_cls.items() if i["bucket"] == b)

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
            row["AI Validation Source"] = "NCBI Taxonomy" + "".join(
                f" + {cl.name}" for cl in self.checklists) + (" + GBIF" if self.gbif else "")
            row["synonyms"] = ",".join(syn_only_all)

            # Hidden TaxIDs: names filed under OTHER TaxIDs that authorities call synonyms.
            hidden_auth = sorted({t for i in auth_cls.values() if i["bucket"] == "other_taxid_confirmed"
                                  for t in i.get("hidden_taxids", [])})
            ai_hidden_buckets = {"other_taxid_confirmed"} | (
                {"other_taxid_unconfirmed"} if self.ai_unverified else set())
            hidden_ai = sorted({t for i in ai_cls.values() if i["bucket"] in ai_hidden_buckets
                                for t in i.get("hidden_taxids", [])})
            hidden_all = sorted(set(hidden_auth) | set(hidden_ai))
            row["Hidden TaxIDs (authority synonyms)"] = ",".join(hidden_auth)
            row["Hidden TaxIDs (AI synonyms)"] = ",".join(hidden_ai)

            # 4) Decide whether to use TaxID or fall back to name-based search
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
            row["NCBI Query (hidden records)"] = hidden_records_query(hidden_all, taxid) if use_taxid else ""

            # Name-based fallback query if TaxID is missing (used for main counts)
            organism_query = organism_query_expanded

            # numeric_counts holds the authoritative (TaxID-based when available)
            # per-datatype count for this species: int, or COUNT_NA on failure.
            numeric_counts = {}
            reconciled_counts = {}
            gained = []
            name_based_total_base = 0
            name_based_total_expanded = 0

            # 5) For each datatype, search nuccore
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

                # Hidden records for this datatype (only when hidden TaxIDs exist).
                if use_taxid and hidden_all:
                    hq = hidden_records_query(hidden_all, taxid, self.FILTERS[dt])
                    if self.debug_queries:
                        print(f"[QUERY] {species} | {dt} | hidden: {hq}")
                    hidden_dt = self._esearch_count(hq, context=f"{species} | {dt} | hidden")
                    if count is COUNT_NA or hidden_dt is COUNT_NA:
                        reconciled_counts[dt] = COUNT_NA
                    else:
                        reconciled_counts[dt] = count + hidden_dt
                        if count == 0 and hidden_dt > 0:
                            gained.append(f"{dt}(+{hidden_dt})")
                else:
                    reconciled_counts[dt] = count

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

            # 6) Aggregate across datatypes.
            #    The SUM double-counts: one record (e.g. a chloroplast genome)
            #    matches chloroplast, WGS, rbcL, matK, ycf1, ndhF at once. It is
            #    reported only as an overlapping tally, never as a record count.
            valid_counts = [c for c in numeric_counts.values() if c is not COUNT_NA]
            na_cells = sum(1 for c in numeric_counts.values() if c is COUNT_NA)
            sum_overlapping = sum(valid_counts)

            # 6b) DEDUPLICATED UNION: records matching AT LEAST ONE selected
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

            # 7) Total NCBI (no filters)
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

            # 8) Hidden records (recall loss) and free-text over-count (precision loss).
            def _count_or_zero(names_taxids, ctx):
                if not use_taxid:
                    return COUNT_NA      # no TaxID -> hidden records are not measurable
                if not names_taxids:
                    return 0
                return self._esearch_count(hidden_records_query(names_taxids, taxid), context=ctx)

            hidden_auth_n = _count_or_zero(hidden_auth, f"{species} | hidden (authority)")
            hidden_ai_n = _count_or_zero(hidden_ai, f"{species} | hidden (AI)")
            hidden_all_n = _count_or_zero(hidden_all, f"{species} | hidden (all)")
            row["Hidden Records (authority synonyms)"] = self._fmt(hidden_auth_n)
            row["Hidden Records (AI synonyms)"] = self._fmt(hidden_ai_n)
            row["Hidden Records (all)"] = self._fmt(hidden_all_n)
            if isinstance(total_all, int) and isinstance(hidden_all_n, int) and (total_all + hidden_all_n) > 0:
                row["Hidden Share of Reconciled Total (%)"] = (
                    f"{100.0 * hidden_all_n / (total_all + hidden_all_n):.3g}")
            else:
                row["Hidden Share of Reconciled Total (%)"] = NA_LABEL
            row["Markers Gained via Synonyms"] = ",".join(gained)

            if use_taxid:
                af = self._esearch_count(all_fields_query(species), context=f"{species} | all-fields")
                oc = self._esearch_count(overcount_query(species, taxid), context=f"{species} | over-count")
                noexp = self._esearch_count(f"txid{taxid}[Organism:noexp]", context=f"{species} | noexp")
                row["Free-text Name Search (All Fields)"] = self._fmt(af)
                row["Free-text Hits Not Under TaxID (over-count)"] = self._fmt(oc)
                row["TaxID Records Excl. Infraspecific (noexp)"] = self._fmt(noexp)
                if isinstance(af, int) and isinstance(oc, int) and af > 0:
                    row["Free-text Precision (%)"] = str(round(100.0 * (af - oc) / af, 2))
                else:
                    row["Free-text Precision (%)"] = NA_LABEL

            # 9) Lumping: names NCBI includes under this TaxID that an authority
            #    treats as separate species (records wrongly included).
            if use_taxid and self.checklists:
                from missmap_classify import authority_verdict, taxid_union
                lumped_syn = [n for n in manual
                              if "different_taxon" in authority_verdict(n, species, self.checklists).values()]
                lumped_children = [d for d in descendants
                                   if "different_taxon" in authority_verdict(
                                       d.scientific_name, species, self.checklists).values()]
                row["NCBI Synonyms Authorities Treat As Separate Species"] = ",".join(lumped_syn)
                row["NCBI Child Taxa Authorities Treat As Separate Species"] = ",".join(
                    f"{d.scientific_name} ({d.taxid})" for d in lumped_children)
                if lumped_children:
                    lc = self._esearch_count(taxid_union([d.taxid for d in lumped_children]),
                                             context=f"{species} | lumped children")
                    row["Records Under Those Child Taxa"] = self._fmt(lc)
                else:
                    row["Records Under Those Child Taxa"] = "0"

            # Record-keeping for the heatmap / coverage summary.
            self._matrix_records.append({
                "Species": species,
                "taxid": taxid,
                "counts": dict(numeric_counts),
                "counts_reconciled": dict(reconciled_counts),
                "union": union_count,
                "total_all": total_all,
                "hidden_all": hidden_all_n,
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

    def save_heatmap(self, ts, key="counts", tag=""):
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
                c = rec.get(key, rec["counts"]).get(dt, COUNT_NA)
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
        ax.set_title("MissMap: sequence data presence / absence"
                     + (" (after synonym reconciliation)" if key == "counts_reconciled" else ""),
                     fontsize=11)
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
        path1 = os.path.join(RESULTS_DIR, f"missmap_presence{tag}_{ts}.png")
        fig.savefig(path1, dpi=150, bbox_inches="tight")
        plt.close(fig)
        saved.append(path1)
        print(f"Saved presence/absence map to {path1}")
        if key == "counts_reconciled":
            return saved  # only the presence map is drawn for the reconciled view

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
        ai_model = f"{self.llm.name}:{self.llm.model}" if self.llm is not None else None
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
            "ai_prompt_version": PROMPT_VERSION if self.llm is not None else None,
            "ai_settings": getattr(self.llm, "settings", None),
            "ai_calls": self._llm_records,
            "expected_lineage_group": self.group,
            "authorities": {cl.name: cl.version for cl in self.checklists},
            "gbif_backbone_api": self.gbif is not None,
            "ncbi_taxonomy_source": (self.taxdump.version if self.taxdump is not None
                                     else "Entrez (live, run date)"),
            "taxid_resolution": {sp: {"flags": i.get("flags", []), "candidates": i.get("candidates", [])}
                                 for sp, i in self._taxon_info.items()},
            "entrez_cache": ({"path": self.cache.path, "max_age_hours": self.cache.max_age,
                              "hits": self.cache.hits, "misses": self.cache.misses}
                             if self.cache is not None else None),
            "failed_cells": self._failed_cells,
            "note": "Counts reflect NCBI nuccore on the run date and change daily.",
        }
        path = os.path.join(RESULTS_DIR, f"run_manifest_{ts}.json")
        with open(path, "w") as f:
            json.dump(manifest, f, indent=2)
        print(f"Saved run manifest to {path}")
        return path

    def save_name_classification(self, ts):
        """One row per (species, candidate name): bucket, NCBI TaxIDs, authority verdicts."""
        if not self._name_rows:
            return None
        path = os.path.join(RESULTS_DIR, f"name_classification_{ts}.csv")
        pd.DataFrame(self._name_rows).to_csv(path, index=False)
        counts = pd.DataFrame(self._name_rows).groupby(["Source", "Bucket"]).size()
        print("\n#==== Candidate names by source and bucket ====#\n")
        print(counts.to_string())
        print(f"\nSaved name classification to {path}")
        return path

    def save_comparisons(self, ts):
        """
        Aim 1 data: per-species cross-authority comparison (JSON) and a Darwin
        Core crosswalk (CSV) with each authority's usage for every name.
        """
        if not self._comparisons:
            return None
        jpath = os.path.join(RESULTS_DIR, f"synonymy_comparison_{ts}.json")
        with open(jpath, "w") as f:
            json.dump(self._comparisons, f, indent=1)
        rows = []
        for cmp_ in self._comparisons:
            sp = cmp_["species"]
            names = [sp] + [r["name"] for r in cmp_["names"]]
            for n in names:
                ncbi_ids = self._taxonomy_name_search(n)
                for cl in self.checklists:
                    for u in cl.find(n):
                        d = u.to_dwc()
                        d.update({"querySpecies": sp, "ncbiTaxIDs": ",".join(ncbi_ids),
                                  "matchMethod": "exact canonical name"})
                        rows.append(d)
                if not any(cl.find(n) for cl in self.checklists):
                    rows.append({"querySpecies": sp, "scientificName": n, "source": "(none)",
                                 "ncbiTaxIDs": ",".join(ncbi_ids), "matchMethod": "unmatched"})
        cpath = os.path.join(RESULTS_DIR, f"crosswalk_dwc_{ts}.csv")
        pd.DataFrame(rows).to_csv(cpath, index=False)
        print(f"Saved synonymy comparison to {jpath}")
        print(f"Saved Darwin Core crosswalk to {cpath}")
        return jpath, cpath

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
        print("\n#==== GenBank data availability matrix (summary; full matrix in CSV/XLSX) ====#\n")
        if not df.empty:
            summary_cols = [c for c in (
                ["Species", "TaxID", "TaxID Flags"] + list(self.dt_list) + [
                    "Records w/ >=1 Selected Marker (deduplicated)", "Total NCBI (no filters)",
                    "Hidden Records (all)", "Markers Gained via Synonyms",
                    "Free-text Hits Not Under TaxID (over-count)", "Free-text Precision (%)",
                    "Cells Failed (NA)"]) if c in df.columns]
            print(df[summary_cols].to_string(index=False))
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
        if any(r.get("counts_reconciled", r["counts"]) != r["counts"] for r in self._matrix_records):
            self.save_heatmap(ts, key="counts_reconciled", tag="_reconciled")
        self.save_name_classification(ts)
        self.save_comparisons(ts)
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

    # authorities
    checklists = []
    if config.wcvp_path:
        print(f"[INFO] Loading WCVP from {config.wcvp_path} ...")
        checklists.append(load_wcvp(config.wcvp_path))
    if config.wfo_path:
        print(f"[INFO] Loading WFO from {config.wfo_path} ...")
        checklists.append(load_wfo(config.wfo_path))
    for name, path in config.dwc_paths.items():
        print(f"[INFO] Loading {name} from {path} ...")
        checklists.append(load_dwc(path, name))
    for cl in checklists:
        print(f"[INFO] {cl.name}: {len(cl):,} name usages loaded.")
    if not checklists and not config.use_gbif:
        print("[INFO] No authority checklists loaded (--wcvp/--wfo/--dwc/--gbif); hidden-record "
              "and disagreement columns need at least one.")
    taxdump = NcbiTaxdump(config.taxdump_dir) if config.taxdump_dir else None
    cache = (EntrezCache(os.path.join(RESULTS_DIR, "entrez_cache.sqlite"), config.cache_hours)
             if config.use_cache else None)

    # build program & init AI client
    program = MissMap(
        exclude_infraspecific=config.exclude_infraspecific,
        ai_unverified=config.ai_unverified,
        verbose=config.verbose,
        debug_queries=config.debug_queries,
        ai_validated_only=config.ai_validated_only,
        non_interactive=config.non_interactive,
        checklists=checklists,
        gbif=GBIFBackbone() if config.use_gbif else None,
        group=config.group,
        taxdump=taxdump,
        cache=cache,
    )
    if not config.no_ai:
        program.init_ai_client(config.config_file, config.llm_spec)
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
