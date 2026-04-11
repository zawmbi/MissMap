# MissMap.py
# Integrated AI-based synonym enrichment in a single file (Groq version + optional subspecies exclusion)

import os                   # file handling
import sys                  # error stuff
import time                 # to delay API calls
import argparse             # accepting terminal input
import pandas as pd         # dataframe for the output matrix
from Bio import Entrez      # accessing NCBI Entrez API
from urllib.error import HTTPError  # handles API request errors
from xml.etree import ElementTree as ET  # parse XML from NCBI
import requests             # Groq API calls
import json                 # parse JSON from LLM

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
--verbose        Print per-datatype progress and count details.
--ai-unverified  Allow AI synonyms not validated in NCBI taxonomy.
--debug-queries  Print exact NCBI query strings used for counts.
--ai-validated-only  Use only NCBI-validated AI synonyms in searches (strict).

AI Disclaimer:
AI output can vary by model and input context. Synonym lists may differ across runs,
models, or prompt updates. Treat AI synonyms as suggestions to be verified.

================================================================================#
"""

RESULTS_DIR = "missmap_results"

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
        self.config_file       = config_file   # used by Terminal inputs and Entrez
        self.azure_config_file = config_file   # reuse same for LLM creds
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
        args = parser.parse_args()

        # Optional override of config file
        if args.configfile:
            self.config_file       = args.configfile
            self.azure_config_file = args.configfile
            # reload config with new file path
            self.load_config()

        # Email
        if args.email:
            self.email = args.email
        elif not self.email and os.getenv("ENTREZ_EMAIL"):
            self.email = os.getenv("ENTREZ_EMAIL")
        elif not self.email:
            e = input("NCBI Email Address (press Enter to skip): ").strip()
            if e:
                self.email = e

        # API key
        if args.apikey:
            self.api_key = args.apikey
        elif not self.api_key and os.getenv("ENTREZ_API_KEY"):
            self.api_key = os.getenv("ENTREZ_API_KEY")
        elif not self.api_key:
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

    def __init__(self, exclude_infraspecific=False, ai_unverified=False, verbose=False, debug_queries=False, ai_validated_only=False):
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
        self._tax_name_cache = {}  # cache taxonomy lookups by name

        # Be a bit more conservative if no API key is set.
        self.delay = 0.12 if getattr(Entrez, "api_key", None) else 0.4

    def init_ai_client(self, azure_config_file):
        """Initialize Groq LLM config for synonym enrichment (reusing config file path)."""
        try:
            self.llm = get_llm_client(azure_config_file)
            print("[INFO] Groq client config initialized for synonym enrichment.")
        except Exception as e:
            print(f"[WARN] Could not initialize Groq client: {e}")
            print("[WARN] Continuing without AI-based synonym enrichment.")
            self.llm = None

    def lookup_synonyms_with_ai(self, species, manual_list):
        """Return extra synonyms using the LLM."""
        return lookup_synonyms_with_ai(species, manual_list, self.llm)

    def ask_species(self):
        inp = input("Enter species (comma-separated): ")
        for sp in inp.split(","):
            sp = sp.strip()
            if sp:
                self.species_list.append(sp)
        if not self.species_list:
            print("No species provided. Exiting.")
            sys.exit(1)

    def ask_data_types(self):
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


    def _filter_infraspecific(self, names, species):
        """
        Optionally drop infraspecific names (var., subsp.) from the list.
        Keeps the original species name regardless.
        """
        if not self.exclude_infraspecific:
            return names

        filtered = []
        for n in names:
            low = n.lower()
            # Always keep the canonical species name
            if low == species.lower():
                filtered.append(n)
                continue
            # Drop infraspecific ranks if they look like 'var.' or 'subsp.'
            if " var. " in low or " subsp. " in low:
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

            total_count = 0
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

                count = 0
                try:
                    with Entrez.esearch(
                        db="nuccore",
                        term=term,
                        retmax=0
                    ) as q:
                        rec = Entrez.read(q)
                    count = int(rec.get("Count", 0))
                    total_count += count
                except HTTPError as e:
                    self.handle_errors(
                        f"HTTP error for species '{species}', datatype '{dt}': {e}"
                    )
                except Exception as e:
                    self.handle_errors(
                        f"Error for species '{species}', datatype '{dt}': {e}"
                    )

                row[dt] = str(count)
                if self.verbose:
                    print(f"[INFO] {species} | {dt} | TaxID/primary count: {count}")
                time.sleep(self.delay)

                # Name-based comparison counts (NCBI synonyms vs NCBI+AI)
                try:
                    with Entrez.esearch(
                        db="nuccore",
                        term=f"{organism_query_base} AND {self.FILTERS[dt]}",
                        retmax=0
                    ) as q:
                        rec = Entrez.read(q)
                    base_count = int(rec.get("Count", 0))
                except Exception as e:
                    self.handle_errors(
                        f"Error for name-based NCBI synonyms '{species}', datatype '{dt}': {e}"
                    )
                    base_count = 0
                time.sleep(self.delay)

                try:
                    with Entrez.esearch(
                        db="nuccore",
                        term=f"{organism_query_expanded} AND {self.FILTERS[dt]}",
                        retmax=0
                    ) as q:
                        rec = Entrez.read(q)
                    expanded_count = int(rec.get("Count", 0))
                except Exception as e:
                    self.handle_errors(
                        f"Error for name-based NCBI+AI '{species}', datatype '{dt}': {e}"
                    )
                    expanded_count = 0
                time.sleep(self.delay)

                name_based_total_base += base_count
                name_based_total_expanded += expanded_count
                if self.verbose:
                    print(
                        f"[INFO] {species} | {dt} | Name-based NCBI: {base_count} | NCBI+AI: {expanded_count}"
                    )
                if self.debug_queries:
                    print(
                        f"[QUERY] {species} | {dt} | name-based NCBI: {organism_query_base} AND {self.FILTERS[dt]}"
                    )
                    print(
                        f"[QUERY] {species} | {dt} | name-based NCBI+AI: {organism_query_expanded} AND {self.FILTERS[dt]}"
                    )

            # 5) Sum across datatypes
            row["Total Count of Data Found"] = str(total_count)
            row["Name-based Total (NCBI Synonyms)"] = str(name_based_total_base)
            row["Name-based Total (NCBI+AI)"] = str(name_based_total_expanded)
            row["Name-based Total Gain (AI - NCBI)"] = str(
                name_based_total_expanded - name_based_total_base
            )

            # 6) Total NCBI (no filters)
            total_all = 0
            total_all_base = 0
            total_all_expanded = 0
            try:
                if use_taxid:
                    base_term = f"txid{taxid}[Organism:exp]"
                else:
                    base_term = organism_query

                with Entrez.esearch(
                    db="nuccore",
                    term=base_term,
                    retmax=0
                ) as q:
                    rec = Entrez.read(q)
                total_all = int(rec.get("Count", 0))
            except Exception as e:
                self.handle_errors(
                    f"Error getting total unfiltered count for '{species}': {e}"
                )

            row["Total NCBI (no filters)"] = str(total_all)
            if self.debug_queries:
                print(f"[QUERY] {species} | no-filter primary: {base_term}")

            # Name-based unfiltered comparison (NCBI synonyms vs NCBI+AI)
            try:
                with Entrez.esearch(
                    db="nuccore",
                    term=organism_query_base,
                    retmax=0
                ) as q:
                    rec = Entrez.read(q)
                total_all_base = int(rec.get("Count", 0))
            except Exception as e:
                self.handle_errors(
                    f"Error getting name-based total (NCBI synonyms) for '{species}': {e}"
                )
                total_all_base = 0
            time.sleep(self.delay)

            try:
                with Entrez.esearch(
                    db="nuccore",
                    term=organism_query_expanded,
                    retmax=0
                ) as q:
                    rec = Entrez.read(q)
                total_all_expanded = int(rec.get("Count", 0))
            except Exception as e:
                self.handle_errors(
                    f"Error getting name-based total (NCBI+AI) for '{species}': {e}"
                )
                total_all_expanded = 0
            time.sleep(self.delay)

            row["Name-based NCBI (no filters)"] = str(total_all_base)
            row["Name-based NCBI+AI (no filters)"] = str(total_all_expanded)
            row["Name-based Gain (no filters)"] = str(
                total_all_expanded - total_all_base
            )
            if self.debug_queries:
                print(f"[QUERY] {species} | no-filter name-based NCBI: {organism_query_base}")
                print(f"[QUERY] {species} | no-filter name-based NCBI+AI: {organism_query_expanded}")

            if total_count == 0:
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


    def display_and_save_matrix(self):
        df = pd.DataFrame(self.rows)
        print("\n#==== GenBank data availability matrix ====#\n")
        if not df.empty:
            print(df.to_string(index=False))
        else:
            print("[INFO] No rows to display.")
        print("\n#========================================#\n")
        ts = time.strftime("%Y%m%d-%H%M%S")
        fn_csv = f"data_availability_matrix_{ts}.csv"
        out_csv = os.path.join(RESULTS_DIR, fn_csv)
        df.to_csv(out_csv, index=False)
        print(f"Saved results to {out_csv}")

        # Also save Excel with auto-fit column widths
        fn_xlsx = f"data_availability_matrix_{ts}.xlsx"
        out_xlsx = os.path.join(RESULTS_DIR, fn_xlsx)
        try:
            with pd.ExcelWriter(out_xlsx, engine="xlsxwriter") as writer:
                df.to_excel(writer, index=False, sheet_name="MissMap")
                worksheet = writer.sheets["MissMap"]
                for i, col in enumerate(df.columns):
                    # set width based on max of header/data
                    series = df[col].astype(str)
                    max_len = max([len(str(col))] + [len(v) for v in series]) + 2
                    worksheet.set_column(i, i, min(max_len, 120))
            print(f"Saved results to {out_xlsx}")
        except Exception as e:
            print(f"[WARN] Could not write Excel file: {e}")

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
        ai_validated_only=config.ai_validated_only
    )
    if not config.no_ai:
        program.init_ai_client(config.azure_config_file)
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
