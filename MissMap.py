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

Storing email & API key and Azure settings in your environment (config.txt):
# email yourncbiemail@address.com
# apikey yourapikeylettersandnumbers
# AZURE_API_KEY_AZURE yourazureapikey
# AZURE_API_VERSION 2024-02-01
# AZURE_API_BASE https://your-azure-endpoint.openai.azure.com/

You can also optionally set:
# species Arabidopsis thaliana,Zea mays
# datatypes chloroplast,mitochondrion,rbcL

================================================================================#
"""

RESULTS_DIR = "missmap_results"

# ======================================================================
# === AI-BASED SYNONYM ENRICHMENT HELPERS (Groq) =======================
# ======================================================================

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
    model = cfg.get("GROQ_MODEL", "meta-llama/llama-4-scout-17b-16e-instruct")

    if not api_key:
        raise KeyError(
            "Missing GROQ_API_KEY. Set it in config.txt or as an environment variable."
        )

    return {
        "api_key": api_key,
        "model": model,
    }


# Template for enriching synonyms via AI (unchanged)
SYN_TEMPLATE = (
    "You are a botanical taxonomy assistant. Given a species name and a list of known synonyms, "
    "return any additional synonyms, common names, or spelling variants in JSON format only.\n"
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
        "temperature": 0.2,
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
        "WGS":           "wgs[filter]",
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

    def __init__(self, exclude_infraspecific=False):
        os.makedirs(RESULTS_DIR, exist_ok=True)
        self.species_list = []
        self.dt_list      = []
        self.rows         = []
        self.llm          = None  # will be set by init_ai_client
        self.exclude_infraspecific = exclude_infraspecific

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
                return syns

            with Entrez.efetch(db="taxonomy", id=ids[0], retmode="xml") as q:
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
        return syns

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

    def data_counts(self):
        """Query NCBI for each species × datatype using species + synonyms, and collect counts + synonyms."""
        for species in self.species_list:
            row = {"Species": species}

            # 1) Get manual + AI synonyms FIRST
            manual = self.lookup_synonyms_manually(species)
            ai_extra = self.lookup_synonyms_with_ai(species, manual)
            # include the original species name in the set
            all_names = [species] + (manual or []) + (ai_extra or [])
            # deduplicate and strip
            all_names = sorted(set(n.strip() for n in all_names if n.strip()))

            # optionally drop subspecies/varieties
            all_names = self._filter_infraspecific(all_names, species)

            # save synonyms (without duplicating the original species name)
            syn_only = [n for n in all_names if n.lower() != species.lower()]
            row["synonyms"] = ",".join(syn_only)

            # 2) Build OR query over all names for the Organism field
            if all_names:
                organism_query = " OR ".join(
                    [f'"{name}"[Organism]' for name in all_names]
                )
            else:
                # fallback to just the original species if something weird happened
                organism_query = f'"{species}"[Organism]'

            organism_query = f"({organism_query})"

            total_count = 0

            # 3) For each datatype, search using species + synonyms
            for dt in self.dt_list:
                term = f'{organism_query} AND {self.FILTERS[dt]}'
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
                time.sleep(self.delay)

            # 4) Sum across datatypes
            row["Total Count of Data Found"] = str(total_count)

            # 5) Total NCBI (no filters), also using species + synonyms
            total_all = 0
            try:
                with Entrez.esearch(
                    db="nuccore",
                    term=organism_query,
                    retmax=0
                ) as q:
                    rec = Entrez.read(q)
                total_all = int(rec.get("Count", 0))
            except Exception as e:
                self.handle_errors(
                    f"Error getting total unfiltered count for '{species}' (with synonyms): {e}"
                )

            row["Total NCBI (no filters)"] = str(total_all)

            if total_count == 0:
                print(
                    f'[WARN] No sequence data found for "{species}" (including synonyms) – '
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
        fn = f"data_availability_matrix_{ts}.csv"
        out = os.path.join(RESULTS_DIR, fn)
        df.to_csv(out, index=False)
        print(f"Saved results to {out}")

# ======================================================================
# === Main execution ===================================================
# ======================================================================
if __name__ == "__main__":
    print(MISSMAP_INTRO)

    # parse config & credentials
    config = Configuration("config.txt")
    config.load_config()
    config.terminal_inputs()

    # build program & init AI client
    program = MissMap(exclude_infraspecific=config.exclude_infraspecific)
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
