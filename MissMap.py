import os # file handling
import sys # error stuff
import time # to delay API calls
import argparse # accepting terminal input
import configparser # made it easier to accept config file
import pandas as pd # dataframe for the output matrix
from datetime import datetime # timestamp for files / logs
from Bio import Entrez # accessing NCBI Entrez API
from datetime import datetime # timestamp for files / logs
from urllib.error import HTTPError # handles API request errors
from xml.etree import ElementTree as ET # looks through XML on NCBI

# decide on using Entrez or OpenAI for synonym recognition - not working in my environment idek whatever
# from langchain_openai import AzureChatOpenAI # 
# from langchain.schema import SystemMessage, HumanMessage #
# from langchain_core.output_parsers import JsonOutputParser #

# SETUP 
# NCBI Log In --> Settings --> API Key Management --> Copy API Key
MISSMAP_INTRO = """
================================================================================#
  MissMap: A pipeline for visualizing sequence data availability in plant clades
================================================================================#

MissMap utilizes the Entrez API from NCBI's nucleotide database, thus recommends that you 
input your NCBI registered email address and API key in order to increase your API request 
rate from 3 requests per second to 10 requests per second. If you do not have an NCBI 
account/API, you can register on NCBI for free. 

Manual entry with email and API key: 
$ python MissMap.py -e yourncbiemail@address.com -k yourapikeylettersandnumbers -s "species1,species2" -dt "chloroplast","rbcL",...

Manual entry without email and API key using config file (recommended): 
$ python MissMap.py -s "species1" -dt "rbcL,matK,ITS"

Storing email & API key in your environment (config.txt):
# email yourncbiemail@address.com
# apikey yourapikeylettersandnumbers
# groq_api_key yourgroqapikey

================================================================================#
""" 

# only nucleotide data types plus common barcoding genes 
# add more based on what labmates say & papers
# trnH-psbA perhaps
DEFAULT_DATA_TYPES = [
    "chloroplast", "mitochondrion", "nuclear", "transcriptome", "RefSeq", 
    "WGS","rbcL", "matK", "trnL", "psbA-trnH", "atpF-atpH", "psbK-psbI",
    "rpoB", "rpoC1", "ycf1", "ndhF"
]

RESULTS_DIR = "missmap_results" 


# Configuration is an object to deal with various edgecases 
class Configuration:
    def __init__(self, config_file): # constructor
        self.config_file = "config.txt"
        self.email = None
        self.api_key = None
        self.groq_api_key = None  # optional Groq key
        self.species = None
        self.data_types = None
        self.disable_enhance = False
        
# reads & loads the config.txt file manually
    def load_config(self):
        try:
            with open(self.config_file, "r") as r: # attempts to open config if it exists
                for line in r:
                    line = line.strip() # remove whitespace
                    if not line or line.startswith("#"): # if there is no input, then skip
                        continue # goes back to the for loop (looking for next line)

                    split_line = line.split(maxsplit = 1) # only splits once, the list has 2 items in it (cfg_name, cfg_input)
                    # print(split_line) # prints as 2 lists containing 4 items total
                    if len(split_line) != 2: # if its not split into 2, 
                        continue # go to next line in file
                    cfg_requirement, cfg_input = split_line # splits list into separate variables

                    cfg_requirement = cfg_requirement.lower() 
                    if cfg_requirement in ("email", "entrez_email"): # if the email is found next to one of these keys
                        self.email = cfg_input # set the self object to user's email from config
                    elif cfg_requirement in ("apikey", "api_key", "entrez_api_key"): # if the api key is found next to one of these keys
                        self.api_key = cfg_input # set the self object to user's api key
                    elif cfg_requirement in ("groq_api_key", "groqapikey"): # config Groq key
                        self.groq_api_key = cfg_input
        except FileNotFoundError: # if the config file is deleted/not found
            pass # no config file, proceed 


# basically the same thing as config class, except it's manual and deals with species & data inputs too
    def terminal_inputs(self):
        parser = argparse.ArgumentParser(prog = "MissMap", description = "Pipeline for NCBI nucleotide data availability in plant clades.")

        parser.add_argument("-e", "--email",  help = "NCBI email", required = False)
        parser.add_argument("-k", "--apikey", help = "NCBI API key", required = False)
        parser.add_argument("-s", "--species", help = "Comma-separated list of species", required = False)
        parser.add_argument("-dt","--datatypes", help = "Comma-separated list of data types", required = False)
        parser.add_argument("--enhanceoff", dest = "disable_enhance", action = "store_true",
                            help = "Disable enhanced synonym recognition via Groq AI")
        
        args = parser.parse_args()

        self.disable_enhance = args.disable_enhance

        # dealing with manual email input 
        if args.email:
            self.email = args.email
        elif self.email:  # already from config.txt
            pass
        elif os.getenv("ENTREZ_EMAIL"): # lets user import from os in terminal for privacy reasons if user doesn't want config setup
            self.email = os.getenv("ENTREZ_EMAIL") # export ENTREZ_EMAIL = "youremail@address.com" at the top of the code
        else: # if no config or environment var is set up then just ask in terminal
            ncbi_email = input("NCBI Email Address (press Enter to skip): ").strip()
            if ncbi_email:
                self.email = ncbi_email

        # same manual process as above, just with API key
        if args.apikey:
            self.api_key = args.apikey
        elif self.api_key:  # config.txt
            pass
        elif os.getenv("ENTREZ_API_KEY") or os.getenv("API_KEY"):
            self.api_key = os.getenv("ENTREZ_API_KEY") or os.getenv("API_KEY") # export API_KEY = "yourapikey" at the top of the code
        else:
            user_api_key = input("API Key (press Enter to skip): ").strip()
            if user_api_key:
                self.api_key = user_api_key


        if self.email: # if we finally have an email then assign it to the Entrez API email object
            Entrez.email   = self.email
        if self.api_key: # same with api key
            Entrez.api_key = self.api_key


# asks for species and data inputs and runs missmap
class MissMap:
    FILTERS = {
        # organelle genomes
        "chloroplast":       "organelle:chloroplast",
        "mitochondrion":     "mitochondrion[Filter]",

        # transcripts
        "mRNA":          "mRNA[filter]",

        # nuclear genomes (excluding whole-genome shotgun)
        # "nuclear":           "biomol_genomic[PROP] NOT srcdb_refseq[PROP] NOT wgs[filter]",  # nuclear dna NOT from RefSeq or WGS

        # whole-genome shotgun
        "WGS":               "wgs[filter]",  # whole genome shotgun

        # # transcripts (excluding RefSeq)
        # "transcriptome":     "biomol_mrna[PROP] NOT srcdb_refseq[PROP]",  # mRNA NOT from RefSeq

        # # RefSeq transcripts only
        # "refseq_transcripts":"srcdb_refseq[PROP] AND biomol_mrna[PROP]",  # anything from RefSeq and mRNA

        # gene markers
        "rbcL":              "\"rbcL\"[Gene]",  # gene
        "matK":              "\"matK\"[Gene]",  # gene
        "ITS":               "internal transcribed spacer[All Fields]",  
        "trnL":              "\"trnL\"[Gene]",  # gene
        "psbA-trnH":         "\"psbA-trnH\"[All Fields]",   # gene but not in gene field?
        "atpF-atpH":         "\"atpF\"[Gene] AND \"atpH\"[Gene]",  # gene
        "psbK-psbI":         "\"psbK-psbI\"[All Fields]",  # gene but not in gene field?
        "rpoB":              "\"rpoB\"[Gene]",  # gene
        "rpoC1":             "\"rpoC1\"[Gene]",  # gene
        "ycf1":              "\"ycf1\"[Gene]",  # gene
        "ndhF":              "\"ndhF\"[Gene]",  # gene
    }

    # keys to sum into the total (facets only)
    FACET_KEYS = [
        "chloroplast", "mitochondrion", "mRNA", "WGS"
    ]

    # initializer / constructor for MissMap
    def __init__(self):
        os.makedirs(RESULTS_DIR, exist_ok = True) # makes results directory for potential csv file output and needed exist_ok so we dont get the FileExistsError
        self.species_list = [] # empty lists for inputs
        self.dt_list = []      # and data types
        self.rows = []         # and output rows
        self.disable_enhance = False
        self.groq_api_key = None

# asks user for species if it was not provided by argparse
    def ask_species(self):
        user_species_input = input("Enter species (comma-separated): ")
        for species in user_species_input.split(","): # for each comma separated species
            species.strip()
            if species:
                self.species_list.append(species) # if the species is provided, append it to the new list
        if not self.species_list: # if no species are provided then exit the program
            print("No species provided. Exiting MissMap.")
            sys.exit()

# asks user for data types if they don't want the whole output given by default
    def ask_data_types(self):
        print("Predefined data types:", ", ".join(self.FILTERS.keys())) # works to put filters automaticlaly for the print statement
        user_data_input = input("Enter data types or press Enter for all: ")
        for data_type in user_data_input.split(","): # split by comma ideally
            data_type.strip()
            if data_type:
                self.dt_list.append(data_type)
            else:
                self.dt_list = list(self.FILTERS.keys()) # defaults to all data types & genes available (harcoded though) 

    def handle_synonyms_with_groq(self, taxon):
        if not self.groq_api_key:
            return []
        url = "https://api.groq.com/openai/v1/chat/completions"
        model = "meta-llama/llama-4-scout-17b-16e-instruct"
        messages = [
            {"role": "system", "content": "You are a helpful assistant."},
            {"role": "user", "content": f"What are previous synonyms for {taxon}? Provide only a comma separated list as the response"}
        ]
        payload = {"model": model, "messages": messages, "temperature": 0.0}
        headers = {"Authorization": f"Bearer {self.groq_api_key}", "Content-Type": "application/json"}
        try:
            resp = requests.post(url, headers=headers, json=payload)
            resp.raise_for_status()
            text = resp.json()["choices"][0]["message"]["content"]
            return [syn.strip() for syn in text.split(",") if syn.strip()]
        except Exception:
            return []

    # actually counts the number of data per type
    def data_counts(self):
        for species in self.species_list:
            row = {"Species": species} # establishing matrix/dictionary format for the row list initialized earlier
            # total_count = 0  # iniitalizes count to 0 for a singular species before parsing

            #TaxonID‐based total ( this DOES include data types not listed in MissMap )
            handle = Entrez.esearch(db="taxonomy", term=f"{species}[Scientific Name]", retmax=1)
            result = Entrez.read(handle); ids = result.get("IdList", []); handle.close()
            if not ids:
                handle = Entrez.esearch(db="taxonomy", term=f"{species}[All Fields]", retmax=1)
                result = Entrez.read(handle); ids = result.get("IdList", []); handle.close()
            if not ids:
                print(f'Warning: No TaxID found for "{species}". Skipping.')
                continue

            taxid = ids[0]
            base  = f"txid{taxid}[Organism:exp]"
            handle = Entrez.esearch(db="nuccore", term=base, retmax=0)
            true_total = int(Entrez.read(handle)["Count"]); handle.close()
            row["Total sequences"] = str(true_total)

            # ——— Facet counts ———
            facet_sum = 0
            for dt in self.dt_list:
                term = f"{base} AND {self.FILTERS[dt]}" # term is the Entrez search term, this is what's use to search their database 
                try:
                    q = Entrez.esearch(db="nuccore", term=term, retmax=0) # searching ONLY in nucleotide database, retmax just asks for only one UID from XML output (default = 20)
                    cnt = int(Entrez.read(q)["Count"]); q.close()
                except HTTPError:
                    cnt = 0 # need more practice to identify errors and their causes- just go to skeleton function for now and keep count @ 0
                row[dt] = str(cnt) # sets up each row's respective column value as each query is made (counts for each species per datatype) 
                if dt in self.FACET_KEYS:
                    facet_sum += cnt
                time.sleep(0.12) # api rest time 
            row["Total Count of Data Found"] = str(facet_sum)  # adds total count of NCBI data counted for the data types provided ONLY

            # ——— Synonym‐only free‐text count ———
            ncbi_syns = self.lookup_synonyms_manually(species)
            groq_syns = [] if self.disable_enhance else self.handle_synonyms_with_groq(species)
            row["NCBI Taxonomy DB Synonyms"] = ",".join(ncbi_syns)
            row["Groq Enhanced Synonyms"]    = ",".join(groq_syns)

            # build the free‐text term
            clean = [s.split("(")[0].strip()
                    for s in ncbi_syns + groq_syns
                    if s and s.lower() != species.lower()]
            if clean:
                syn_term  = " OR ".join(f'"{s}"[Organism]' for s in clean)
                syn_q     = Entrez.esearch(db="nuccore", term=syn_term, retmax=0)
                syn_count = int(Entrez.read(syn_q)["Count"]); syn_q.close()
            else:
                syn_count = 0

            row["Synonym‐only count"] = str(syn_count)
            row["% Synonym‐only"]     = f"{syn_count/true_total*100:.1f}" if true_total else "0.0"

            self.rows.append(row) # appends the dictionary to the rows list 

# displays and saves the matrix of species to data availability in terminal and in a csv file
    def display_and_save_matrix(self):
        df = pd.DataFrame(self.rows) # build DataFrame
        print("\n#====================================== GenBank (nuccore ) data availability ======================================#\n")
        print(df.to_string(index = False)) # output table
        print("\n#==================================================================================================================#\n")
        time_completed = datetime.now().strftime("%Y%m%d-%H%M%S")
        filename = f"data_availability_matrix_{time_completed}.csv"
        output_file = os.path.join(RESULTS_DIR, filename)
        df.to_csv(output_file) # save CSV
        print(f"Saved results to a CSV file named {output_file} in the {RESULTS_DIR} directory.")

    def lookup_synonyms_manually(self, taxon):
        synonym_list = [] # will store each species' found synonym to display at the end
        try:
            query = Entrez.esearch(db = "taxonomy", term = f"{taxon.lower()}[Scientific Name]") # searches the scientific name using entrez in taxonomy database
            result = Entrez.read(query) # returns the query (its XML and put into a python dictionary)
            # print(result) # useful to check what data returns from this query 
            query.close()
            ids = result.get("IdList", []) # IdList gives all matching taxonomy IDs of the selection species, puts it in the ids list 
            if not ids: # if there are none, then search might work with all names
                query = Entrez.esearch(db = "taxonomy", term = f"{taxon.lower()}[All Names]")
                result = Entrez.read(query)
                query.close()
                ids = result.get("IdList", []) # get these results too if needed
            if not ids: # if there are no more results, return the synonym list for this species
                return synonym_list
            query = Entrez.efetch(db = "taxonomy", id = ids[0], retmode = "xml") # first id in the list is fetched in XML format
            xml = query.read() # reads the query
            query.close()

            root = ET.fromstring(xml) # uses ElementTree to parse the XML output
            for element in root.findall(".//OtherNames/Synonym"): # find all synonyms (element represents one <Synonym> tag)
                synonym = (element.text or "").strip() # 
                if synonym:
                    synonym_list.append(synonym) # append to final synonym list 
            for element in root.findall(".//OtherNames/GenbankCommonName"): # looks thru GenBankCommonName tag
                synonym = (element.text or "").strip()
                if synonym and synonym not in synonym_list:
                    synonym_list.append(synonym) # append to final synonym list 
        except Exception as e:
            self.handle_errors(e)

        return synonym_list # returns the list of synonyms for all taxon 

    def handle_errors(self, e):
        print(f"Error: {e}")


if __name__ == "__main__":
    print(MISSMAP_INTRO)

    config = Configuration("config.txt") # creates config object using Configuration class
    config.load_config() # loading in config file & its info if it exists
    config.terminal_inputs() # accepts terminal inputs for the configuration
    program = MissMap() # initialize MissMap
    program.disable_enhance = config.disable_enhance
    program.groq_api_key = config.groq_api_key

    if config.species:
        program.species_list = config.species # user passed -s, use it
    else:
        program.ask_species() # asks for species if not provided

    if config.data_types:
        program.dt_list = config.data_types # user passed -dt, use it
    else:
        program.ask_data_types() # asks for data types if not provided

    program.data_counts() # actually counts the number of data per type
    program.display_and_save_matrix() # displays and saves the matrix