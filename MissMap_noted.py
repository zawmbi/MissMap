import os # file handling
import sys # error stuff
import time # to delay API calls
import math # math
import json # for formatting data
import argparse # accepting terminal input
import numpy as np # for
import configparser # made it easier to accept config file
import pandas as pd # dataframe for the output matrix
from Bio import Entrez # accessing NCBI Entrez API
from datetime import datetime # timestamp for files / logs
from urllib.error import HTTPError # handles API request errors
from xml.etree import ElementTree as ET # looks through XML on NCBI

# decide on using Entrez or OpenAI for synonym recognition 
# from langchain_openai import AzureChatOpenAI # 
# from langchain.schema import SystemMessage, HumanMessage #
# from langchain_core.output_parsers import JsonOutputParser #


# SETUP 
# NCBI Log In --> Settings --> API Key Management --> Copy API Key
MISSMAP_INTRO = """
================================================================================#
  MissMap: A pipeline for visualizing sequence data availability in plant clades
================================================================================#

MissMap utilizes the Entrez API from NCBI GenBank, thus recommends that you input 
your NCBI registered email address and API key in order to increase your API request 
rate from 3 requests per second to 10 requests per second. If you do not have an NCBI 
account/API, you can register on NCBI for free, or simply skip this step if you are not
performing a large number of requests. 

Manual entry with email and API key: 
$ python MissMap.py -e yourncbiemail@address.com -k yourapikeylettersandnumbers -s "species1,species2" -dt "chloroplast","rbcL",...

Manual entry without email and API key using config file (recommended): 
$ python MissMap.py -s "species1" -dt "rbcL,matK,ITS"

Storing email & API key in your environment (config.txt):
# email yourncbiemail@address.com
# apikey yourapikeylettersandnumbers

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

# def Config(): dont think this will be enough
# had to make config a whole object to deal with various edgecases 
class Configuration:
    def __init__(self, config_file): # constructor
        self.config_file = "config.txt"
        self.email = None
        self.api_key = None
        self.species = None
        self.data_types = None

        # for ai implementation 
        # self.model_name = model_name
        # self.llm = AzureChatOpenAI(
        #     deployment_name = 'gpt-4o',
        #     openai_api_version = data['AZURE_API_VERSION'],
        #     openai_api_key = data['AZURE_API_KEY_AZURE'],
        #     azure_endpoint = data['AZURE_API_BASE'],
        #     openai_organization = data['AZURE_ORGANIZATION'],
        #     temperature = 1.0,
        #     max_tokens = 256,
        #     top_p = 0,
        #     frequency_penalty = 0,
        #     presence_penalty = 0
        # )

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
        
        except FileNotFoundError: # if the config file is deleted/not found
            pass # idk if I should add an error message like "config.txt not found, proceeding with default api request value"

# uses configparser to load config.txt 
    # def load_config(self):
    #         config = configparser.ConfigParser()
    #         read = config.read(self.config_file)
    #         if not read:
    #             return # file wasn’t found or was empty

    #         self.email   = config.get("DEFAULT", "email", fallback = None)
    #         self.api_key = config.get("DEFAULT", "entrez_api_key", fallback = None)  


# basically the same thing as config class, except it's manual and deals with species & data inputs too
    def terminal_inputs(self):
        parser = argparse.ArgumentParser(prog = "MissMap", description = "Pipeline for GenBank nucleotide data availability in plant clades.")

        # script friendly
        parser.add_argument("-e", "--email",  help = "NCBI email", required = False)
        parser.add_argument("-k", "--apikey", help = "NCBI API key", required = False)
        parser.add_argument("-s", "--species", help = "Comma-separated list of species", required = False)
        parser.add_argument("-dt","--datatypes", help = "Comma-separated list of data types", required = False)
        args = parser.parse_args()

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



# might need to isolate it for testing but idk yet
    # def load_Entrez(self):
    #     if self.email:
    #         Entrez.email   = self.email
    #     if self.api_key:
    #         Entrez.api_key = self.api_key

# END OF CONFIG CLASS!




# actually asks for species and data inputs and runs missmap
class MissMap:

    # this is probably going to change a lot still - user feedback needed 
    FILTERS = {
        "chloroplast":    "chloroplast[Filter]",
        "mitochondrion":  "mitochondrion[Filter]", 
        "nuclear":        "biomol_genomic[PROP] NOT srcdb_refseq[PROP]", # nuclear dna NOT from RefSeq
        "transcriptome":  "biomol_mrna[PROP] NOT srcdb_refseq[PROP]", # mRNA NOT from RefSeq
        "RefSeq":         "srcdb_refseq[PROP]", # anything from RefSeq
        "WGS":            "wgs[Filter]", # whole genome shotgun 
        "rbcL":           "\"rbcL\"[Gene]", # gene
        "matK":           "\"matK\"[Gene]", # gene
        "ITS": "\"internal transcribed spacer\"[All Fields]", # not sure but somewhere said it was important
        "trnL":           "\"trnL\"[Gene]", # gene
        "psbA-trnH":      "\"psbA-trnH\"[All Fields]",  # gene but not in gene field?
        "atpF-atpH":      "\"atpF\"[Gene] AND \"atpH\"[Gene]", # gene
        "psbK-psbI":      "\"psbK-psbI\"[All Fields]", # gene but not in gene field?
        "rpoB":           "\"rpoB\"[Gene]", # gene
        "rpoC1":          "\"rpoC1\"[Gene]", # gene
        "ycf1":           "\"ycf1\"[Gene]", # gene
        "ndhF":           "\"ndhF\"[Gene]", # gene
    }



 # initializer / constructor for MissMap
    def __init__(self):
        os.makedirs(RESULTS_DIR, exist_ok = True) # makes results directory for potential csv file output and needed exist_ok so we dont get the FileExistsError
        self.species_list = [] #  empty lists to initialize for inputs and output rows
        self.dt_list = []
        self.rows = []



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


# actually counts the number of data per type
# admittedly used chatgpt a lot in this function
    def data_counts(self):
        for species in self.species_list: 
            row = {"Species": species} # establishing matrix/dictionary format for the row list initialized earlier
            total_count = 0  # iniitalizes count to 0 for a singular species before parsing

            for dt in self.dt_list: # for each data type in the accepted list
                term = f'"{species}"[Organism] AND {self.FILTERS[dt]}'
                print("DEBUG QUERY:", term)
                query = Entrez.esearch(db="nuccore", term=term, retmax=0)

                term = f'"{species}"[Organism] AND {self.FILTERS[dt]}' # term is the Entrez search term, this is what's use to search their database 
                query = None
                count = 0

                # ask for forgiveness not permission
                try:
                    # searching ONLY in nucleotide database, retmax just asks for only one UID from XML output (default = 20)
                    query = Entrez.esearch(db = "nuccore", term = term, retmax = 0) # which will be the count
                    record = Entrez.read(query)
                    count = int(record["Count"]) # obtains count of species found in nucleotide db
                    query.close()

                    total_count += count  # increments for the total even if there is no match 

                # need more practice to identify errors and their causes- just go to skeleton function for now and keep count @ 0
                except HTTPError as e:
                    self.handle_errors(e)
                    count = 0

                # this ALWAYS executes just in case no species are given or something goes wrong, it still closes the query request
                finally: 
                    if query != None:
                        query.close()
                
                row[dt] = str(count) # sets up each row's respective column value as each query is made (counts for each species per datatype) 
                time.sleep(0.12) # api rest

            row["Total Count of Data Found"] = str(total_count)  # adds total count of NCBI data counted for the data types provided ONLY
            row["synonyms"] = ",".join(self.lookup_synonyms_manually(species)) # comma joined synonym list for the synonym outputs
            # row["synonyms"] = ",".join(self.handle_synonyms_with_AI(species))

            # only happens if nothing returns- probably indicative of an error or misinput 
            if total_count == 0:
                print(f'"Warning: No sequence data found for "{species}" - check the spelling or try a known synonym."')

            self.rows.append(row)  # appends the dictionary to the rows list 




def data_counts(self):
    for species in self.species_list:
        row = {"Species": species}
        total = 0

        for dt in self.dt_list:
            base_term = f'"{species}"[Organism] AND {self.FILTERS[dt]}'
            query = None
            count = 0

            try:
                print("TRYING ▶", base_term)
                query = Entrez.esearch(db="nuccore", term=base_term, retmax=0)
                rec   = Entrez.read(query)
                count = int(rec.get("Count", 0))

            except HTTPError:
                # fallback: drop the qualifier entirely and just search by the raw dt string
                fallback = f'"{species}"[Organism] AND {dt}'
                print("FALLBACK ▶", fallback)
                try:
                    if query: query.close()
                    query = Entrez.esearch(db="nuccore", term=fallback, retmax=0)
                    rec   = Entrez.read(query)
                    count = int(rec.get("Count", 0))
                except:
                    count = 0

            finally:
                if query: 
                    try: query.close()
                    except: pass

            row[dt] = str(count)
            total  += count
            time.sleep(0.12)

        row["Total Count of Data Found"] = str(total)
        row["synonyms"] = ",".join(self.lookup_synonyms_manually(species))
        if total == 0:
            print(f'Warning: No data found for "{species}".')
        self.rows.append(row)



# displays and saves the matrix of species to data availability in terminal and in a csv file
    def display_and_save_matrix(self):
        df = pd.DataFrame(self.rows) # rows list from earlier is finally used
        print("\n#====================================== GenBank (nuccore ) data availability ======================================#\n")
        print(df.to_string(index = False)) # keep the index off unless you care about seeing numbers next to the species and matrix size
        print("\n#==================================================================================================================#\n")
        
        # need option to ask if user wants csv saved- maybe add to config function ? or is it too many inputs
        # also probably add date/time options so that the file doesn't overwrite itself
        time_completed = time.strftime("%Y%m%d-%H%M%S")
        # print(time_completed)
        filename = f"data_availability_matrix_{time_completed}.csv"
        output_file = os.path.join(RESULTS_DIR, filename)
        df.to_csv(output_file)
        print(f"Saved results to a CSV file named {output_file} in the {RESULTS_DIR} directory.")


# adding option to ask Azure to query for species synonyms or not is good for specificity that user's may want
# doesn't use OpenAI, but searchs through the taxonomy database using Entrez --> could be more accurate than AI 
# I did have to use chat for this function's search terms but it works really well since it's accurate to NCBI
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

        # chat helped a lot with this- 
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

        # sending it off to this incomplete function again            
        except Exception as e:
            self.handle_errors(e)


        return synonym_list # returns the list of synonyms for all taxon 


# # uses the same methods used in doomharvest.py- not sure if this is correct at all 
#     def handle_synonyms_with_AI(self, taxon):
#         parser = JsonOutputParser()


    def handle_errors(self, e):
        print(f"Error: {e}")

if __name__ == "__main__":
    print(MISSMAP_INTRO)

    config = Configuration("config.txt") # creates config object using Configuration class
    config.load_config() # loading in config file & its info if it exists
    config.terminal_inputs() # accepts terminal inputs for the configuration
    # config.load_Entrez()
    program = MissMap() # calls MissMap constructor 


    # program.ask_species() we don't need thi anymore since we have if args in ask_creds
    # program.ask_data_types()

    # if user passed -s, use it, otherwise prompt in the active terminal 
    if config.species:
        program.species_list = config.species
        print(program.species_list)
    else:
        program.ask_species()
        print(program.species_list)


    # same for data types
    if config.data_types:
        program.dt_list = config.data_types
    else:
        program.ask_data_types()
    program.data_counts()

    # displays the matrix and saves it in csv file in missmap_results directory
    program.display_and_save_matrix()


# notes
# wanna add a timer to see how long run is taking

# throw error for invalid key (it says) 

# if a family is submitted (like solanaceae) then the table outputs total for it, 
# in addition to the subsequent species under it (can be over 1000)

# adding option to ask Azure to query for species synonyms or not is good for specificity that user's may want

# might need edgecase for when user does not want to enter email/api key and it won't ask them again,
# or just replace email & api values in the config file with any invalid input so it doesnt ask 

# need to state that the default request value is being used when this happens

# need to add details like semblans has in the github repo decsription for ease of use
# dependencies !!!! important af
# need to add -c config_file option 
# include test data 