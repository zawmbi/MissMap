from Bio import Entrez
import pandas as pd
import os
from datetime import date
import argparse 
import re
import urllib.error import HTTPError


# SETUP 
# NCBI Log In --> Settings --> Scroll down to API Key Management --> Copy API Key
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
$ python MissMap.py -e yourNCBIemail@address.com -k 1a2b3c4d5e6f7g8h9i -s1 "species name" -dt1 "data type"

Manual entry without email and API key: 
$ python MissMap.py -s1 "species name" -dt1 "data type" -s2 "species2 name" -dt2 "data type2"

Storing email & API key in your environment using bash or zsh (No need for re-entry):
# gotta figure this out 

If you are storing the email and API information, restart your shell after initialization to save.
 
More help can be found using -h or --help. 

================================================================================#

"""
DEFAULT_DATA_TYPES = [""] # unsure about how many / which ones to add
RESULTS_DIR = "missmap_results"
DEFAULT_DATA_TYPES = ["chloroplast", "mitochondrion", "nuclear", 
                      "transcriptome", "RefSeq", "WGS", "EST", 
                      "GSS", "HTG", "TSA"]



# def read_config(config_file):
#     data = {} # initializes empty dictionary
#     f = open(config_file, "r") # open & read file
    
#     for i in f: # for each line in the file
#         i = i.strip() # remove unecessary whitespace
#         if len(i) > 1: # if the length of the line has at least one char
#             split_lines = i.split(" ") #  split the line by spaces
#             data[split_lines[0]] = split_lines[1] # the dict 
            
#     f.close()
#     return data 

# obtains the email & API key so calls can be performed responsibly. 
def get_user_credentials():
    parser = argparse.ArgumentParser(prog = "MissMap", description = "A pipeline for visualizing sequence data availability in plant clades")
    parser.add_argument("-e", "--email", help = "The email you use to sign into NCBI with.", required = False)
    parser.add_argument("-k", "--apikey", help = "The API key can be obtained by going to NCBI Log In --> Settings --> Scroll down to API Key Management --> Copy API Key. **Note: The API key must be associated with the same email address you are using for your email input in order to unlock the higher request rate for the Entrez API.", required = False)
    
    args = parser.parse_args()
    email = (args.email or os.getenv("ENTREZ_EMAIL") or input("NCBI Email Address: "))
    api_key = (args.apikey or os.getenv("API_KEY") or input ("API Key: "))
    

    return email, api_key


if __name__ == "__main__":
    
    print(MISSMAP_INTRO, end = "\n\n")
    email, api_key = get_user_credentials()
    Entrez.email = email
    Entrez.api_key = api_key
    
    import time
    from xml.etree import ElementTree as ET

    # ─── collect species interactively ──────────────────────────────────────────
    species_input = input(
        "\nEnter one or more species names separated by commas: "
    )
    species_list = [s.strip() for s in species_input.split(",") if s.strip()]
    if not species_list:
        print("No species supplied → exiting."); exit()

    # optional custom data-types, else use the defaults already defined
    dt_input = input(
        
        f"Enter data types from {DEFAULT_DATA_TYPES} separated by commas "
        "(or press Enter to use all): " # definitely change this eventually, it's not cute
    )
    dt_list = [d.strip() for d in dt_input.split(",") if d.strip()] or DEFAULT_DATA_TYPES


# This is the most confusing part so far
    # Entrez query fragments for each MissMap data-type
    
    # implement synonym search
    FILTERS = {
        "chloroplast"  : "organelle:chloroplast",
        "mitochondrion": "organelle:mitochondrion",
        "nuclear"      : "(organelle:noexist) AND biomol_genomic[PROP]",
        "transcriptome": "biomol_mrna[PROP] NOT srcdb_refseq[PROP]",
        "RefSeq"       : "srcdb_refseq[PROP]",
        "WGS"          : "wgs[filter]",
        "EST"          : "est[filter]",
        "GSS"          : "gss[filter]",
        "HTG"          : "htgs[filter]",
        "TSA"          : "tsa[filter]",
    }

    rows = []
    for sp in species_list:
        row = {"species": sp}
        for dt in dt_list:
            term = f'"{sp}"[Organism]'
            if dt in FILTERS:
                term += f" AND {FILTERS[dt]}"
                # entrex esearch function using nuccore database
            h = Entrez.esearch(db = "nuccore", term = term, rettype = "count")
            response = Entrez.read(h)
            count = int(response.get("Count", 0)) if isinstance(response, dict) else 0
            h.close()
            row[dt] = str(count)
            time.sleep(0.12)  # keeps under 10 req/s with API key
        rows.append(row)

    # matrix printing (maybe make its own function ?)
    df = pd.DataFrame(rows)
    print("\n#===================== GenBank data availability ======================#\n")
    print(df.to_string(index = False))
    print("\n#=======================================================================#\n")


# for later, when working with species and data type inputs
# parser.add_argument("-s", "--species")
# parser.add_argument("-dt", "--datatype")

# API
# f96b1da3333784c245eb48fbea93966eeb08

# notes
# contact emily mctavish

# python MissMap.py -e ljmansour02@gmail.com -k f96b1da3333784c245eb48fbea93966eeb08

# throw error for invalid key (it says) 

# if a family is submitted (like solanaceae) then the table outputs total for it, 
# in addition to the subsequent species under it (can be over 1000)

# use Steven's lab query for inspo to find synonyms (like what else is solanaceae called?)
# incorporate chatgpt for this potentially (or Open Azure)
# https://github.com/FePhyFoFum/doomharvest/blob/main/doomharvest.py

# config file probably necessary

# eventually query over certain genes & barcoding genes

# adding option to ask Azure to query for species synonyms or not is good for specificity that user's may want
