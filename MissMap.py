from Bio import Entrez
import pandas as pd
import os
from datetime import datetime
import argparse


# SETUP 
# NCBI Log In --> Settings --> Scroll down to API Key Management --> Copy API Key
MISSMAP_INTRO = """
================================================================================#
  MissMap: A pipeline for visualizing sequence data availability in plant clades
================================================================================#


MissMap utilizes the Entrez API from NCBI GenBank, thus recommends that you input 
your NCBI registered email address and API key in order to increase your API request 
rate from 3 requests per second to 10 requests per second. If you do not have an NCBI 
account/API, you can register on NCBI for free, or simply skip this step. 

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


# obtains the email & API key so calls can be performed responsibly. 
def get_user_credentials():
    parser = argparse.ArgumentParser(prog = "MissMap", description = "A pipeline for visualizing sequence data availability in plant clades")
    parser.add_argument("-e", "--email", help = "The email you use to sign into NCBI with.")
    parser.add_argument("-k", "--apikey", help = "The API key can be obtained by going to NCBI Log In --> Settings --> Scroll down to API Key Management --> Copy API Key. **Note: The API key must be associated with the same email address you are using for your email input in order to unlock the higher request rate for the Entrez API.")
    
    
    
    args = parser.parse_args()
    email = (args.email or os.getenv("ENTREZ_EMAIL") or input("NCBI Email Address: "))
    api_key = (args.apikey or os.getenv("API_KEY") or input ("API Key: "))

    return email, api_key


if __name__ == "__main__":
    
    print(MISSMAP_INTRO, end = "\n\n")
    email, api_key = get_user_credentials()
    Entrez.email = email
    Entrez.api_key = api_key





# for later, when working with species and data type inputs
# parser.add_argument("-s", "--species")
# parser.add_argument("-dt", "--datatype")
