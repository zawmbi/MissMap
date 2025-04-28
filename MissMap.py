from Bio import Entrez
import time
import pandas as pd

# Must have an NCBI account to obtain your API key ( theres also an email option ) 
# Log In --> Settings --> Scroll down to API Key Management --> Copy API Key
Entrez.api_key = "f96b1da3333784c245eb48fbea93966eeb08"

# NCBI can contact you if there is a problem
Entrez.email = "ljmansour02@gmail.com" 

# to be filled with species inputted by the user
species_list = ["Arabidopsis thaliana",
                "Zea mays",
                "Oryza sativa",
                "Lemna gibba"
                ] # arabidopsis, maize, rice, duckweed

# depending on the type that the user wants to search for, maybe add an option to add all types;
# ok this is super messy, theres a bunch of types so I suppose an "all" option would not work that well
data_types = {
    "chloroplast": "chloroplast",
    "mitochondrial": "mitochondrion",
    "nuclear": "nuclear",
    "transcriptome": "transcriptome",
    "barcode": "barcode",
} # barcode genes will all be different
#  splitting up search would be optimal 
 
stream = Entrez.einfo(db = "nucleotide") # tells us the fields we have available to search from- cannot use all of them!!
# it's just for listing
# consult with joe to see what we should use
record = Entrez.read(stream)
stream.close()

for field in record["DbInfo"]["FieldList"]:
    print(f"{field['Name']}: {field['Description']}")