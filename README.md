# MissMap
**Author:** Linda Mansour, *Walker Lab, UIC Department of Biological Sciences*

MissMap is a Python pipeline for visualizing GenBank nucleotide data availability across plant species.  It queries NCBI’s Entrez API to count sequences by data type (e.g. transcriptome, nuclear, barcoding genes), assembles the results into a matrix, and prints a table to the terminal and a creates a CSV saved in the RESULTS_DIR folder.  It also enriches each identifies the common species synonyms provided by the NCBI taxonomy database.

## Dependencies

MissMap requires:

- Python 3.8+  
- [Biopython](https://biopython.org/)  
- [pandas](https://pandas.pydata.org/)  
- [argparse](https://docs.python.org/3/library/argparse.html) (standard library)

  ### Currently Optional 
- [langchain](https://github.com/langchain-ai/langchain) & [langchain_core](https://github.com/langchain-ai/langchain-core) — **for future AI implementation**  
- [openai](https://pypi.org/project/openai/) or [azure-ai-openai](https://pypi.org/project/azure-ai-openai/) — **for future AI implementation**  

On Ubuntu you can install most prerequisites with:

```bash
sudo apt update
sudo apt install python3-pip python3-biopython
pip3 install pandas langchain langchain_core openai azure-ai-openai
```

## Installation

To install MissMap, clone the repository or copy and paste this command in the console while the current directory is the folder you wish to install MissMap to.

```bash
git clone https://github.com/zawmbi/MissMap/
```


