# MissMap
A python pipeline for visualizing sequence data availability in plant clades.

## Config
MissMap reads credentials and settings from `config.txt`, environment variables, or a `.env` file.
Supported keys: `ENTREZ_EMAIL`, `ENTREZ_API_KEY`, `GROQ_API_KEY`, `GROQ_MODEL`, plus optional
`species` and `datatypes` in `config.txt`.

## AI Disclaimer
AI output can vary by model and input context. Synonym lists may differ across runs,
models, or prompt updates. Treat AI synonyms as suggestions to be verified.

## Columns
Species: input species name.
NCBI Synonyms: synonyms from NCBI Taxonomy OtherNames.
AI Proposed Synonyms: AI output list (unverified).
AI-validated Synonyms: AI synonyms that resolve to the same NCBI TaxID as the species.
AI Rejected Synonyms: AI synonyms that do not validate to the same TaxID.
AI Validation Source: source used for validation.
synonyms: NCBI Synonyms plus AI Proposed Synonyms used in name-based search (unless --ai-validated-only).
NCBI Query (no filters, primary): TaxID query if TaxID exists, else name-based NCBI+AI query.
NCBI Query (no filters, NCBI synonyms): name-based query using species + NCBI synonyms.
NCBI Query (no filters, NCBI+AI): name-based query using species + NCBI + AI synonyms.
chloroplast, mitochondrion, nuclear, transcriptome, RefSeq, WGS, rbcL, matK, ITS, trnL, psbA-trnH, atpF-atpH, psbK-psbI, rpoB, rpoC1, ycf1, ndhF: counts from nuccore using the datatype filter and the primary query (TaxID if available).
Total Count of Data Found: sum of datatype counts above.
Name-based Total (NCBI Synonyms): sum of datatype counts using name-based NCBI synonyms query.
Name-based Total (NCBI+AI): sum of datatype counts using name-based NCBI+AI query.
Name-based Total Gain (AI - NCBI): Name-based Total (NCBI+AI) minus Name-based Total (NCBI Synonyms).
Total NCBI (no filters): unfiltered count using primary query.
Name-based NCBI (no filters): unfiltered count using name-based NCBI synonyms query.
Name-based NCBI+AI (no filters): unfiltered count using name-based NCBI+AI query.
Name-based Gain (no filters): Name-based NCBI+AI minus Name-based NCBI.
