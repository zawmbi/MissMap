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
**Species:** the inputted species name(s).

**NCBI Synonyms:** synonyms from NCBI Taxonomy Database.

**AI Proposed Synonyms:** AI output list (unverified).

**AI-validated Synonyms:** AI synonyms that resolve to the same NCBI TaxID as the species.

**AI Rejected Synonyms:** AI synonyms that do not validate to the same TaxID.

**AI Validation Source:** source used for validation.

**Synonyms:** NCBI Synonyms plus AI Proposed Synonyms used in name-based search (unless --ai-validated-only).

**NCBI Query (no filters, primary):** TaxID query if TaxID exists, else name-based NCBI+AI query.

**NCBI Query (no filters, NCBI synonyms):** name-based query using species + NCBI synonyms.

**NCBI Query (no filters, NCBI+AI):** name-based query using species + NCBI + AI synonyms.

**chloroplast, mitochondrion, nuclear, transcriptome, RefSeq, WGS, rbcL, matK, ITS, trnL, psbA-trnH, atpF-atpH, psbK-psbI, rpoB, rpoC1, ycf1, ndhF:** counts from nuccore using the datatype filter and the primary query (TaxID if available). A cell shows `NA` (not `0`) when the query failed after retries — a genuine "no data" and a network error are never conflated.

**Records w/ >=1 Selected Marker (deduplicated):** the true record count of sequences matching *at least one* selected datatype, obtained by OR-combining all filters in a single query so NCBI deduplicates. This is the headline "how much data" number and is always ≤ *Total NCBI (no filters)*.

**Sum of Marker Counts (overlapping, not a record count):** the naive sum of the per-datatype cells. Because one record (e.g. a chloroplast genome) matches many filters at once, this **double-counts** and is *not* a record count — it is kept only for reference. Compare it to the deduplicated union above.

**Cells Failed (NA):** number of datatype cells for this species that failed and were recorded as NA.

**Name-based Total (NCBI Synonyms):** sum of datatype counts using name-based NCBI synonyms query.

**Name-based Total (NCBI+AI):** sum of datatype counts using name-based NCBI+AI query.

**Name-based Total Gain (AI - NCBI):** Name-based Total (NCBI+AI) minus Name-based Total (NCBI Synonyms).

**Total NCBI (no filters):** unfiltered count using primary query.

**Name-based NCBI (no filters):** unfiltered count using name-based NCBI synonyms query.

**Name-based NCBI+AI (no filters):** unfiltered count using name-based NCBI+AI query.

**Name-based Gain (no filters):** Name-based NCBI+AI minus Name-based NCBI.

## Outputs
Each run writes a timestamped set of files to `missmap_results/`:

- `data_availability_matrix_<ts>.csv` / `.xlsx` — the full matrix.
- `missmap_presence_<ts>.png` — species × datatype **presence / absence** map (the "map" in MissMap). Blue = data present, grey = absent, red = NA. Genomic categories and barcode markers are grouped with a divider so they are never read as an additive whole.
- `missmap_counts_<ts>.png` — the same grid coloured by `log10(count+1)` so genome-scale categories don't visually swamp single-amplicon barcodes.
- `coverage_summary_<ts>.csv` — per-datatype coverage (% of species with data) plus a printed gap list of which species are missing each standard barcode.
- `run_manifest_<ts>.json` — run date, Entrez db/context, options, AI model, and any failed cells, so a result is reproducible and citable (NCBI counts change daily).

## Reliability
- Every count flows through one helper that **retries transient HTTP 429/5xx and network errors with exponential backoff** and records a failure as `NA`, never `0`.
- Counts are grouped as *genomic resources* vs *barcode markers*; the deduplicated union count replaces the old summed "total" as the headline figure.

## Flags
`--no-ai`, `--no-subspecies`, `--ai-unverified`, `--ai-validated-only`, `--verbose`, `--debug-queries` (print exact query strings), and `--non-interactive` (never prompt — for HPC/batch/CI; also auto-enabled when stdin is not a TTY).

## Testing
`pip install -r requirements.txt` then `pytest -q`. Tests mock NCBI Entrez and cover NA-vs-0, retry/backoff, the union-vs-sum distinction, the heatmaps, and the infraspecific filter.

## Security
Do **not** commit real credentials. `config.txt` is gitignored, but prefer environment variables or a gitignored `.env` (see `.env.example`) so secrets are never written to the working tree. If a key has ever been committed or shared, rotate it (NCBI: Account Settings → API Key Management; Groq: console API keys).

