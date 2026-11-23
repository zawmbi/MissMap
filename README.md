# MissMap
How much NCBI sequence data does a name hide? MissMap counts GenBank nucleotide records per
species, reconciles each species' synonymy across taxonomic authorities (WCVP, WFO, GBIF, any
Darwin Core checklist), and measures the records a TaxID search **misses** (filed under names
NCBI treats as separate taxa) and the records a free-text name search **wrongly includes**.
An optional LLM step proposes extra synonyms, and a separate benchmark compares LLMs at it.

## Quick start
```bash
pip install -r requirements.txt
python MissMap.py -s "Solanum lycopersicum,Lampranthus emarginatoides" -dt "rbcL,matK,ITS" \
    --wcvp wcvp_names.csv --wfo classification.csv --gbif --group Viridiplantae
```

## Config
MissMap reads credentials and settings from `config.txt`, environment variables, or a `.env` file
(see `.env.example`). Keys: `ENTREZ_EMAIL`, `ENTREZ_API_KEY`, `LLM_MODEL` (a `provider:model` spec,
below) plus the matching API key (`ANTHROPIC_API_KEY`, `GROQ_API_KEY`, `OPENROUTER_API_KEY`,
`OPENAI_API_KEY`). The legacy `GROQ_MODEL` is still read if `LLM_MODEL` is unset. Optional
`config.txt` keys: `species`, `datatypes`, `wcvp`, `wfo`, `group`, `taxdump`.

LLM specs: `anthropic:claude-sonnet-5-5`, `anthropic:claude-opus-5-5+web` (with web search),
`groq:llama-3.3-70b-versatile`, `openrouter:<vendor/model>`, `ollama:<model>` (local).

## Authorities and IDs
Each source has its own identifiers; they are **not** interchangeable:

| Source | ID | Bridge used by MissMap |
|---|---|---|
| NCBI Taxonomy | `taxid` | exact canonical name (NCBI stores no IPNI/WFO IDs) |
| WCVP / POWO | `plant_name_id`, `powo_id`, `ipni_id` | IPNI ID links to WFO |
| WFO | `taxonID` (`wfo-...`), IPNI in `scientificNameID` | IPNI ID links to WCVP |
| GBIF backbone | `usageKey` | name match API |
| Wikidata | item `Q...` holding NCBI/GBIF/IPNI/POWO/WFO IDs | independent check (`verify_crosswalk`) |

Every crosswalk row records how it was matched (`matchMethod`). TaxIDs are resolved by name **with
rank and lineage checks** (`--group`), so cross-kingdom homonyms (*Prunella*, *Morus*, *Oenanthe*)
are flagged instead of silently taking the first hit. Download the WCVP release (Kew) and the WFO
backbone (`classification.csv`) and keep the files: their fingerprints go into the run manifest.

## Name buckets
Every AI or authority synonym is classified (`name_classification_<ts>.csv`, `AI: <bucket>` columns):

| Bucket | Meaning |
|---|---|
| `same_taxid` | NCBI already files the name under this TaxID (nothing new to find) |
| `other_taxid_confirmed` | NCBI files it under **another** TaxID and an authority calls it a synonym: **hidden records** |
| `other_taxid_wrong_species` | NCBI files it elsewhere and an authority calls it a different species: adding it would contaminate |
| `other_taxid_unconfirmed` | NCBI files it elsewhere; no loaded authority knows it |
| `valid_not_in_ncbi` | a real name (authority/GBIF) with no NCBI taxon |
| `not_found` | not found anywhere checked (possible hallucination) |
| `not_a_scientific_name` | common name, cultivar, malformed |

## AI Disclaimer
AI output can vary by model and input context. Synonym lists may differ across runs,
models, or prompt updates. Treat AI synonyms as suggestions to be verified.

## Columns
**Species:** the inputted species name(s).

**NCBI Synonyms:** synonyms from NCBI Taxonomy Database.

**AI Proposed Synonyms:** AI output list (unverified).

**AI-validated Synonyms:** AI synonyms in bucket `same_taxid` or `other_taxid_confirmed` (see Name buckets).

**AI Rejected Synonyms:** all other AI synonyms.

**AI Validation Source:** sources used for validation (NCBI Taxonomy + loaded authorities + GBIF).

**AI: \<bucket\>:** AI names in each name bucket.

**TaxID / TaxID Flags:** the resolved NCBI TaxID and any warnings (`homonym_other_group`, `matched_via_synonym`, `rank_mismatch`, `ambiguous`; ambiguous names get no TaxID rather than a guess).

**NCBI Common Names:** vernacular names from NCBI, kept out of the synonym search.

**\<Authority\> Accepted Name / \<Authority\> Synonyms:** what each loaded checklist accepts and lists as synonyms.

**Accepted Name Agrees Across Authorities:** whether all checklists resolve the species to the same accepted name.

**Authority Synonyms Missing From NCBI:** synonyms any authority lists that NCBI's synonym list lacks.

**Hidden TaxIDs (authority synonyms / AI synonyms):** NCBI TaxIDs of confirmed synonyms that are not this TaxID or below it.

**Hidden Records (authority synonyms / AI synonyms / all):** records under those TaxIDs that a search on this TaxID does **not** return (`(txidA OR txidB)[Organism:exp] NOT txidP[Organism:exp]`, so overlaps are never double counted).

**Hidden Share of Reconciled Total (%):** hidden / (TaxID total + hidden).

**Markers Gained via Synonyms:** datatypes that are 0 under the TaxID but present under hidden TaxIDs, e.g. `rbcL(+2)`.

**Free-text Name Search (All Fields):** what typing the name into NCBI returns.

**Free-text Hits Not Under TaxID (over-count):** free-text hits that are not records of this taxon (they only mention it).

**Free-text Precision (%):** share of free-text hits that are records of this taxon.

**TaxID Records Excl. Infraspecific (noexp):** records on the species TaxID itself, without subspecies/varieties.

**NCBI Synonyms / Child Taxa Authorities Treat As Separate Species, Records Under Those Child Taxa:** what NCBI lumps into this TaxID that an authority splits (records wrongly included).

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
- `missmap_presence_reconciled_<ts>.png` — the presence map after adding hidden records (only when something changed).
- `name_classification_<ts>.csv` — every AI/authority name with its bucket, NCBI TaxIDs and authority verdicts.
- `synonymy_comparison_<ts>.json` — per species: each authority's accepted name and synonyms, pairwise Jaccard, and a class per name (`agreed`, `lump_split`, `name_missing`, `status_differs`).
- `crosswalk_dwc_<ts>.csv` — Darwin Core rows (`taxonID`, `acceptedNameUsageID`, `taxonomicStatus`, `scientificNameAuthorship`, external IDs) per authority with the NCBI TaxIDs each name maps to.
- `llm_calls.jsonl` — every LLM call with raw response, exact model ID, prompt version, latency and tokens. Also a cache: a finished call is never re-billed.
- `entrez_cache.sqlite` — same-day Entrez answers (`--cache-hours 0` freezes them for a thesis dataset).
- `run_manifest_<ts>.json` — run date, Entrez db/context, options, AI model + prompt version, authority file fingerprints, TaxID resolution flags, cache stats, and any failed cells, so a result is reproducible and citable (NCBI counts change daily).

## Reliability
- Every count flows through one helper that **retries transient HTTP 429/5xx and network errors with exponential backoff** and records a failure as `NA`, never `0`.
- Counts are grouped as *genomic resources* vs *barcode markers*; the deduplicated union count replaces the old summed "total" as the headline figure.

## Flags
`--no-ai`, `--no-subspecies`, `--ai-unverified`, `--ai-validated-only`, `--verbose`, `--debug-queries` (print exact query strings), and `--non-interactive` (never prompt — for HPC/batch/CI; also auto-enabled when stdin is not a TTY).

Reconciliation: `--wcvp PATH`, `--wfo PATH`, `--dwc NAME=PATH` (repeatable; e.g. MDD, ASW, COL exports), `--gbif`, `--group Viridiplantae|Mammalia|...`, `--taxdump DIR` (offline NCBI `names.dmp`/`nodes.dmp`), `--llm provider:model`, `--no-cache`, `--cache-hours`.

## Study tools
- `missmap_sample.py` — seeded, stratified random sample of accepted species (synonym-count strata, per-family cap, design weights).
- `missmap_benchmark.py run` — query several LLMs `k` times per species (logged, cached).
- `missmap_benchmark.py score` — precision / recall / F1 / hallucination rate with Wilson CIs, run-to-run Jaccard, Cochran's Q and pairwise McNemar (Holm-adjusted), against each authority, a consensus set, and an optional hand-checked gold set.
- `missmap_stats.py` — the statistics above plus Fleiss' kappa and bootstrap CIs.

## Testing
`pip install -r requirements.txt` then `pytest -q`. Tests mock NCBI Entrez, GBIF and the LLM APIs and cover NA-vs-0, retry/backoff, the union-vs-sum distinction, the heatmaps, the infraspecific filter, name parsing/matching, authority loaders, homonym-safe TaxID resolution, the hidden-record scenario (*Lampranthus*), the LLM harness, the statistics, the benchmark scorer and the sampler.

If the heatmap tests fail with `numpy.core.multiarray failed to import`, your environment has numpy 2.x with a matplotlib built for numpy 1.x; install the pinned versions in a fresh virtual environment.

## Security
Do **not** commit real credentials. `config.txt` is gitignored, but prefer environment variables or a gitignored `.env` (see `.env.example`) so secrets are never written to the working tree. If a key has ever been committed or shared, rotate it (NCBI: Account Settings → API Key Management; Groq: console API keys).

