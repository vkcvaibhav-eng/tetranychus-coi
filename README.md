# tetranychus-coi

Reproducible pipeline for collecting and auditing public COI (cytochrome c oxidase subunit I) barcode data for *Tetranychus urticae* and other *Tetranychus* species.

The pipeline gathers COI records from GenBank and BOLD, occurrence records from GBIF, and climate values from WorldClim. It then merges the records, removes duplicates, scores how well each record is documented, and checks every sequence against a reference COX1 gene for problems such as stop codons and frameshifts.

## Pipeline

| Step | Script | What it does |
|---|---|---|
| 1 | `fetch_ncbi.py`, `fetch_bold.py`, `fetch_gbif.py` | Download sequences and occurrence records |
| 1 | `climate_extract.py` | Attach WorldClim BIO1–BIO19 to every record with coordinates |
| 2 | `fetch_reference.py` | Download the reference mitogenome and extract its COX1 gene |
| 2 | `merge_records.py` | Merge GenBank and BOLD, remove duplicates, score metadata (0–7) |
| 3 | `sequence_qc.py` | Orientation, COX1 region covered, identity, ambiguous bases, internal stop codons (genetic code 5), frameshifts |
| 4 | `species_analysis.py` | K2P distances on shared sites, barcode gap per species, shared haplotypes, nearest-neighbour label conflicts, distance clusters |
| 5 | `anchor_check.py` | Compares clusters and records with anchor sequences from morphologically identified, published specimens |

## Running it on GitHub

1. Add repository secrets under **Settings → Secrets and variables → Actions**: `NCBI_EMAIL` (required) and `NCBI_API_KEY` (optional).
2. **Actions → Fetch sequences and occurrence data → Run workflow** downloads the data into `data/`.
3. **Actions → Run analysis → Run workflow** runs steps 2–4 and saves the output in `results/`.

BOLD may refuse requests from GitHub's servers. In that case download the records from https://portal.boldsystems.org as TSV (BCDM) and add them to `data/` as `bold_<name>.tsv`.

## Running it locally

```
pip install -r requirements.txt
export NCBI_EMAIL=you@example.org
python scripts/fetch_ncbi.py 'Tetranychus urticae[Organism] AND (COI[Gene] OR COX1[Gene])' ncbi_tetranychus_urticae
python scripts/fetch_bold.py "Tetranychus urticae" bold_tetranychus_urticae
python scripts/fetch_reference.py NC_010526 "Tetranychus urticae"
python scripts/merge_records.py
python scripts/sequence_qc.py
python scripts/species_analysis.py --results results --out results/species
pytest -q
```

## Outputs

`results/urticae/` (focal species audit) and `results/genus/` (all *Tetranychus*) each contain:

- `merged_records.csv`: one row per unique record with metadata and a completeness score
- `merged_sequences.fasta`: all unique sequences
- `qc_records.csv`: quality-check result and flags for each sequence
- `qc_passed.fasta`: sequences that passed, in forward orientation
- `merge_summary.txt`, `qc_summary.txt`: counts for reporting

`results/species/` contains `record_checks.csv`, `haplotypes.csv`, `barcode_gap.csv`, `clusters.csv` and `species_summary.txt`. Step logs are in `results/logs/` and `data/logs/`.

GenBank records longer than 2,500 bp (mitochondrial genomes) contribute only their annotated COX1 gene. RefSeq records that copy an original submission are removed.

## Metadata rules

- **Voucher:** for records BOLD copied from GenBank ("Mined from GenBank, NCBI"), the BOLD sample ID is the GenBank accession and is not counted as a voucher; only a museum ID is.
- **Identifier:** identifications made by the BOLD Identification System are flagged in `identification_automated` and not counted, because they come from matching the sequence against the database itself.
- **Coordinates:** BOLD coordinates placed at a country or province centroid are excluded from the metadata score and from the climate extraction; the reason is kept in `coordinate_source`.

## QC flags

| Flag | Meaning |
|---|---|
| `SHORT` | Shorter than 300 bp |
| `AMBIGUOUS` | More than 1% non-ACGT bases |
| `NOT_ALIGNED` | Does not align to the reference COX1 |
| `LOW_IDENTITY` | Identity to the reference below the set threshold |
| `PARTIAL_ALIGNMENT` | Less than 90% of the sequence aligns |
| `STOP_CODON` | Internal stop codon in the reference reading frame |
| `FRAMESHIFT` | Indel that is not a multiple of three |

Stop codons and frameshifts in a protein-coding mitochondrial gene usually point to sequencing errors or nuclear mitochondrial pseudogenes (NUMTs). A frameshift is flagged only when reading the aligned stretch in one frame produces stop codons; gaps of odd length with an intact frame are alignment artefacts against a divergent reference and are reported in `alignment_frameshift_gaps` without a flag.

## Species analysis notes

- Sequences are placed on the reference COX1 coordinates from pairwise alignments, and insertions relative to the reference are dropped. This is adequate for screening; final phylogenetic analyses should use a codon-aware multiple alignment.
- Distances are computed only over sites two sequences share, with a minimum overlap of 300 bp, so 5′ and 3′ fragments are not compared.
- A mislabelled record also changes the result for correctly labelled neighbours: they can be flagged `NEAREST_IS_OTHER_SPECIES`, and the species carrying the wrong label loses its barcode gap. Flags mark places to check, not proven errors.

## Anchor records

`data/anchors/anchor_studies.csv` lists published studies whose *Tetranychus* specimens were identified morphologically, with separate fields for morphology, males examined, named identifier, voucher and crossing tests. `anchor_records.csv` gives their COI accessions, each confirmed to exist in the downloaded data. Evidence levels:

- **A**: morphology with males examined and vouchers in a named collection
- **B**: morphology stated, vouchers not confirmed
- **C**: lead only, not used for decisions

Within a cluster, level A anchors outrank level B. Names are compared after the synonymies in `synonyms.csv`, each with its published source. Results are in `results/anchors/`.

## Data sources

GenBank (NCBI), BOLD Systems, GBIF and WorldClim 2.1. Cite each source as its terms require; GBIF downloads need a citation following https://www.gbif.org/citation-guidelines.

## Author

V. K. Chaudhari, AINP on Agricultural Acarology, Department of Entomology, N. M. College of Agriculture, Navsari Agricultural University, Navsari, Gujarat, India.
