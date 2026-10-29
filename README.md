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

## Running it on GitHub

1. Add repository secrets under **Settings → Secrets and variables → Actions**: `NCBI_EMAIL` (required) and `NCBI_API_KEY` (optional).
2. **Actions → Fetch sequences and occurrence data → Run workflow** downloads the data into `data/`.
3. **Actions → Run analysis → Run workflow** runs steps 2–3 and saves the output in `results/`.

## Running it locally

```
pip install -r requirements.txt
export NCBI_EMAIL=you@example.org
python scripts/fetch_ncbi.py 'Tetranychus urticae[Organism] AND (COI[Gene] OR COX1[Gene])' ncbi_tetranychus_urticae
python scripts/fetch_bold.py "Tetranychus urticae" bold_tetranychus_urticae
python scripts/fetch_reference.py NC_010526 "Tetranychus urticae"
python scripts/merge_records.py
python scripts/sequence_qc.py
pytest -q
```

## Outputs (`results/`)

- `merged_records.csv`: one row per unique record with metadata and a completeness score
- `merged_sequences.fasta`: all unique sequences
- `qc_records.csv`: quality-check result and flags for each sequence
- `qc_passed.fasta`: sequences that passed, in forward orientation
- `merge_summary.txt`, `qc_summary.txt`: counts for reporting

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

Stop codons and frameshifts in a protein-coding mitochondrial gene usually point to sequencing errors or nuclear mitochondrial pseudogenes (NUMTs).

## Data sources

GenBank (NCBI), BOLD Systems, GBIF and WorldClim 2.1. Cite each source as its terms require; GBIF downloads need a citation following https://www.gbif.org/citation-guidelines.

## Author

V. K. Chaudhari, AINP on Agricultural Acarology, Department of Entomology, N. M. College of Agriculture, Navsari Agricultural University, Navsari, Gujarat, India.
