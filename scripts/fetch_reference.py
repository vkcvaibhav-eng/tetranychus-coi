"""Download a reference mitochondrial genome from NCBI and save its COI (COX1)
coding sequence as the alignment reference for the quality check.

Usage: python scripts/fetch_reference.py <accession> [expected organism]
Example: python scripts/fetch_reference.py NC_010526 "Tetranychus urticae"
Needs NCBI_EMAIL in the environment; NCBI_API_KEY is optional.
"""
import os
import sys

from Bio import Entrez, SeqIO

OUT_DIR = os.path.join("data", "reference")


def is_cox1(feature):
    q = feature.qualifiers
    names = [n.lower() for n in q.get("gene", []) + q.get("product", [])]
    return any(n in ("cox1", "coi", "co1", "coxi", "cytochrome c oxidase subunit i",
                     "cytochrome oxidase subunit i", "cytochrome c oxidase subunit 1")
               for n in names)


def main():
    if len(sys.argv) < 2:
        sys.exit("Usage: python scripts/fetch_reference.py <accession> [organism]")
    acc = sys.argv[1]
    expected = sys.argv[2] if len(sys.argv) > 2 else None

    email = (os.environ.get("NCBI_EMAIL") or "").strip()
    if not email:
        sys.exit("ERROR: the NCBI_EMAIL secret is empty or missing.")
    Entrez.email = email
    api_key = (os.environ.get("NCBI_API_KEY") or "").strip()

    os.makedirs(OUT_DIR, exist_ok=True)
    gb_path = os.path.join(OUT_DIR, f"{acc}.gb")
    text = None
    for key_value in ([api_key, None] if api_key else [None]):
        Entrez.api_key = key_value
        try:
            with Entrez.efetch(db="nuccore", id=acc, rettype="gb", retmode="text") as h:
                text = h.read()
            break
        except Exception as err:
            print(f"NCBI fetch failed ({'with' if key_value else 'without'} API key): {err}")
    if not text:
        sys.exit("ERROR: could not download the reference record.")
    with open(gb_path, "w") as out:
        out.write(text)

    rec = SeqIO.read(gb_path, "genbank")
    organism = rec.annotations.get("organism", "")
    print(f"{rec.id}: {organism}, {len(rec.seq)} bp")
    if expected and expected.lower() not in organism.lower():
        sys.exit(f"Reference organism is '{organism}', expected '{expected}'")

    cds = [f for f in rec.features if f.type == "CDS" and is_cox1(f)]
    if not cds:
        sys.exit("No COX1 CDS annotated in the reference record")
    seq = cds[0].extract(rec.seq)
    path = os.path.join(OUT_DIR, "coi_reference.fasta")
    with open(path, "w") as f:
        f.write(f">{rec.id}_COX1 {organism}\n{seq}\n")
    print(f"COX1 CDS: {len(seq)} bp saved to {path}")


if __name__ == "__main__":
    main()
