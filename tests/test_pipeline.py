"""Tests for the merge and sequence-QC modules using synthetic records.

Run: pytest -q
"""
import csv
import random
import subprocess
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, Reference, SeqFeature
from Bio.SeqRecord import SeqRecord

ROOT = Path(__file__).resolve().parent.parent
STOPS = {"TAA", "TAG"}  # stop codons in the invertebrate mitochondrial code
random.seed(1)


def random_orf(n_codons):
    codons = []
    while len(codons) < n_codons:
        c = "".join(random.choice("ACGT") for _ in range(3))
        if c not in STOPS:
            codons.append(c)
    return "ATG" + "".join(codons) + "TAA"


def mutate_third_positions(seq, rate):
    s = list(seq)
    for i in range(2, len(s), 3):
        if random.random() < rate:
            old = s[i]
            s[i] = random.choice([b for b in "ACGT" if b != old])
            if "".join(s[i - 2:i + 1]) in STOPS:
                s[i] = old
    return "".join(s)


REF = random_orf(510)  # ~1.5 kb, similar to COX1


def gb_record(acc, seq, organism="Tetranychus urticae", country="India: Gujarat",
              lat_lon=None, voucher=None, published=False):
    rec = SeqRecord(Seq(seq), id=acc + ".1", name=acc, description="test COI")
    rec.annotations["molecule_type"] = "DNA"
    rec.annotations["organism"] = organism
    q = {"organism": [organism], "geo_loc_name": [country], "host": ["okra"],
         "collection_date": ["2024"]}
    if lat_lon:
        q["lat_lon"] = [lat_lon]
    if voucher:
        q["specimen_voucher"] = [voucher]
    rec.features = [SeqFeature(FeatureLocation(0, len(seq)), type="source", qualifiers=q),
                    SeqFeature(FeatureLocation(0, len(seq)), type="gene", qualifiers={"gene": ["COI"]})]
    ref = Reference()
    ref.title = "Population genetics of Tetranychus" if published else "Direct Submission"
    ref.journal = "Exp Appl Acarol 1:1-10" if published else "Submitted (01-JAN-2024)"
    rec.annotations["references"] = [ref]
    return rec


def run(script, *args, cwd):
    res = subprocess.run([sys.executable, str(ROOT / "scripts" / script), *args],
                         cwd=cwd, capture_output=True, text=True)
    assert res.returncode == 0, res.stderr + res.stdout
    return res.stdout


def read_csv(path):
    with open(path, newline="") as f:
        return {r["record_id"]: r for r in csv.DictReader(f)}


def build_dataset(tmp):
    data = tmp / "data"
    (data / "reference").mkdir(parents=True)
    (data / "reference" / "coi_reference.fasta").write_text(f">ref_COX1\n{REF}\n")

    good = mutate_third_positions(REF[90:750], 0.08)
    stop = good[:300] + "TAA" + good[303:]
    shift = good[:400] + good[401:]
    three_prime = mutate_third_positions(REF[800:1450], 0.08)
    ambiguous = good[:200] + "N" * 20 + good[220:]
    foreign = "".join(random.choice("ACGT") for _ in range(650))
    reverse = str(Seq(good).reverse_complement())

    file1 = [gb_record("AB000001", good, lat_lon="20.95 N 72.92 E", published=True),
             gb_record("AB000002", stop),
             gb_record("AB000003", shift),
             gb_record("AB000004", three_prime),
             gb_record("AB000005", good[:200]),
             gb_record("AB000006", ambiguous),
             gb_record("AB000007", foreign),
             gb_record("AB000008", reverse)]
    file2 = [gb_record("AB000001", good, lat_lon="20.95 N 72.92 E", published=True)]
    SeqIO.write(file1, data / "ncbi_t1.gb", "genbank")
    SeqIO.write(file2, data / "ncbi_t2.gb", "genbank")

    header = ["processid", "sampleid", "museumid", "voucher_type", "species",
              "country/ocean", "coord", "coord_source", "nuc", "marker_code",
              "insdc_acs", "identified_by", "bin_uri", "inst"]
    mined = "Mined from GenBank, NCBI"
    rows = [
        ["MITE001-20", "V-1", "V-1", "Vouchered:Registered Collection", "Tetranychus urticae",
         "India", "[21.1, 72.8]", "GPS", good, "COI-5P", "AB000002.1", "A. Taxonomist",
         "BOLD:AAA0001", "NAU"],
        ["MITE002-20", "V-2", "", "Vouchered:Registered Collection", "Tetranychus urticae",
         "Spain", "[40.4, -3.7]", "GPS", good, "COI-5P", "", "BOLD Identification System [2023]",
         "BOLD:AAA0001", "Museum"],
        ["MITE003-20", "V-3", "", "", "Tetranychus urticae", "Spain", "", "", "ACGTACGT",
         "ITS2", "", "", "", ""],
        # Copied from GenBank: sampleid is the accession, coordinate is a country centroid.
        ["GBMIN001-13", "AB000003", "", "", "Tetranychus urticae", "Russia",
         "[59.039433, 98.6705]", "Coordinates from country centroid", shift, "COI-5P",
         "AB000003.1", "", "", mined],
        ["GBMIN002-13", "XY999999", "", "", "Tetranychus urticae", "China",
         "[38.073254, 104.69114]", "Coordinates from country centroid", good, "COI-5P",
         "XY999999.1", "", "", mined],
    ]
    with open(data / "bold_t.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(header)
        w.writerows(rows)
    return tmp


def test_merge_and_qc(tmp_path):
    build_dataset(tmp_path)
    run("merge_records.py", cwd=tmp_path)
    merged = read_csv(tmp_path / "results" / "merged_records.csv")

    # Duplicate GenBank accession across two files is kept once.
    assert len([k for k in merged if k.startswith("AB000001")]) == 1
    # BOLD copy of a GenBank record is merged, not counted twice.
    assert merged["AB000002.1"]["source"] == "GenBank+BOLD"
    assert merged["AB000002.1"]["voucher"] == "V-1"
    assert merged["AB000002.1"]["identified_by"] == "A. Taxonomist"
    assert "MITE001-20" not in merged
    # BOLD-only COI record kept; non-COI marker skipped.
    assert merged["MITE002-20"]["source"] == "BOLD"
    assert "MITE003-20" not in merged
    assert len(merged) == 10
    # GenBank accession used as BOLD sampleid is not a voucher.
    assert merged["AB000003.1"]["voucher"] == ""
    assert merged["GBMIN002-13"]["voucher"] == ""
    # Country-centroid coordinates are not used as locations.
    assert merged["GBMIN002-13"]["latitude"] == ""
    assert "excluded" in merged["GBMIN002-13"]["coordinate_source"]
    assert merged["AB000003.1"]["latitude"] == ""
    # A directly submitted, vouchered BOLD record keeps its sample ID as voucher.
    assert merged["MITE002-20"]["voucher"] == "V-2"
    # Automated BOLD identification does not count as an identifier.
    assert merged["MITE002-20"]["identification_automated"] == "yes"
    assert "identifier" in merged["MITE002-20"]["metadata_missing"]
    assert "identifier" not in merged["AB000002.1"]["metadata_missing"]
    # Metadata score: the published record with coordinates scores higher.
    assert int(merged["AB000001.1"]["metadata_score"]) > int(merged["AB000003.1"]["metadata_score"])
    assert float(merged["AB000001.1"]["latitude"]) == 20.95

    run("sequence_qc.py", cwd=tmp_path)
    qc = read_csv(tmp_path / "results" / "qc_records.csv")
    assert qc["AB000001.1"]["qc_status"] == "PASS"
    assert qc["AB000001.1"]["cox1_half"] == "5prime"
    assert "STOP_CODON" in qc["AB000002.1"]["qc_flags"]
    assert "FRAMESHIFT" in qc["AB000003.1"]["qc_flags"]
    assert qc["AB000004.1"]["qc_status"] == "PASS"
    assert qc["AB000004.1"]["cox1_half"] == "3prime"
    assert "SHORT" in qc["AB000005.1"]["qc_flags"]
    assert "AMBIGUOUS" in qc["AB000006.1"]["qc_flags"]
    assert qc["AB000007.1"]["qc_status"] == "FAIL"
    assert qc["AB000008.1"]["orientation"] == "reverse"
    assert qc["AB000008.1"]["qc_status"] == "PASS"

    passed = {r.id: str(r.seq) for r in SeqIO.parse(tmp_path / "results" / "qc_passed.fasta", "fasta")}
    # Reverse-complemented input is stored in the forward orientation.
    assert passed["AB000008.1"] == passed["AB000001.1"]
