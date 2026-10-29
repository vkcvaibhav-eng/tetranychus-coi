"""Download public BOLD Systems records for a taxon (BOLD Portal API) and
write the raw table plus a FASTA of the COI-5P barcode sequences.

Usage: python fetch_bold.py "<taxon name>" <output_name>
Example: python fetch_bold.py "Tetranychus urticae" bold_tetranychus_urticae
"""
import csv
import io
import json
import os
import sys
import time
import urllib.parse
import urllib.request

BASE = "https://portal.boldsystems.org/api"


def get(url, retries=3, fatal=True):
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "tetranychus-coi-study"})
            with urllib.request.urlopen(req, timeout=300) as r:
                return r.read().decode("utf-8")
        except Exception as err:
            print(f"Retry {attempt + 1}: {err}")
            time.sleep(10)
    if fatal:
        sys.exit(f"Failed: {url}\nBOLD refused the request. Download the records from "
                 "https://portal.boldsystems.org yourself (TSV format) and upload the file "
                 "to data/ as bold_<name>.tsv; the merge step reads it the same way.")
    return None


def find_query_id(obj):
    if isinstance(obj, dict):
        for k in ("query_id", "queryId", "id"):
            if isinstance(obj.get(k), str):
                return obj[k]
        for v in obj.values():
            found = find_query_id(v)
            if found:
                return found
    return None


def main():
    if len(sys.argv) != 3:
        sys.exit('Usage: python fetch_bold.py "<taxon>" <output_name>')
    taxon, name = sys.argv[1], sys.argv[2]
    rank = "species" if " " in taxon.strip() else "genus"
    query = f"tax:{rank}:{taxon}"
    os.makedirs("data", exist_ok=True)

    check = get(f"{BASE}/query/preprocessor?query={urllib.parse.quote(query)}", retries=1, fatal=False)
    print("Preprocessor:", (check or "skipped")[:500])

    raw = get(f"{BASE}/query?query={urllib.parse.quote(query)}&extent=full")
    qid = find_query_id(json.loads(raw))
    if not qid:
        sys.exit(f"No query id in response: {raw[:500]}")
    print("Query id:", qid)

    text = get(f"{BASE}/documents/{urllib.parse.quote(qid, safe='')}/download?format=tsv")
    rows = []
    if "\t" in text.split("\n", 1)[0]:
        table_path = os.path.join("data", f"{name}.tsv")
        with open(table_path, "w") as f:
            f.write(text)
        rows = list(csv.DictReader(io.StringIO(text), delimiter="\t"))
    else:
        text = get(f"{BASE}/documents/{urllib.parse.quote(qid, safe='')}/download?format=json")
        table_path = os.path.join("data", f"{name}.jsonl")
        with open(table_path, "w") as f:
            f.write(text)
        stripped = text.strip()
        if stripped.startswith("["):
            rows = json.loads(stripped)
        else:
            rows = [json.loads(line) for line in stripped.splitlines() if line.strip()]
    print(f"{len(rows)} records saved to {table_path}")

    n = 0
    with open(os.path.join("data", f"{name}_COI.fasta"), "w") as fa:
        for r in rows:
            seq = (r.get("nuc") or "").replace("-", "").strip()
            marker = r.get("marker_code") or ""
            if seq and marker.startswith("COI"):
                label = "|".join(str(r.get(k) or "NA") for k in
                                 ("processid", "species", "country/ocean", "insdc_acs"))
                fa.write(f">{label.replace(' ', '_')}\n{seq}\n")
                n += 1
    print(f"{n} COI sequences written")

    with open(os.path.join("data", f"{name}_query.txt"), "w") as f:
        f.write(f"BOLD query: {query}\nRecords: {len(rows)}\nCOI sequences: {n}\n"
                f"Downloaded: {time.strftime('%Y-%m-%d %H:%M UTC', time.gmtime())}\n")


if __name__ == "__main__":
    main()
