"""Download GBIF occurrence records for a species or genus.

Usage: python fetch_gbif.py "<taxon name>" <output_name>
Example: python fetch_gbif.py "Tetranychus urticae" gbif_tetranychus_urticae
"""
import csv
import json
import os
import sys
import time
import urllib.parse
import urllib.request

API = "https://api.gbif.org/v1"
FIELDS = ["gbifID", "scientificName", "acceptedScientificName", "decimalLatitude",
          "decimalLongitude", "coordinateUncertaintyInMeters", "country", "countryCode",
          "stateProvince", "locality", "year", "month", "eventDate", "basisOfRecord",
          "institutionCode", "catalogNumber", "datasetName", "recordedBy",
          "identifiedBy", "associatedTaxa", "habitat", "issues"]


def get_json(url):
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                return json.loads(r.read().decode("utf-8"))
        except Exception as err:
            print(f"Retry {attempt + 1}: {err}")
            time.sleep(10)
    sys.exit(f"Failed: {url}")


def main():
    if len(sys.argv) != 3:
        sys.exit('Usage: python fetch_gbif.py "<taxon>" <output_name>')
    taxon, name = sys.argv[1], sys.argv[2]

    match = get_json(f"{API}/species/match?name={urllib.parse.quote(taxon)}")
    key = match.get("usageKey")
    if not key:
        sys.exit(f"GBIF found no taxon for {taxon}: {match}")
    print(f"GBIF taxon: {match.get('scientificName')} (key {key}, {match.get('matchType')})")

    os.makedirs("data", exist_ok=True)
    path = os.path.join("data", f"{name}.csv")
    offset, limit, total = 0, 300, 0
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        w.writeheader()
        while True:
            page = get_json(f"{API}/occurrence/search?taxonKey={key}"
                            f"&occurrenceStatus=PRESENT&limit={limit}&offset={offset}")
            for rec in page.get("results", []):
                row = {k: rec.get(k, "") for k in FIELDS}
                row["issues"] = ";".join(rec.get("issues", []))
                w.writerow(row)
                total += 1
            print(f"{total}/{page.get('count')}")
            if page.get("endOfRecords") or offset + limit >= 100000:
                break
            offset += limit
            time.sleep(0.3)
    print(f"Saved {total} occurrences to {path}")

    with open(os.path.join("data", f"{name}_query.txt"), "w") as f:
        f.write(f"GBIF taxon: {match.get('scientificName')} (taxonKey {key})\n"
                f"Records: {total}\n"
                f"Downloaded: {time.strftime('%Y-%m-%d %H:%M UTC', time.gmtime())}\n"
                "Cite GBIF as required at https://www.gbif.org/citation-guidelines\n")


if __name__ == "__main__":
    main()
