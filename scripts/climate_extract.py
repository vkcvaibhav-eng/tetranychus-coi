"""Attach WorldClim 2.1 bioclimatic values (BIO1-BIO19) to every record with
coordinates in the data folder (GenBank .gb, BOLD .tsv/.jsonl, GBIF .csv).

Usage: python climate_extract.py [resolution]   (10m, 5m or 2.5m; default 10m)
The WorldClim rasters are downloaded to a temporary folder and not stored
in the repository; only the extracted table is saved.
"""
import csv
import glob
import json
import os
import re
import sys
import tempfile
import urllib.request
import zipfile

import rasterio
from Bio import SeqIO

RES = sys.argv[1] if len(sys.argv) > 1 else "10m"
URL = f"https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_{RES}_bio.zip"
NUM = re.compile(r"-?\d+(?:\.\d+)?")


def genbank_points():
    for path in glob.glob("data/*.gb"):
        for rec in SeqIO.parse(path, "genbank"):
            src = next((f for f in rec.features if f.type == "source"), None)
            if not src or "lat_lon" not in src.qualifiers:
                continue
            m = re.match(r"([\d.]+)\s*([NS])\s+([\d.]+)\s*([EW])", src.qualifiers["lat_lon"][0])
            if not m:
                continue
            lat = float(m.group(1)) * (1 if m.group(2) == "N" else -1)
            lon = float(m.group(3)) * (1 if m.group(4) == "E" else -1)
            q = src.qualifiers
            yield {"source": "GenBank", "record_id": rec.id,
                   "taxon": q.get("organism", [""])[0], "country": q.get("country", [""])[0],
                   "host": q.get("host", [""])[0], "lat": lat, "lon": lon}


def bold_points():
    rows = []
    for path in glob.glob("data/bold_*.tsv"):
        with open(path) as f:
            rows += list(csv.DictReader(f, delimiter="\t"))
    for path in glob.glob("data/bold_*.jsonl"):
        with open(path) as f:
            text = f.read().strip()
        rows += json.loads(text) if text.startswith("[") else \
            [json.loads(l) for l in text.splitlines() if l.strip()]
    for r in rows:
        if "centroid" in (r.get("coord_source") or "").lower():
            continue  # country/province centre, not a collection site
        coord = r.get("coord")
        if isinstance(coord, list) and len(coord) == 2:
            lat, lon = coord
        else:
            nums = NUM.findall(str(coord or ""))
            if len(nums) != 2:
                continue
            lat, lon = map(float, nums)
        yield {"source": "BOLD", "record_id": r.get("processid", ""),
               "taxon": r.get("species", ""), "country": r.get("country/ocean", ""),
               "host": r.get("associated_taxa", ""), "lat": float(lat), "lon": float(lon)}


def gbif_points():
    for path in glob.glob("data/gbif_*.csv"):
        with open(path) as f:
            for r in csv.DictReader(f):
                if r["decimalLatitude"] and r["decimalLongitude"]:
                    yield {"source": "GBIF", "record_id": r["gbifID"],
                           "taxon": r["scientificName"], "country": r["country"],
                           "host": r["associatedTaxa"], "lat": float(r["decimalLatitude"]),
                           "lon": float(r["decimalLongitude"])}


def main():
    points = list(genbank_points()) + list(bold_points()) + list(gbif_points())
    print(f"{len(points)} records with coordinates")
    if not points:
        return

    tmp = tempfile.mkdtemp()
    zpath = os.path.join(tmp, "wc.zip")
    print("Downloading", URL)
    urllib.request.urlretrieve(URL, zpath)
    with zipfile.ZipFile(zpath) as z:
        z.extractall(tmp)

    coords = [(p["lon"], p["lat"]) for p in points]
    for i in range(1, 20):
        tif = glob.glob(os.path.join(tmp, "**", f"wc2.1_{RES}_bio_{i}.tif"), recursive=True)[0]
        with rasterio.open(tif) as src:
            nodata = src.nodata
            for p, val in zip(points, src.sample(coords)):
                v = float(val[0])
                p[f"bio{i}"] = "" if nodata is not None and v == nodata else round(v, 3)

    out = f"data/climate_worldclim_{RES}.csv"
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(points[0].keys()))
        w.writeheader()
        w.writerows(points)
    print(f"Saved {out}")


if __name__ == "__main__":
    main()
