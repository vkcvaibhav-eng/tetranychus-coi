"""Download nucleotide records from NCBI GenBank as FASTA and GenBank (full) files.

Usage: python fetch_ncbi.py "<NCBI query>" <output_name>
Needs NCBI_EMAIL in the environment; NCBI_API_KEY is optional (raises the rate limit).
"""
import os
import sys
import time

from Bio import Entrez

BATCH = 500


def main():
    if len(sys.argv) != 3:
        sys.exit('Usage: python fetch_ncbi.py "<query>" <output_name>')
    query, name = sys.argv[1], sys.argv[2]

    Entrez.email = os.environ.get("NCBI_EMAIL") or sys.exit("Set NCBI_EMAIL")
    api_key = os.environ.get("NCBI_API_KEY")
    if api_key:
        Entrez.api_key = api_key
    pause = 0.12 if api_key else 0.4

    handle = Entrez.esearch(db="nuccore", term=query, usehistory="y", retmax=0)
    result = Entrez.read(handle)
    handle.close()
    count = int(result["Count"])
    webenv, key = result["WebEnv"], result["QueryKey"]
    print(f"{count} records found for: {query}")
    if count == 0:
        return

    os.makedirs("data", exist_ok=True)
    for rettype, ext in (("fasta", "fasta"), ("gb", "gb")):
        path = os.path.join("data", f"{name}.{ext}")
        with open(path, "w") as out:
            for start in range(0, count, BATCH):
                for attempt in range(3):
                    try:
                        h = Entrez.efetch(db="nuccore", rettype=rettype, retmode="text",
                                          retstart=start, retmax=BATCH,
                                          webenv=webenv, query_key=key)
                        out.write(h.read())
                        h.close()
                        break
                    except Exception as err:
                        print(f"Retry {attempt + 1} at {start}: {err}")
                        time.sleep(5)
                else:
                    sys.exit(f"Failed at record {start}")
                print(f"{rettype}: {min(start + BATCH, count)}/{count}")
                time.sleep(pause)
        print(f"Saved {path}")

    with open(os.path.join("data", f"{name}_query.txt"), "w") as f:
        f.write(f"Query: {query}\nRecords: {count}\nDownloaded: {time.strftime('%Y-%m-%d %H:%M UTC', time.gmtime())}\n")


if __name__ == "__main__":
    main()
