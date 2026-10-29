"""Download nucleotide records from NCBI GenBank as FASTA and GenBank (full) files.

Usage: python fetch_ncbi.py "<NCBI query>" <output_name> [output_dir]
output_dir defaults to data; use a subfolder (e.g. data/its2) for non-COI markers
so the COI merge step does not read them.
Needs NCBI_EMAIL in the environment; NCBI_API_KEY is optional (raises the rate limit).
"""
import os
import sys
import time

from Bio import Entrez

BATCH = 500


def main():
    if len(sys.argv) not in (3, 4):
        sys.exit('Usage: python fetch_ncbi.py "<query>" <output_name> [output_dir]')
    query, name = sys.argv[1], sys.argv[2]
    out_dir = sys.argv[3] if len(sys.argv) == 4 else "data"

    email = (os.environ.get("NCBI_EMAIL") or "").strip()
    if not email:
        sys.exit("ERROR: the NCBI_EMAIL secret is empty or missing. "
                 "Add it under Settings > Secrets and variables > Actions.")
    Entrez.email = email
    api_key = (os.environ.get("NCBI_API_KEY") or "").strip()
    print(f"Email set: yes. API key set: {'yes' if api_key else 'no'}")

    result = None
    for key_value in ([api_key, None] if api_key else [None]):
        Entrez.api_key = key_value
        try:
            handle = Entrez.esearch(db="nuccore", term=query, usehistory="y", retmax=0)
            result = Entrez.read(handle)
            handle.close()
            break
        except Exception as err:
            body = ""
            if hasattr(err, "read"):
                try:
                    body = err.read().decode("utf-8", "replace")[:300]
                except Exception:
                    pass
            print(f"NCBI search failed ({'with' if key_value else 'without'} API key): {err} {body}")
            if key_value:
                print("Retrying without the API key. Check that NCBI_API_KEY was pasted "
                      "correctly, with no spaces.")
    if result is None:
        sys.exit("ERROR: NCBI search failed. See messages above.")
    api_key = Entrez.api_key
    pause = 0.12 if api_key else 0.4
    count = int(result["Count"])
    webenv, key = result["WebEnv"], result["QueryKey"]
    print(f"{count} records found for: {query}")
    if count == 0:
        return

    os.makedirs(out_dir, exist_ok=True)
    for rettype, ext in (("fasta", "fasta"), ("gb", "gb")):
        path = os.path.join(out_dir, f"{name}.{ext}")
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

    with open(os.path.join(out_dir, f"{name}_query.txt"), "w") as f:
        f.write(f"Query: {query}\nRecords: {count}\nDownloaded: {time.strftime('%Y-%m-%d %H:%M UTC', time.gmtime())}\n")


if __name__ == "__main__":
    main()
