#! /usr/bin/env python3
#
# Get all protein names, gene names, and gene ontology terms from a given organism from ncbi
#
# Usage:
#   ncbi.py <organism_id>
#
# Example:
#   ncbi.py 160488
#
#   Output:
#   pputida160488_protein_names_only.txt
#   pputida160488_gene_names_only.txt
#   pputida160488_gene_ontology_terms_only.txt
#

import requests, argparse, json, time

parser = argparse.ArgumentParser(description="Get all protein names, gene names, and gene ontology terms from a given organism from ncbi")
parser.add_argument("taxonomy_id", type=int, help="The NCBI taxonomy ID (e.g. 9606 for human)")
args = parser.parse_args()

# Helper to call NCBI eUtils
def ncbi_request(url: str, params: dict, method: str = "get") -> dict:
    for attempt in range(3):  # simple retry
        if method == "post":
            resp = requests.post(url, data=params, timeout=30)
        else:
            resp = requests.get(url, params=params, timeout=30)
        if resp.status_code == 200:
            return resp.json()
        time.sleep(1)
    resp.raise_for_status()


# 1. Use ESearch to get all Gene IDs for the organism
esearch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
esearch_params = {
    "db": "gene",
    "term": f"txid{args.taxonomy_id}[Organism:exp]",
    "retmax": 100000,
    "retmode": "json",
}
esearch_data = ncbi_request(esearch_url, esearch_params, "get")
gene_ids = esearch_data.get("esearchresult", {}).get("idlist", [])

if not gene_ids:
    print("No genes found for taxonomy ID", args.taxonomy_id)
    exit()

# 2. Fetch summaries in batches (ESummary limit ~500 IDs per call)
esummary_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
records = []
batch_size = 500

for i in range(0, len(gene_ids), batch_size):
    batch_ids = ",".join(gene_ids[i : i + batch_size])
    esummary_params = {
        "db": "gene",
        "id": batch_ids,
        "retmode": "json",
    }
    data = ncbi_request(esummary_url, esummary_params, "post")
    result = data.get("result", {})
    # Each gene ID key maps to a dict; skip 'uids'
    for gid in result.get("uids", []):
        item = result.get(gid, {})
        records.append(item)

# Convert to DataFrame and write CSV with all columns
import pandas as pd

if not records:
    print("No summaries retrieved for taxonomy", args.taxonomy_id)
    exit()

df = pd.json_normalize(records)

out_file = f"gene_records_{args.taxonomy_id}.csv"
df.to_csv(out_file, index=False)

print(len(df), "gene records saved to", out_file)