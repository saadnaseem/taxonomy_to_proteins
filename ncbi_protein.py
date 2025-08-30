#!/usr/bin/env python3
"""
Fetch protein information by NCBI Taxonomy ID using E-utilities.

Outputs TSV with:
  accession, protein_name, gene_id, gene_symbol, locus_tag, protein_sequence

Usage:
  python ncbi_protein.py --taxid 562 --email you@lab.org --api_key YOUR_KEY --out proteins.tsv
  # E. coli = 562, Human = 9606, Pseudomonas putida = 160488
"""

import argparse, time, math, re
from collections import defaultdict
import requests
import pandas as pd

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
THROTTLE_NO_KEY = 0.34   # ~3 rps
THROTTLE_WITH_KEY = 0.11 # ~9 rps headroom for network jitter

def throttle(sleep):
    time.sleep(sleep)

def polite_params(tool, email=None, api_key=None):
    p = {"tool": tool}
    if email: p["email"] = email
    if api_key: p["api_key"] = api_key
    return p

def esearch_protein_ids_by_taxon(taxid, email=None, api_key=None, tool="pp_protein"):
    """Return all protein UIDs for txid via ESearch (usehistory enabled)."""
    sleep = THROTTLE_WITH_KEY if api_key else THROTTLE_NO_KEY

    params = {
        "db": "protein",
        "term": f"txid{taxid}[Organism:exp]",
        "retmode": "json",
        "usehistory": "y",
        "retmax": 0  # we only need WebEnv + count
    }
    params.update(polite_params(tool, email, api_key))
    throttle(sleep)
    r = requests.get(f"{EUTILS}/esearch.fcgi", params=params, timeout=60)
    r.raise_for_status()
    j = r.json()["esearchresult"]
    count = int(j.get("count", 0))
    webenv = j.get("webenv")
    query_key = j.get("querykey")
    if count == 0:
        return [], webenv, query_key
    print(f"[DEBUG] ESearch reports {count} protein records for taxid {taxid}")
    return count, webenv, query_key

def esummary_protein_meta_history(count, webenv, query_key, email=None, api_key=None, tool="pp_protein"):
    """ESummary using WebEnv history to avoid huge ID lists."""
    sleep = THROTTLE_WITH_KEY if api_key else THROTTLE_NO_KEY
    meta = {}
    retmax = 500
    for retstart in range(0, count, retmax):
        params = {
            "db": "protein",
            "query_key": query_key,
            "WebEnv": webenv,
            "retstart": retstart,
            "retmax": retmax,
            "retmode": "json"
        }
        params.update(polite_params(tool, email, api_key))
        throttle(sleep)
        r = requests.post(f"{EUTILS}/esummary.fcgi", data=params, timeout=90)
        r.raise_for_status()
        res = r.json().get("result", {})
        for k, v in res.items():
            if k == "uids": continue
            acc = v.get("caption")
            title = v.get("title")
            meta[k] = {"accession": acc, "protein_name": title}
    return meta

def elink_protein_to_gene_map(uids, email=None, api_key=None, tool="pp_protein"):
    """ELink: protein -> gene (GeneID list) with by-id mapping preserved."""
    # We don't specify linkname; Entrez supports protein->gene linking and ELink preserves by-id in JSON.
    # See ELink docs and sample apps. 
    sleep = THROTTLE_WITH_KEY if api_key else THROTTLE_NO_KEY
    prot2gene = defaultdict(list)

    for i in range(0, len(uids), 500):
        chunk = uids[i:i+500]
        params = {
            "dbfrom": "protein",
            "db": "gene",
            "id": ",".join(chunk),
            "retmode": "json",
            "cmd": "neighbor_history"  # keep mapping; JSON returns per-uid linksets
        }
        params.update(polite_params(tool, email, api_key))
        throttle(sleep)
        r = requests.get(f"{EUTILS}/elink.fcgi", params=params, timeout=90)
        r.raise_for_status()
        j = r.json()
        for linkset in j.get("linksets", []):
            src_ids = linkset.get("ids", [])
            if not src_ids:
                continue
            src_uid = src_ids[0]  # protein uid
            for ldb in (linkset.get("linksetdbs") or []):
                # Accept any protein->gene link name (e.g., protein_gene)
                if "gene" in ldb.get("linkname", ""):
                    prot2gene[src_uid].extend(ldb.get("links", []))
    return prot2gene

def efetch_gene_symbol_locus(gene_ids, email=None, api_key=None, tool="pp_protein"):
    """EFetch gene (XML) -> symbol + locus_tag per GeneID."""
    sleep = THROTTLE_WITH_KEY if api_key else THROTTLE_NO_KEY
    out = {}
    import xml.etree.ElementTree as ET

    for i in range(0, len(gene_ids), 500):
        chunk = gene_ids[i:i+500]
        params = {
            "db": "gene",
            "id": ",".join(chunk),
            "retmode": "xml"
        }
        params.update(polite_params(tool, email, api_key))
        throttle(sleep)
        r = requests.get(f"{EUTILS}/efetch.fcgi", params=params, timeout=120)
        r.raise_for_status()
        root = ET.fromstring(r.text)
        for eg in root.findall(".//Entrezgene"):
            gid = eg.findtext(".//Gene-track/Gene-track_geneid")
            symbol = eg.findtext(".//Gene-ref/Gene-ref_locus")
            locus_tag = eg.findtext(".//Gene-ref/Gene-ref_locus-tag")
            out[gid] = {"gene_symbol": symbol, "locus_tag": locus_tag}
    return out

def efetch_protein_fasta(uids, email=None, api_key=None, tool="pp_protein"):
    """EFetch protein FASTA -> AA sequence per protein UID."""
    sleep = THROTTLE_WITH_KEY if api_key else THROTTLE_NO_KEY
    seqs = {}
    for i in range(0, len(uids), 500):
        chunk = uids[i:i+500]
        params = {
            "db": "protein",
            "id": ",".join(chunk),
            "rettype": "fasta",
            "retmode": "text"
        }
        params.update(polite_params(tool, email, api_key))
        throttle(sleep)
        r = requests.get(f"{EUTILS}/efetch.fcgi", params=params, timeout=180)
        r.raise_for_status()
        # parse FASTA
        cur = None; buf = []
        for line in r.text.splitlines():
            if line.startswith(">"):
                if cur is not None:
                    seqs[cur] = "".join(buf)
                cur = line[1:].strip()  # header
                buf = []
            else:
                buf.append(line.strip())
        if cur is not None:
            seqs[cur] = "".join(buf)
    # map protein UID to sequence by extracting "gi|uid" or "uid" is NOT guaranteed in header now,
    # so we do a second pass: ESummary gave us accession; we'll index by accession instead.
    # Most FASTA headers start with the accession and version.
    acc_to_seq = {}
    acc_re = re.compile(r"^([A-Z]{2,}_[0-9]+\.[0-9]+)\s")
    for hdr, seq in seqs.items():
        m = acc_re.search(hdr)
        if m:
            acc_to_seq[m.group(1)] = seq
    return acc_to_seq

def main():
    ap = argparse.ArgumentParser(description="Retrieve protein info by NCBI Taxonomy ID (E-utilities).")
    ap.add_argument("--taxid", required=True, help="NCBI Taxonomy ID (e.g., 562 for E. coli, 9606 for human)")
    ap.add_argument("--email", default=None, help="Your contact email for NCBI (recommended).")
    ap.add_argument("--api_key", default=None, help="NCBI E-utilities API key (optional but faster).")
    ap.add_argument("--out", default="proteins_by_taxon.tsv", help="Output TSV file.")
    args = ap.parse_args()

    print(f"[INFO] Starting protein retrieval for taxonomy ID: {args.taxid}")
    print(f"[INFO] Known working taxonomy IDs: 562 (E. coli), 9606 (Human), 160488 (P. putida)")

    try:
        # 1) Get all protein UIDs for the taxon
        count, webenv, qk = esearch_protein_ids_by_taxon(args.taxid, args.email, args.api_key)
        
        if count == 0:
            pd.DataFrame(columns=["accession","protein_name","gene_id","gene_symbol","locus_tag","protein_sequence"]).to_csv(args.out, sep="\t", index=False)
            print(f"[WARNING] No protein records found for taxid {args.taxid}. Wrote empty file {args.out}")
            return
            
        print(f"[SUCCESS] Found {count} protein records for txid {args.taxid}")

        # 2) Protein metadata (accession, title)
        pmeta = esummary_protein_meta_history(count, webenv, qk, args.email, args.api_key)
        # 3) Map each protein -> GeneID(s)
        p2g = elink_protein_to_gene_map(pmeta.keys(), args.email, args.api_key)
        all_gene_ids = sorted({gid for lst in p2g.values() for gid in lst})
        # 4) Gene symbol + locus_tag
        gmeta = efetch_gene_symbol_locus(all_gene_ids, args.email, args.api_key)
        # 5) Protein sequences (FASTA)
        acc_to_seq = efetch_protein_fasta(pmeta.keys(), args.email, args.api_key)

        # Build rows
        rows = []
        for puid, meta in pmeta.items():
            acc = meta.get("accession")
            pname = meta.get("protein_name")
            seq = None
            # try accession.version first; if ESummary caption lacks version, accept without version
            if acc and acc in acc_to_seq:
                seq = acc_to_seq[acc]
            else:
                # try without version
                if acc:
                    core = acc.split(".")[0]
                    # find any acc_to_seq key that starts with core + "."
                    for k in acc_to_seq:
                        if k.startswith(core + "."):
                            seq = acc_to_seq[k]; break

            # Gene mapping (there can be multiple genes; list the first to keep table tidy)
            gene_ids = p2g.get(puid, [])
            gid = gene_ids[0] if gene_ids else None
            gsym = gmeta.get(gid, {}).get("gene_symbol") if gid else None
            ltag = gmeta.get(gid, {}).get("locus_tag") if gid else None

            rows.append({
                "accession": acc,
                "protein_name": pname,
                "gene_id": gid,
                "gene_symbol": gsym,
                "locus_tag": ltag,
                "protein_sequence": seq
            })

        df = pd.DataFrame(rows)
        # Optional: drop rows with no sequence (keep if you want metadata-only)
        # df = df[df["protein_sequence"].notna()]
        df.to_csv(args.out, sep="\t", index=False)
        print(f"[ok] Wrote {args.out} with {len(df)} rows")

    except Exception as e:
        print(f"[FATAL ERROR] Script failed: {e}")
        print(f"[SUGGESTION] Try with a known working taxonomy ID like 562 (E. coli)")
        raise

if __name__ == "__main__":
    main()