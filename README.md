# UniProt / NCBI Retrieval Toolkit

This mini-toolkit lets you download protein or gene information for a given organism (taxonomy ID) from two major biological databases.

## 1. Required Files

| Script | Purpose |
|--------|---------|
| `uniprot.py` | Fetches **protein-centric** annotations from UniProt REST API (TSV/CSV + name list). |
| `ncbi.py` | Fetches **gene-centric** metadata from NCBI (ESearch / ESummary) into a CSV file. |
| `ncbi_protein.py` | Fetches **protein-centric** data from NCBI (ESearch / ESummary / EFetch) and writes a TSV with sequences and gene mappings. |

> Make sure these three Python files are in the same directory (this repo’s root). No other project files are required.

## 2. Prerequisites

* Python 3.8+
* `requests` and `pandas` (install with `pip install requests pandas`)

Optional but recommended:
* An NCBI **API key** (free – increases rate-limit to 10 req/s).

## 3. Usage

```bash
# activate your virtualenv / conda env first

# 3.1 UniProt – download protein table for an organism
python uniprot.py <organism_taxid>
# Example: Pseudomonas putida (160488)
python uniprot.py 160488
# → outputs:
#   - protein_only_<taxid>.csv         (full table)
#   - protein_names_only_<taxid>.txt   (unique names list)

# 3.2 NCBI – gene metadata for an organism
python ncbi.py <taxonomy_id>
# Example: Human (9606)
python ncbi.py 9606
# → outputs gene_records_<taxid>.csv

# 3.3 NCBI – protein + gene mapping + sequences (advanced)
python ncbi_protein.py --taxid <taxonomy_id> [--email you@lab.org] [--api_key YOUR_NCBI_KEY] [--out output.tsv]
# Example: E. coli (562)
python ncbi_protein.py --taxid 562 --email you@lab.org --api_key ABC123 --out ecoli_proteins.tsv
```

### Flags Explained (ncbi_protein.py)

* `--taxid`   NCBI taxonomy ID (required)
* `--email`   Your contact email (optional but recommended – NCBI courtesy)
* `--api_key` NCBI API key to raise rate-limit (optional)
* `--out`     Output TSV filename (default `proteins_by_taxon.tsv`)

### Output Columns (ncbi_protein.py)

* `accession` – protein accession (with version)
* `protein_name` – description / defline
* `gene_id` – linked NCBI Gene ID
* `gene_symbol` – gene symbol (if available)
* `locus_tag` – locus tag (if available)
* `protein_sequence` – amino-acid sequence

## 4. Notes / Troubleshooting

* NCBI occasionally returns **HTTP 500** for very large queries. `ncbi_protein.py` now uses ESearch history + paginated ESummary POST requests to avoid this.
* If you hit rate-limits, lower `--taxid` scope (e.g., subspecies) or supply an API key.

## 5. Testing

A minimal smoke-test (requires internet):

```bash
python -m unittest tests/test_ncbi_scripts.py
```

See `tests/test_ncbi_scripts.py` for simple assertions that the scripts run and write non-empty tables for small model organisms (E. coli 562).

---
Feel free to open issues or PRs for improvements.
