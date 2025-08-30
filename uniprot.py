#! /usr/bin/env python3
#
# Get all protein names from a given organism from UniProt
#
# Usage:
#   uniprot.py <organism_id>
#
# Example:
#   uniprot.py 160488
#
# Output:
#   pputida160488_protein_names_only.txt
#

import requests, gzip, io, pandas as pd, argparse

parser = argparse.ArgumentParser(description="Get all protein names from a given organism from UniProt")
parser.add_argument("organism_id", type=int, help="The organism ID from UniProt")
args = parser.parse_args()

url = ("https://rest.uniprot.org/uniprotkb/stream"
       f"?query=organism_id:{args.organism_id}"
       "&fields=gene_oln, accession,protein_name,	protein_families,id,gene_names,sequence,cc_catalytic_activity,cc_cofactor,ec,	ft_binding,	cc_activity_regulation,	ft_act_site,absorption,cc_function,	kinetics,cc_pathway,redox_potential,rhea,	ft_site,reviewed,cc_interaction,cc_subunit,go_p,	go_c,go,go_f,go_id,ft_intramem,cc_subcellular_location,ft_topo_dom,	ft_transmem,lit_pubmed_id,ft_motif,ft_region,	ft_repeat,ft_zn_fing  "
       "&format=tsv&compressed=true")

# https://www.uniprot.org/help/return_fields 

r = requests.get(url, stream=True)
r.raise_for_status()

# Load TSV directly from gzipped bytes into a DataFrame
df = pd.read_csv(io.BytesIO(r.content), compression="gzip", sep="\t")

#df to csv
df.to_csv(f"protein_only_{args.organism_id}.csv", index=False)

# If you only want the names:
names = df["Protein names"].dropna().unique()
print(len(names), "unique names")
# Save just the names if desired
pd.Series(names).to_csv(f"protein_names_only_{args.organism_id}.txt", index=False)

