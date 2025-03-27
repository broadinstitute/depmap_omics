# pip install gcsfs
# pip install fsspec
# pip install biomart
#fissfc entity_tsv -w depmap-omics-rna -p broad-firecloud-ccle -t sample | sed "s/\[\"//g" | sed 's/\"\]//g' > terra_data_table_rnaseq_25q2.tsv

from collections import defaultdict
from taigapy import create_taiga_client_v3
import argparse
import pandas as pd
import numpy as np
from concurrent.futures import ThreadPoolExecutor, as_completed

parser = argparse.ArgumentParser()
parser.add_argument("--date_on_or_before", type=str, required=True)

args = parser.parse_args()
release_date = args.date_on_or_before

terra_samples = pd.read_table("/localstuff/terra_data_table_rnaseq_25q2.tsv")
samples_to_process_all = pd.read_csv("/localstuff/merged_table.20250325_2319selcols.csv")
# Ensure the column is datetime
samples_to_process_all["internal_release_date"] = pd.to_datetime(samples_to_process_all["internal_release_date"], errors='coerce')

samples_to_process_release = samples_to_process_all.loc[
	(samples_to_process_all["is_default_entry"]) & 
	(samples_to_process_all["datatype"] == "rna") & 
	(samples_to_process_all["internal_release_date"] != "None") &  # Ensure the date is not NA
	(samples_to_process_all["internal_release_date"] <= pd.Timestamp('2025-05-01'))
].copy()



samples = pd.merge(terra_samples, samples_to_process_release, left_on ="entity:sample_id", right_on="sequencing_id", how="inner")
tc = create_taiga_client_v3()

hgnc_table = tc.get("hgnc-gene-table-e250.3/hgnc_complete_set")
hgnc_table = hgnc_table[
		(~hgnc_table["entrez_id"].isna())
	]
hgnc_table = hgnc_table[
		(~hgnc_table["symbol"].isna())
	]
hgnc_table["hugo_entrez"] = (
		hgnc_table["symbol"].astype(str)
		+ " ("
		+ hgnc_table["entrez_id"].astype("Int64").astype(str)
		+ ")"
	)


ens_to_hugo_entrez = dict(zip(hgnc_table["ensembl_gene_id"], hgnc_table["hugo_entrez"]))
ens_gene_biotype = dict(zip(hgnc_table["ensembl_gene_id"], hgnc_table["locus_group"]))

tpm_dict_human_all_genes = {}
tpm_dict_human_pc_genes = {}
tpm_dict_virus = {}

counts_dict_human_all_genes = {}
counts_dict_human_pc_genes = {}
counts_dict_virus = {}

profile_cds_dict = dict(zip(samples["sequencing_id"],samples["profile_id"]))
model_cds_dict = dict(zip(samples["sequencing_id"],samples["model_id"]))
cell_line_cds_dict = dict(zip(samples["sequencing_id"],samples["stripped_cell_line_name"]))
code_cds_dict = dict(zip(samples["sequencing_id"],samples["depmap_code"]))
lineage_cds_dict = dict(zip(samples["sequencing_id"],samples["lineage"]))

for sample_id, sample_data in list(samples.iterrows())[1:100]:
	sample_key = sample_data["entity:sample_id"]
	if pd.notna(sample_data["quant_genes"]):
		sample_key = sample_key
		print(sample_id, sample_key)
		df = pd.read_table(sample_data["quant_genes"])
		df_human_all_genes = df[df["Name"].str.startswith("ENSG")].copy()
		profile_id = profile_cds_dict[sample_key]
		cellline = cell_line_cds_dict[sample_key]
		model_id = model_cds_dict[sample_key]
		depmap_code = code_cds_dict[sample_key]
		lineage = lineage_cds_dict[sample_key]
		sample_profile_model = model_id + "@" + cellline + "@" + depmap_code + "@" + lineage + "@" +  profile_id + "@" +sample_key 
		ens_no_suffix = df_human_all_genes["Name"].str.split(".", n=1).str[0]
		df_human_all_genes["ens_no_suffix"] = ens_no_suffix
		df_human_all_genes["hugo_entrez"] = ens_no_suffix.map(ens_to_hugo_entrez).fillna(df_human_all_genes["Name"])
		df_human_all_genes["gene_biotype"] = ens_no_suffix.map(ens_gene_biotype)
		df_human_all_genes.set_index("hugo_entrez", inplace=True)
		df_human_pc_genes = df_human_all_genes[df_human_all_genes["gene_biotype"] == "protein-coding gene"]
		tpm_dict_human_all_genes[sample_profile_model] = df_human_all_genes["TPM"]
		tpm_dict_human_pc_genes[sample_profile_model] = df_human_pc_genes["TPM"]
		counts_dict_human_all_genes[sample_profile_model] = df_human_all_genes["NumReads"].round().astype(int)
		counts_dict_human_pc_genes[sample_profile_model] = df_human_pc_genes["NumReads"].round().astype(int)
		df_virus = df[~df["Name"].str.startswith("ENSG")].copy()
		df_virus.set_index("Name", inplace=True)
		tpm_dict_virus[sample_profile_model] = df_virus["TPM"]
		counts_dict_virus[sample_profile_model] = df_virus["NumReads"].round().astype(int)
	else:
		print(sample_key + " does not have salmon output")

tpm_df_human_all_genes = pd.DataFrame(tpm_dict_human_all_genes).fillna(0)
tpm_df_human_pc_genes = pd.DataFrame(tpm_dict_human_pc_genes).fillna(0)
tpm_df_virus = pd.DataFrame(tpm_dict_virus).fillna(0)

tpm_df_human_all_genes = tpm_df_human_all_genes.T
log2tpm_df_human_all_genes = np.log2(tpm_df_human_all_genes + 1)
tpm_df_human_pc_genes = tpm_df_human_pc_genes.T
log2tpm_df_human_pc_genes = np.log2(tpm_df_human_pc_genes + 1)
tpm_df_virus = tpm_df_virus.T
log2tpm_df_virus = np.log2(tpm_df_virus + 1)

counts_df_human_all_genes = pd.DataFrame(counts_dict_human_all_genes).fillna(0)
counts_df_human_pc_genes = pd.DataFrame(counts_dict_human_pc_genes).fillna(0)
counts_df_virus = pd.DataFrame(counts_dict_virus).fillna(0)

counts_df_human_all_genes = counts_df_human_all_genes.T
counts_df_human_pc_genes = counts_df_human_pc_genes.T
counts_df_virus = counts_df_virus.T


tpm_df_human_all_genes.to_csv("tpm_df_human_all_genes.csv")
log2tpm_df_human_all_genes.to_csv("log2tpm_df_human_all_genes.csv")
tpm_df_human_pc_genes.to_csv("tpm_df_human_pc_genes.csv")
log2tpm_df_human_pc_genes.to_csv("/localstuff/log2tpm_df_human_pc_genes.csv")
tpm_df_virus.to_csv("tpm_df_virus.csv")
log2tpm_df_virus.to_csv("log2tpm_df_virus.csv")



def process_sample(sample_data):
	result = defaultdict(dict)
	if pd.notna(sample_data["quant_genes"]) and sample_data["entity:sample_id"] in profile_cds_dict.keys():
		#print(sample_data["entity:sample_id"])
		df = pd.read_table(sample_data["quant_genes"])
		df_human_all_genes = df[df["Name"].str.startswith("ENSG")].copy()
		profile_id = profile_cds_dict[sample_data["entity:sample_id"]]
		cellline = cell_line_cds_dict[sample_data["entity:sample_id"]]
		model_id = model_cds_dict[sample_data["entity:sample_id"]]
		depmap_code = code_cds_dict[sample_data["entity:sample_id"]]
		lineage = lineage_cds_dict[sample_data["entity:sample_id"]]
		sample_profile_model = f"{model_id}@{cellline}@{depmap_code}@{lineage}@{profile_id}@{sample_data['entity:sample_id']}"
		df_human_all_genes["ens_no_suffix"] = df_human_all_genes["Name"].str.split(".").str[0]
		df_human_all_genes["hugo_entrez"] = df_human_all_genes["ens_no_suffix"].map(ens_to_hugo_entrez).fillna(df_human_all_genes["Name"])
		df_human_all_genes["gene_biotype"] = df_human_all_genes["ens_no_suffix"].map(ens_gene_biotype)
		df_human_all_genes.set_index("hugo_entrez", inplace=True)
		df_human_pc_genes = df_human_all_genes[df_human_all_genes["gene_biotype"] == "protein-coding gene"]
		result['tpm_dict_human_all_genes'][sample_profile_model] = df_human_all_genes["TPM"]
		result['tpm_dict_human_pc_genes'][sample_profile_model] = df_human_pc_genes["TPM"]
		result['counts_dict_human_all_genes'][sample_profile_model] = df_human_all_genes["NumReads"].round().astype(int)
		result['counts_dict_human_pc_genes'][sample_profile_model] = df_human_pc_genes["NumReads"].round().astype(int)
		df_virus = df[~df["Name"].str.startswith("ENSG")].copy()
		df_virus.set_index("Name", inplace=True)
		result['tpm_dict_virus'][sample_profile_model] = df_virus["TPM"]
		result['counts_dict_virus'][sample_profile_model] = df_virus["NumReads"].round().astype(int)
	return result

def merge_results(results):
	merged = defaultdict(dict)
	for result in results:
		for key, sub_dict in result.items():
			merged[key].update(sub_dict)
	return merged

# Parallel Execution
results = []
starttime = time.time()
with ThreadPoolExecutor() as executor:
	futures = [executor.submit(process_sample, sample_data) for _, sample_data in samples.iterrows()]
	for future in as_completed(futures):
		results.append(future.result())

# Final Merge
merged_results = merge_results(results)
endtime = time.time()
merged_results.to_csv("/localstuff/merged_results.csv")
print(f"Time taken: {time.time() - starttime}")