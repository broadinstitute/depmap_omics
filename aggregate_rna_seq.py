# pip install gcsfs
# pip install fsspec
# pip install biomart
#fissfc entity_tsv -w depmap-omics-rna -p broad-firecloud-ccle -t sample | sed "s/\[\"//g" | sed 's/\"\]//g' > terra_data_table_rnaseq_25q2.tsv

from collections import defaultdict
from taigapy import create_taiga_client_v3
from taigapy.client_v3 import UploadedFile, LocalFormat
import argparse
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("--terra_table", type=str, required=True, help='Terra table to use')
parser.add_argument("--sample_metadata", type=str, required=True, help='Sample metadata file to use')

args = parser.parse_args()
terra_samples_file = args.terra_table # test with terra_samples_file="/localstuff/terra_data_table_rnaseq_25q2.tsv"
sample_metadata_file = args.sample_metadata # test with sample_metadata_file="/localstuff/2025-05-01-master-mapping-table_v4-internal-release-date-2025-05-01-master-mapping-table.csv"
#export GOOGLE_APPLICATION_CREDENTIALS=/localstuff/depmap-omics-9764fbdbe040.json
terra_samples = pd.read_table(terra_samples_file)
samples_to_process_all = pd.read_csv(sample_metadata_file)

samples_to_process = samples_to_process_all.loc[(samples_to_process_all["datatype"] == "rna")]
samples = pd.merge(terra_samples, samples_to_process, left_on ="entity:sample_id", right_on="sequencing_id", how="inner")

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
hgnc_gene_biotype = dict(zip(hgnc_table["hugo_entrez"], hgnc_table["locus_group"]))

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
is_default_cds_dict = dict(zip(samples["sequencing_id"],samples["is_default_entry"]))

df_all_tpms = pd.DataFrame()
df_all_counts = pd.DataFrame()
df_all_lengths = pd.DataFrame()
# Establish the order of genes in output files (by alphabetical order) by anchoring on the first sample
geneorder = pd.read_table(samples.loc[0,"quant_genes"]).sort_values(by="Name")["Name"].reset_index(drop=True)
ens_no_suffix = geneorder.str.split(".", n=1).str[0]
ens_no_suffix.loc[geneorder.str.endswith("_PAR_Y")] = ens_no_suffix.loc[geneorder.str.endswith("_PAR_Y")] + "_PAR_Y"
gene_biotype = pd.Series(ens_no_suffix.map(ens_gene_biotype))
hgnc_names = pd.Series(ens_no_suffix.map(ens_to_hugo_entrez))
hgnc_names = hgnc_names.fillna(geneorder)
human_genes_all_indexes = geneorder[geneorder.str.startswith("ENSG")].index
human_genes_pc_indexes = geneorder[(geneorder.str.startswith("ENSG")) & (ens_no_suffix.map(ens_gene_biotype) == "protein-coding gene")].index
virus_genes_all_indexes = geneorder[~geneorder.str.startswith("ENSG")].index

all_tpms_list = []
all_counts_list = []
all_lengths_list = []
sample_labels = []
bigfusiontable = pd.DataFrame()

for sample_id, sample_data in list(samples.iterrows())[:20]:
	sample_key = sample_data["entity:sample_id"]
	print(sample_id)
	if pd.notna(sample_data["quant_genes"]):
		df = pd.read_table(sample_data["quant_genes"]).sort_values(by="Name").reset_index(drop=True)
		tpms = df["TPM"].values  # Use NumPy arrays for speed
		counts = df["NumReads"].round().astype(int).values
		lengths = df["EffectiveLength"].values
		genenames = df["Name"].values
		if (genenames != geneorder).any():
			raise ValueError("Gene order is not the same for " + sample_key)
		sample_labels.append(sample_key)
		# Append data to lists
		all_tpms_list.append(tpms)
		all_counts_list.append(counts)
		all_lengths_list.append(lengths)
	else:
		raise ValueError(sample_key + " does not have salmon output")
	if pd.notna(sample_data["fusions"]):
		fusiondf = pd.read_table(sample_data["fusions"])
		# Use STAR read output to get total reads for normalization later
		star_read_counts_df = pd.read_table(sample_data["reads_per_gene"], sep='\t',names=["colname","unstranded_counts","forward_stranded_counts","reverse_stranded_counts"])
		total_reads = pd.Series.sum(star_read_counts_df["unstranded_counts"][4:])
		if fusiondf.shape[0] > 0:
			#fusiondf["sample_id"] = sample_key
			fusiondf["profile_id"] = sample_data["profile_id"]
			fusiondf["model_id"] = sample_data["model_id"]
			fusiondf["cellline"] = sample_data["stripped_cell_line_name"]
			fusiondf["oncotree_code"] = sample_data["depmap_code"]
			fusiondf["lineage"] = sample_data["lineage"]
			fusiondf["total_reads"] = total_reads
			bigfusiontable = pd.concat([bigfusiontable, fusiondf], ignore_index=True)
		else:
			print("No fusions for " + sample_key)


# Canonicalize gene fusion names, making it easier to group by

bigfusiontable['Canonical_FusionName'] = bigfusiontable.apply(lambda row: '--'.join(sorted([row['#gene1'], row['gene2']])), axis=1)
# Calculate FFPM by summing (split_reads1+split_reads2)/total reads for sample and multiplying by 1e6
bigfusiontable["fusion_reads"] = bigfusiontable["split_reads1"] + bigfusiontable["split_reads2"] + bigfusiontable["discordant_mates"]
bigfusiontable["FFPM"] = 1e6*bigfusiontable["fusion_reads"]/bigfusiontable["total_reads"]

# Convert lists to DataFrame at the end
df_all_tpms = pd.DataFrame(np.column_stack(all_tpms_list), columns=sample_labels)
df_all_tpms.insert(0, 'Name', hgnc_names)
df_all_counts = pd.DataFrame(np.column_stack(all_counts_list), columns=sample_labels)
df_all_counts.insert(0, 'Name', hgnc_names)
df_all_lengths = pd.DataFrame(np.column_stack(all_lengths_list), columns=sample_labels)
df_all_lengths.insert(0, 'Name', hgnc_names)


df_all_tpms_cp = df_all_tpms.copy()
df_all_counts_cp = df_all_counts.copy()
df_all_lengths_cp = df_all_lengths.copy()

upload_files = []
df_dict = {
	"df_tpms": df_all_tpms_cp,
	"df_counts": df_all_counts_cp,
	"df_lengths": df_all_lengths_cp,
}
df_outputs = {
	"human_all_genes":human_genes_all_indexes,
	"human_pc_genes":human_genes_pc_indexes,
	"virus_genes":virus_genes_all_indexes
}
for thisdfname, thisdf in df_dict.items():
	thisdf["hgnc_name"] = hgnc_names
	thisdf["hgnc_name"] = thisdf["hgnc_name"].fillna(thisdf["Name"])
	print(thisdfname)
	thisdf.set_index("hgnc_name", inplace=True)
	thisdf = thisdf.drop(['Name'], axis = 1)
	if thisdfname == "df_tpms":
		thisdf = np.log2(thisdf + 1)
	thisdf = thisdf.T
	profile_id = thisdf.index.map(profile_cds_dict)
	cellline = thisdf.index.map(cell_line_cds_dict)
	model_id = thisdf.index.map(model_cds_dict)
	depmap_code = thisdf.index.map(code_cds_dict)
	lineage = thisdf.index.map(lineage_cds_dict)
	is_default_entry = thisdf.index.map(is_default_cds_dict)
	metadatadf = pd.DataFrame({"seqid":thisdf.index, "profile_id":profile_id, "is_default_entry":is_default_entry, "cellline":cellline, "model_id":model_id, "oncotree_code":depmap_code, "lineage":lineage})
	metadatadf.set_index("seqid", inplace=True)
	# Create an upload file for the full dataset, human and virus genes
	thisdf.to_parquet(thisdfname + ".parquet", engine="pyarrow", index=False)
	upload_files.append(UploadedFile(name=thisdfname, local_path=thisdfname + ".parquet", format=LocalFormat.PARQUET_TABLE))
	# Create an upload file for the human and virus genes
	for df_output_name, df_output_indexes in df_outputs.items():
		globals()[thisdfname + "_"+ df_output_name] = metadatadf.join(thisdf.iloc[:, df_output_indexes])
		globals()[thisdfname + "_"+ df_output_name].to_parquet(thisdfname + "_"+ df_output_name+".parquet", engine="pyarrow", index=False)
		upload_files.append(UploadedFile(name=thisdfname + "_"+ df_output_name, local_path=thisdfname + "_"+ df_output_name+".parquet", format=LocalFormat.PARQUET_TABLE))

bigfusiontable.rename(columns={"#gene1":"gene1", "total_reads":"Total Reads In Sample", "fusion_reads": "Total Reads Supporting Fusion"}, inplace=True)
fusion_output_columns = ['profile_id', 'model_id', 'cellline', 'oncotree_code', 'lineage', 
						 'Canonical_FusionName', 'gene1','gene2', 'Total Reads In Sample', 'Total Reads Supporting Fusion',  'FFPM', 
						 'confidence','split_reads1', 'split_reads2', 'discordant_mates', 
						 'strand1(gene/fusion)', 'strand2(gene/fusion)',  'reading_frame',
						 'breakpoint1', 'breakpoint2', 'site1', 'site2', 'type', 'coverage1', 'coverage2',
						 'tags', 'retained_protein_domains',
						 'closest_genomic_breakpoint1', 'closest_genomic_breakpoint2',
						 'gene_id1', 'gene_id2', 'transcript_id1', 'transcript_id2',
						 'direction1', 'direction2', 'filters', 'fusion_transcript',
						 'peptide_sequence']

bigfusiontable[fusion_output_columns].to_parquet("SR_fusion_table.parquet", engine="pyarrow", index=False) 
upload_files.append(UploadedFile(name="SR_fusion_table", local_path="SR_fusion_table.parquet", format=LocalFormat.PARQUET_TABLE))
tc = create_taiga_client_v3()
tc.create_dataset(name="25Q2 Expression Data Full Set",description="Testing upload of all expr files from new pipeline", files=upload_files)