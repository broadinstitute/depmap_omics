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

def str2bool(v):
    return v.lower() in ("yes", "true", "t", "1")

parser = argparse.ArgumentParser()
parser.add_argument("--terra_table", type=str, required=True, help='Terra table to use')
parser.add_argument("--sample_metadata", type=str, required=True, help='Sample metadata file to use')
parser.add_argument("--release_permaname", type=str, required=True, help='Permaname to use for release')
parser.add_argument("--stranded", type=str2bool, default=False, help="Whether the data is stranded")

args = parser.parse_args()
terra_table = args.terra_table # test with terra_samples_file="/localstuff/terra_data_table_rnaseq_25q2.tsv"
sample_metadata = args.sample_metadata # test with sample_metadata_file="/localstuff/2025-05-01-master-mapping-table_v4-internal-release-date-2025-05-01-master-mapping-table.csv"
release_date = args.release_permaname # Permaname to use for release
stranded = args.stranded # not used for any logic, just used to know whether to add "Stranded" suffix
if (stranded):
	stranded_suffix="Stranded"
	quant_genes_column = "quant_genes_auto"
	quant_trancripts_column = "quant_transcripts_auto"
else:
	stranded_suffix=""
	quant_genes_column = "quant_genes_iu"
	quant_trancripts_column = "quant_transcripts_iu"
#export GOOGLE_APPLICATION_CREDENTIALS=/localstuff/depmap-omics-9764fbdbe040.json
terra_samples = pd.read_table(terra_table)
terra_samples['quant_genes_auto'].fillna(terra_samples['quant_genes_iu'], inplace=True)
samples_to_process_all = pd.read_csv(sample_metadata)

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
is_default_cds_dict = dict(zip(samples["sequencing_id"],samples["is_default_entry"]))

df_all_tpms = pd.DataFrame()
df_all_counts = pd.DataFrame()
df_all_lengths = pd.DataFrame()
# Establish the order of genes in output files (by alphabetical order) by anchoring on the first sample
geneorder = pd.read_table(samples.loc[0,quant_genes_column]).sort_values(by="Name")["Name"].reset_index(drop=True)
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
all_raw_counts_list = []
sample_labels = []
bigfusiontable = pd.DataFrame()

for sample_id, sample_data in list(samples.iterrows()):
	sample_key = sample_data["entity:sample_id"]
	print(sample_id)
	star_read_counts_df = pd.read_table(sample_data["reads_per_gene"], sep='\t',names=["colname","unstranded_counts","forward_stranded_counts","reverse_stranded_counts"])
	all_raw_counts_list.append(star_read_counts_df["unstranded_counts"][4:])
	total_reads_for_sample = pd.Series.sum(star_read_counts_df["unstranded_counts"][4:])
	if pd.notna(sample_data[quant_genes_column]):
		df = pd.read_table(sample_data[quant_genes_column]).sort_values(by="Name").reset_index(drop=True)
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
		if fusiondf.shape[0] > 0:
			#fusiondf["sample_id"] = sample_key
			fusiondf["profile_id"] = sample_data["profile_id"]
			fusiondf["model_id"] = sample_data["model_id"]
			fusiondf["cellline"] = sample_data["stripped_cell_line_name"]
			fusiondf["oncotree_code"] = sample_data["depmap_code"]
			fusiondf["lineage"] = sample_data["lineage"]
			fusiondf["is_default_entry"] = sample_data["is_default_entry"]
			fusiondf["total_reads_in_sample"] = total_reads_for_sample
			bigfusiontable = pd.concat([bigfusiontable, fusiondf], ignore_index=True)
		else:
			print("No fusions for " + sample_key)

# Calculate FFPM by summing (split_reads1+split_reads2)/total reads for sample and multiplying by 1e6
bigfusiontable["total_reads_supporting_fusion"] = bigfusiontable["split_reads1"] + bigfusiontable["split_reads2"] + bigfusiontable["discordant_mates"]
bigfusiontable["FFPM"] = 1e6*bigfusiontable["total_reads_supporting_fusion"]/bigfusiontable["total_reads_in_sample"]
# Canonicalize gene fusion names, making it easier to group by
bigfusiontable['gene1_clean'] = bigfusiontable['#gene1'].str.split(',').str[0].str.extract(r'^([^\(]+)')[0]
bigfusiontable['gene2_clean'] = bigfusiontable['gene2'].str.split(',').str[0].str.extract(r'^([^\(]+)')[0]
bigfusiontable['Canonical_FusionName'] = bigfusiontable.apply(lambda row: '--'.join(sorted([row['gene1_clean'], row['gene2_clean']])), axis=1)
bigfusiontable['total_fusion_coverage'] = bigfusiontable['coverage1'] + bigfusiontable['coverage2']
bigfusiontable_defaults_only = bigfusiontable.loc[(bigfusiontable['is_default_entry']) & (bigfusiontable['confidence'] == "high") & (bigfusiontable['site1'] != 'intergenic') & (bigfusiontable['site2'] != 'intergenic')] 
fusion_grouped_by_sample_and_genes = bigfusiontable_defaults_only.groupby(['model_id', 'profile_id', 'Canonical_FusionName', 'gene1_clean', 'gene2_clean', 'strand1(gene/fusion)', 'strand2(gene/fusion)','reading_frame', 'total_reads_in_sample']).apply(
    lambda g: pd.Series({
		'split_reads1' : g['split_reads1'].sum(),
		'split_reads2' : g['split_reads2'].sum(),
		'discordant_mates' : g['discordant_mates'].sum(),
		'total_reads_supporting_fusion' : g['total_reads_supporting_fusion'].sum(),
		'total_fusion_coverage': g['total_fusion_coverage'].sum()})).reset_index()

fusion_grouped_by_sample_and_genes['FFPM'] = 1e6*fusion_grouped_by_sample_and_genes["total_reads_supporting_fusion"]/fusion_grouped_by_sample_and_genes["total_reads_in_sample"]
fusion_grouped_by_sample_and_genes.rename(columns={"gene1_clean":"gene1", "gene2_clean": "gene2"}, inplace=True)

selcols = [
'Canonical_FusionName',
'gene1',
'gene2',
'profile_id',
'model_id',
'total_reads_supporting_fusion',
'total_fusion_coverage',
'FFPM',
'split_reads1',
'split_reads2',
'discordant_mates',
'strand1(gene/fusion)',
'strand2(gene/fusion)',
'reading_frame']
fusion_by_model_df = fusion_grouped_by_sample_and_genes[selcols]
fusion_by_model_df.to_parquet("OmicsFusionFiltered.parquet", index=True)

# Convert lists to DataFrame at the end
df_all_tpms = pd.DataFrame(np.column_stack(all_tpms_list), columns=sample_labels)
df_all_tpms.insert(0, 'Name', hgnc_names)
df_all_counts = pd.DataFrame(np.column_stack(all_counts_list), columns=sample_labels)
df_all_counts.insert(0, 'Name', hgnc_names)
df_all_raw_counts = pd.DataFrame(np.column_stack(all_raw_counts_list), columns=sample_labels)
df_all_raw_counts.insert(0, 'Name', hgnc_names)
df_all_lengths = pd.DataFrame(np.column_stack(all_lengths_list), columns=sample_labels)
df_all_lengths.insert(0, 'Name', hgnc_names)


df_all_tpms_cp = df_all_tpms.copy()
df_all_counts_cp = df_all_counts.copy()
df_all_lengths_cp = df_all_lengths.copy()
df_all_raw_counts_cp = df_all_raw_counts.copy()

upload_files = []

df_dict = {
	"OmicsExpressionTPMLogp1": df_all_tpms_cp,
	"OmicsExpressionExpectedCount": df_all_counts_cp,
	"OmicsExpressionEffectiveLength": df_all_lengths_cp,
	"OmicsExpressionRawReadCount": df_all_raw_counts_cp
}
df_outputs = {
	"HumanAllGenes"+stranded_suffix:human_genes_all_indexes,
	"HumanProteinCodingGenes"+stranded_suffix:human_genes_pc_indexes,
	"VirusAllGenes"+stranded_suffix:virus_genes_all_indexes
}
for thisdfname, thisdf in df_dict.items():
	thisdf["hgnc_name"] = hgnc_names
	thisdf["hgnc_name"] = thisdf["hgnc_name"].fillna(thisdf["Name"])
	print(thisdfname)
	thisdf.set_index("hgnc_name", inplace=True)
	thisdf = thisdf.drop(['Name'], axis = 1)
	if thisdfname == "HumanProteinCodingGenes"+stranded_suffix:
		thisdf = np.log2(thisdf + 1)
	thisdf = thisdf.T
	profile_id = thisdf.index.map(profile_cds_dict)
	model_id = thisdf.index.map(model_cds_dict)
	is_default_entry = thisdf.index.map(is_default_cds_dict)
	metadatadf = pd.DataFrame({"seqid":thisdf.index, "profile_id":profile_id, "is_default_entry":is_default_entry, "model_id":model_id})
	metadatadf.set_index("seqid", inplace=True)
	# Create an upload file for the full dataset, human and virus genes
	#thisdf.to_parquet(thisdfname + ".parquet", engine="pyarrow", index=False)
	#upload_files.append(UploadedFile(name=thisdfname, local_path=thisdfname + ".parquet", format=LocalFormat.PARQUET_TABLE))
	# Create an upload file for the human and virus genes
	for df_output_name, df_output_indexes in df_outputs.items():
		globals()[thisdfname + df_output_name] = metadatadf.join(thisdf.iloc[:, df_output_indexes])
		globals()[thisdfname + df_output_name].to_parquet(thisdfname + df_output_name+".parquet", engine="pyarrow", index=False)
		upload_files.append(UploadedFile(name=thisdfname + df_output_name, local_path=thisdfname + df_output_name+".parquet", format=LocalFormat.PARQUET_TABLE))
OmicsExpressionTPMLogp1_primary = globals()['OmicsExpressionTPMLogp1HumanProteinCodingGenes'+stranded_suffix].loc[globals()['OmicsExpressionTPMLogp1HumanProteinCodingGenes'+stranded_suffix]["is_default_entry"]]
OmicsExpressionTPMLogp1_primary = OmicsExpressionTPMLogp1_primary.drop(['profile_id', 'is_default_entry'], axis=1)
OmicsExpressionTPMLogp1_primary.to_parquet("OmicsExpressionTPMLogp1_primary.parquet", engine="pyarrow", index=False)
upload_files.append(UploadedFile(name="OmicsExpressionTPMLogp1"+stranded_suffix, local_path="OmicsExpressionTPMLogp1_primary.parquet", format=LocalFormat.PARQUET_TABLE))

bigfusiontable.rename(columns={"#gene1":"gene1(HGNC ID)", "gene2":"gene2(HGNC ID)"}, inplace=True)
fusion_output_columns = ['profile_id', 'model_id',  
						 'Canonical_FusionName', 'gene1(HGNC ID)','gene2(HGNC ID)', 
						 'total_reads_in_sample', 
						 'total_reads_supporting_fusion',  'FFPM', 
						 'confidence','split_reads1', 'split_reads2', 'discordant_mates', 
						 'strand1(gene/fusion)', 'strand2(gene/fusion)',  'reading_frame',
						 'breakpoint1', 'breakpoint2', 'site1', 'site2', 'type', 'coverage1', 'coverage2',
						 'tags', 'retained_protein_domains',
						 'direction1', 'direction2']
upload_files.append(UploadedFile(name="OmicsFusionFiltered", local_path="OmicsFusionFiltered.parquet", format=LocalFormat.PARQUET_TABLE))
bigfusiontable[fusion_output_columns].to_parquet("OmicsFusionFiltered_supplementary.parquet", engine="pyarrow", index=False) 
upload_files.append(UploadedFile(name="OmicsFusionFiltered_supplementary", local_path="OmicsFusionFiltered_supplementary.parquet", format=LocalFormat.PARQUET_TABLE))

tc = create_taiga_client_v3()
tc.update_dataset(permaname=release_date, reason="First upload of all output files from new RNA-seq pipeline", additions=upload_files)
#tc.create_dataset(name="test_all_rna", description="dryrun", files=upload_files)