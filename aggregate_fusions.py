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

profile_cds_dict = dict(zip(samples["sequencing_id"],samples["profile_id"]))
model_cds_dict = dict(zip(samples["sequencing_id"],samples["model_id"]))
is_default_cds_dict = dict(zip(samples["sequencing_id"],samples["is_default_entry"]))

sample_labels = []
bigfusiontable = pd.DataFrame()

for sample_id, sample_data in list(samples.iterrows()):
	sample_key = sample_data["entity:sample_id"]
	print(sample_id)
	star_read_counts_df = pd.read_table(sample_data["reads_per_gene"], sep='\t',names=["colname","unstranded_counts","forward_stranded_counts","reverse_stranded_counts"])
	total_reads_for_sample = pd.Series.sum(star_read_counts_df["unstranded_counts"][4:])
	if pd.notna(sample_data["fusions"]):
		fusiondf = pd.read_table(sample_data["fusions"])
		# Use STAR read output to get total reads for normalization later
		if fusiondf.shape[0] > 0:
			#fusiondf["sample_id"] = sample_key
			fusiondf["ProfileID"] = sample_data["profile_id"]
			fusiondf["ModelID"] = sample_data["model_id"]
			fusiondf["CellLine"] = sample_data["stripped_cell_line_name"]
			fusiondf["OncotreeCode"] = sample_data["depmap_code"]
			fusiondf["lineage"] = sample_data["lineage"]
			fusiondf["IsDefaultEntry"] = sample_data["is_default_entry"]
			fusiondf["TotalReadsInSample"] = total_reads_for_sample
			bigfusiontable = pd.concat([bigfusiontable, fusiondf], ignore_index=True)
		else:
			print("No fusions for " + sample_key)

# Calculate FFPM by summing (SplitReads1+SplitReads2)/total reads for sample and multiplying by 1e6
bigfusiontable["TotalReadsSupportingFusion"] = bigfusiontable["split_reads1"] + bigfusiontable["split_reads2"] + bigfusiontable["discordant_mates"]
bigfusiontable["FFPM"] = 1e6*bigfusiontable["TotalReadsSupportingFusion"]/bigfusiontable["TotalReadsInSample"]
# Canonicalize gene fusion names, making it easier to group by
bigfusiontable['gene1_clean'] = bigfusiontable['#gene1'].str.split(',').str[0].str.extract(r'^([^\(]+)')[0]
bigfusiontable['gene2_clean'] = bigfusiontable['gene2'].str.split(',').str[0].str.extract(r'^([^\(]+)')[0]
bigfusiontable['CanonicalFusionName'] = bigfusiontable.apply(lambda row: '--'.join(sorted([row['gene1_clean'], row['gene2_clean']])), axis=1)
# Special handling for CIC-DUX4 because of large numbers of identical pseudogenes, we manually set the canonical name to CIC--DUX4 when the second partner is "DUX4L* family of pseudogenes"
bigfusiontable['gene1_clean'] = bigfusiontable['gene1_clean'].replace(r'^DUX4L.*', 'DUX4', regex=True)
bigfusiontable['gene2_clean'] = bigfusiontable['gene2_clean'].replace(r'^DUX4L.*', 'DUX4', regex=True)
bigfusiontable['CanonicalFusionName'] = bigfusiontable['CanonicalFusionName'].replace(r'^CIC--DUX4L.*', 'CIC--DUX4', regex=True)
bigfusiontable['TotalFusionCoverage'] = bigfusiontable['coverage1'] + bigfusiontable['coverage2']
bigfusiontable_defaults_only = bigfusiontable.loc[(bigfusiontable['IsDefaultEntry']) & ((bigfusiontable['confidence'] == "high") | (bigfusiontable['confidence'] == "medium")) & ((bigfusiontable['site1'] != 'intergenic') | (bigfusiontable['site2'] != 'intergenic'))] 
fusion_grouped_by_sample_and_genes = bigfusiontable_defaults_only.groupby(['ModelID', 'ProfileID', 'CanonicalFusionName', 'gene1_clean', 'gene2_clean', 'strand1(gene/fusion)', 'strand2(gene/fusion)','reading_frame', 'TotalReadsInSample']).apply(
    lambda g: pd.Series({
		'SplitReads1' : g['split_reads1'].sum(),
		'SplitReads2' : g['split_reads2'].sum(),
		'DiscordantMates' : g['discordant_mates'].sum(),
		'TotalReadsSupportingFusion' : g['TotalReadsSupportingFusion'].sum(),
		'TotalFusionCoverage': g['TotalFusionCoverage'].sum()})).reset_index()

fusion_grouped_by_sample_and_genes['FFPM'] = 1e6*fusion_grouped_by_sample_and_genes["TotalReadsSupportingFusion"]/fusion_grouped_by_sample_and_genes["TotalReadsInSample"]
fusion_grouped_by_sample_and_genes.rename(columns={"gene1_clean":"Gene1", "gene2_clean": "Gene2", "strand1(gene/fusion)":"Strand1", "strand2(gene/fusion)":"Strand2", "reading_frame":"ReadingFrame"}, inplace=True)

selcols = [
'CanonicalFusionName',
'Gene1',
'Gene2',
'ProfileID',
'ModelID',
'TotalReadsSupportingFusion',
'TotalFusionCoverage',
'FFPM',
'SplitReads1',
'SplitReads2',
'DiscordantMates',
'Strand1',
'Strand2',
'ReadingFrame']
fusion_by_model_df = fusion_grouped_by_sample_and_genes[selcols]
fusion_by_model_df.to_parquet("OmicsFusionFiltered.parquet", index=True)
upload_files = []

bigfusiontable.rename(columns={"#gene1":"gene1(HGNC ID)", "gene2":"gene2(HGNC ID)"}, inplace=True)
fusion_output_columns = ['ProfileID', 'ModelID',  
						 'CanonicalFusionName', 'gene1(HGNC ID)','gene2(HGNC ID)', 
						 'TotalReadsInSample', 
						 'TotalReadsSupportingFusion',  'FFPM', 
						 'confidence','split_reads1','split_reads2', 'discordant_mates', 
						 'strand1(gene/fusion)', 'strand2(gene/fusion)',  'reading_frame',
						 'breakpoint1', 'breakpoint2', 'site1', 'site2', 'type', 'coverage1', 'coverage2',
						 'tags', 'retained_protein_domains',
						 'direction1', 'direction2']
upload_files.append(UploadedFile(name="OmicsFusionFiltered", local_path="OmicsFusionFiltered.parquet", format=LocalFormat.PARQUET_TABLE))
bigfusiontable[fusion_output_columns].to_parquet("OmicsFusionFiltered_supplementary.parquet", engine="pyarrow", index=False) 
upload_files.append(UploadedFile(name="OmicsFusionFiltered_supplementary", local_path="OmicsFusionFiltered_supplementary.parquet", format=LocalFormat.PARQUET_TABLE))

tc = create_taiga_client_v3()
tc.update_dataset(permaname=release_date, reason="Changed column names to CamelCase", additions=upload_files)

