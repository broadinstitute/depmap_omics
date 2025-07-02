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
terra_table = args.terra_table # test with terra_table="/localstuff/terra_data_table_rnaseq_25q2_FINAL.tsv"
sample_metadata = args.sample_metadata # test with sample_metadata="/localstuff/2025-05-01-master-mapping-table_v4-internal-release-date-2025-05-01-master-mapping-table.csv"
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

samples_to_process = samples_to_process_all.loc[(samples_to_process_all["sequence_type"] == "rna")]
samples = pd.merge(terra_samples, samples_to_process, left_on ="entity:sample_id", right_on="omics_sequencing_id", how="inner")

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
model_cds_dict = dict(zip(samples["omics_sequencing_id"],samples["ModelID"]))
mc_cds_dict = dict(zip(samples["omics_sequencing_id"],samples["ModelConditionID"]))
is_default_cds_dict_mc = dict(zip(samples["omics_sequencing_id"],samples["isDefaultEntryForMC"]))
is_default_cds_dict_model = dict(zip(samples["omics_sequencing_id"],samples["isDefaultEntryForModel"]))

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

for sample_id, sample_data in list(samples.iterrows()):
	sample_key = sample_data["entity:sample_id"]
	print(sample_id)
	star_read_counts_df = pd.read_table(sample_data["reads_per_gene"], sep='\t',names=["colname","unstranded_counts","forward_stranded_counts","reverse_stranded_counts"])
	all_raw_counts_list.append(star_read_counts_df["unstranded_counts"][4:])
	total_reads_for_sample = pd.Series.sum(star_read_counts_df["unstranded_counts"][4:])
	if pd.notna(sample_data[quant_genes_column]):
		df = pd.read_table(sample_data[quant_genes_column]).sort_values(by="Name").reset_index(drop=True)
		tpms = df["TPM"].values  # Use NumPy arrays for speed
		counts = df["NumReads"].round().astype('int').values
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
all_tables = {}
for thisdfname, thisdf in df_dict.items():
	thisdf["hgnc_name"] = hgnc_names
	thisdf["hgnc_name"] = thisdf["hgnc_name"].fillna(thisdf["Name"])
	print(thisdfname)
	thisdf.set_index("hgnc_name", inplace=True)
	thisdf = thisdf.drop(['Name'], axis = 1)
	if thisdfname == "OmicsExpressionTPMLogp1":
		thisdf = np.log2(thisdf + 1)
	thisdf = thisdf.T
	ModelConditionID = thisdf.index.map(mc_cds_dict)
	ModelID = thisdf.index.map(model_cds_dict)
	isDefaultEntryMC = thisdf.index.map(is_default_cds_dict_mc)
	isDefaultEntryModel = thisdf.index.map(is_default_cds_dict_model)
	#metadatadf.set_index("seqid", inplace=True)
	# Create upload files for the human and virus genes, for counts, TPMs, effective lengths and raw counts
	for df_output_name, df_output_indexes in df_outputs.items():
		all_tables[thisdfname + df_output_name] = thisdf.iloc[:, df_output_indexes].copy()
		all_tables[thisdfname + df_output_name].loc[:,'ModelConditionID'] = ModelConditionID
		all_tables[thisdfname + df_output_name].loc[:,'isDefaultEntryForMC'] = isDefaultEntryMC
		all_tables[thisdfname + df_output_name].set_index(["ModelConditionID","isDefaultEntryForMC"], inplace=True, verify_integrity=True)
		all_tables[thisdfname + df_output_name].to_parquet(thisdfname + df_output_name+".parquet", engine="pyarrow", index=True)
		upload_files.append(UploadedFile(name=thisdfname + df_output_name, local_path=thisdfname + df_output_name+".parquet", format=LocalFormat.PARQUET_TABLE))
# This is for model level data. Select only the default entry for each model (isDefaultEntry =- True)
OmicsExpressionTPMLogp1_primary = all_tables['OmicsExpressionTPMLogp1HumanProteinCodingGenes'+stranded_suffix]
OmicsExpressionTPMLogp1_primary.loc[:,'ModelID'] = ModelID
OmicsExpressionTPMLogp1_primary.loc[:,'isDefaultEntryForModel'] = isDefaultEntryModel
OmicsExpressionTPMLogp1_primary = OmicsExpressionTPMLogp1_primary.loc[OmicsExpressionTPMLogp1_primary['isDefaultEntryForModel'] == "Yes"]
OmicsExpressionTPMLogp1_primary.set_index(["ModelID","isDefaultEntryForModel"], inplace=True, verify_integrity=True)
OmicsExpressionTPMLogp1_primary.to_parquet("OmicsExpressionTPMLogp1_primary.parquet", index=True)
upload_files.append(UploadedFile(name="OmicsExpressionProteinCodingGenesTPMLogp1"+stranded_suffix, local_path="OmicsExpressionTPMLogp1_primary.parquet", format=LocalFormat.PARQUET_TABLE))
OmicsExpressionTPMLogp1_virus = all_tables['OmicsExpressionTPMLogp1VirusAllGenes'+stranded_suffix]
OmicsExpressionTPMLogp1_virus.loc[:,'ModelID'] = ModelID
OmicsExpressionTPMLogp1_virus.loc[:,'isDefaultEntryForModel'] = isDefaultEntryModel
OmicsExpressionTPMLogp1_virus = OmicsExpressionTPMLogp1_virus.loc[OmicsExpressionTPMLogp1_virus['isDefaultEntryForModel'] == "Yes"]
OmicsExpressionTPMLogp1_virus.set_index(["ModelID","isDefaultEntryForModel"], inplace=True, verify_integrity=True)
OmicsExpressionTPMLogp1_virus.to_parquet("OmicsExpressionTPMLogp1Virus.parquet", index=True)
upload_files.append(UploadedFile(name="OmicsExpressionTPMLogp1Virus"+stranded_suffix, local_path="OmicsExpressionTPMLogp1Virus.parquet", format=LocalFormat.PARQUET_TABLE))


tc = create_taiga_client_v3()
tc.update_dataset(permaname="test-all-rna-c62c", reason="Test upload of all output files from new RNA-seq pipeline", additions=upload_files)
