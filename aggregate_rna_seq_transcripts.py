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


ens_gene_biotype = dict(zip(hgnc_table["ensembl_gene_id"], hgnc_table["locus_group"]))

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
geneorder = pd.read_table(samples.loc[0,quant_trancripts_column]).sort_values(by="Name")["Name"].reset_index(drop=True)
human_genes_all_indexes = geneorder[geneorder.str.startswith("ENST")].index
virus_genes_all_indexes = geneorder[~geneorder.str.startswith("ENST")].index

all_tpms_list = []
all_counts_list = []
all_lengths_list = []
all_raw_counts_list = []
sample_labels = []
bigfusiontable = pd.DataFrame()

for sample_id, sample_data in list(samples.iterrows())[:2]:
	sample_key = sample_data["entity:sample_id"]
	print(sample_id)
	if pd.notna(sample_data[quant_trancripts_column]):
		df = pd.read_table(sample_data[quant_trancripts_column]).sort_values(by="Name").reset_index(drop=True)
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

# Convert lists to DataFrame at the end
df_all_tpms = pd.DataFrame(np.column_stack(all_tpms_list), columns=sample_labels)
df_all_tpms.insert(0, 'Name', geneorder)
df_all_counts = pd.DataFrame(np.column_stack(all_counts_list), columns=sample_labels)
df_all_counts.insert(0, 'Name', geneorder)
df_all_raw_counts = pd.DataFrame(np.column_stack(all_raw_counts_list), columns=sample_labels)
df_all_raw_counts.insert(0, 'Name', geneorder)
df_all_lengths = pd.DataFrame(np.column_stack(all_lengths_list), columns=sample_labels)
df_all_lengths.insert(0, 'Name', geneorder)

df_all_tpms_cp = df_all_tpms.copy()
df_all_counts_cp = df_all_counts.copy()
df_all_lengths_cp = df_all_lengths.copy()
df_all_raw_counts_cp = df_all_raw_counts.copy()

upload_files = []

df_dict = {
	"OmicsExpressionTranscriptTPMLogp1": df_all_tpms_cp,
	"OmicsExpressionTranscriptExpectedCount": df_all_counts_cp,
	"OmicsExpressionTranscriptEffectiveLength": df_all_lengths_cp,
}
df_outputs = {
	"HumanAllGenes"+stranded_suffix:human_genes_all_indexes,
	"HumanProteinCodingGenes"+stranded_suffix:human_genes_pc_indexes,
	"VirusAllGenes"+stranded_suffix:virus_genes_all_indexes
}
for thisdfname, thisdf in df_dict.items():
	thisdf["transcript_id"] = geneorder
	print(thisdfname)
	thisdf.set_index("transcript_id", inplace=True)
	thisdf = thisdf.drop(['Name'], axis = 1)
	if thisdfname == "HumanProteinCodingGenes"+stranded_suffix:
		thisdf = np.log2(thisdf + 1)
	thisdf = thisdf.T
	profile_id = thisdf.index.map(profile_cds_dict)
	model_id = thisdf.index.map(model_cds_dict)
	is_default_entry = thisdf.index.map(is_default_cds_dict)
	metadatadf = pd.DataFrame({"seqid":thisdf.index, "profile_id":profile_id, "is_default_entry":is_default_entry, "model_id":model_id})
	metadatadf.set_index("seqid", inplace=True)
	# Create an upload file for the human and virus genes
	for df_output_name, df_output_indexes in df_outputs.items():
		globals()[thisdfname + df_output_name] = metadatadf.join(thisdf.iloc[:, df_output_indexes])
		globals()[thisdfname + df_output_name].to_parquet(thisdfname + df_output_name+".parquet", engine="pyarrow", index=False)
		upload_files.append(UploadedFile(name=thisdfname + df_output_name, local_path=thisdfname + df_output_name+".parquet", format=LocalFormat.PARQUET_TABLE))

tc = create_taiga_client_v3()
tc.update_dataset(permaname=release_date, reason="Add transcript level output files from new RNA-seq pipeline", additions=upload_files)
#tc.create_dataset(name="test_all_rna", description="dryrun", files=upload_files)


