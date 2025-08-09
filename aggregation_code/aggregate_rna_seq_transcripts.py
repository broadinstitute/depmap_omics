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
	quant_genes_column = "quant_genes"
	quant_trancripts_column = "quant_transcripts"
else:
	stranded_suffix=""
	quant_genes_column = "quant_genes_iu"
	quant_trancripts_column = "quant_transcripts_iu"
#export GOOGLE_APPLICATION_CREDENTIALS=/localstuff/depmap-omics-9764fbdbe040.json
terra_samples = pd.read_table(terra_table)
terra_samples['quant_transcripts_auto'].fillna(terra_samples['quant_transcripts_iu'], inplace=True)
samples_to_process_all = pd.read_csv(sample_metadata)

samples_to_process = samples_to_process_all.loc[(samples_to_process_all["DataType"] == "rna")]
samples = pd.merge(terra_samples, samples_to_process, left_on ="entity:sample_id", right_on="SequencingID", how="inner")

mc_cds_dict = dict(zip(samples["SequencingID"],samples["ModelConditionID"]))
model_cds_dict = dict(zip(samples["SequencingID"],samples["ModelID"]))
is_default_cds_dict_mc = dict(zip(samples["SequencingID"],samples["IsDefaultEntryForMC"]))
is_default_cds_dict_model = dict(zip(samples["SequencingID"],samples["IsDefaultEntryForModel"]))

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
sample_labels = []

for sample_id, sample_data in list(samples.iterrows()):
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
df_all_lengths = pd.DataFrame(np.column_stack(all_lengths_list), columns=sample_labels)
df_all_lengths.insert(0, 'Name', geneorder)

df_all_tpms_cp = df_all_tpms.copy()
df_all_counts_cp = df_all_counts.copy()
df_all_lengths_cp = df_all_lengths.copy()

upload_files = []

df_dict = {
	"OmicsExpressionTranscriptTPMLogp1": df_all_tpms_cp,
	"OmicsExpressionTranscriptExpectedCount": df_all_counts_cp,
	"OmicsExpressionTranscriptEffectiveLength": df_all_lengths_cp,
}
df_outputs = {
	"HumanAllGenes"+stranded_suffix:human_genes_all_indexes,
	"VirusAllGenes"+stranded_suffix:virus_genes_all_indexes
}
all_tables = {}

for thisdfname, thisdf in df_dict.items():
	print(thisdfname)
	thisdf = thisdf.set_index('Name')
	if thisdfname == "OmicsExpressionTranscriptTPMLogp1"+stranded_suffix:
		thisdf = np.log2(thisdf + 1)
	thisdf = thisdf.T
	SequencingID = thisdf.index.to_series()
	ModelConditionID = thisdf.index.map(mc_cds_dict)
	ModelID = thisdf.index.map(model_cds_dict)
	isDefaultEntryMC = thisdf.index.map(is_default_cds_dict_mc)
	isDefaultEntryModel = thisdf.index.map(is_default_cds_dict_model)
	id_columns = ['SequencingID','ModelID','IsDefaultEntryForModel','ModelConditionID','IsDefaultEntryForMC']
	# Create an upload file for the human and virus genes
	for df_output_name, df_output_indexes in df_outputs.items():
		all_tables[thisdfname + df_output_name] = thisdf.iloc[:, df_output_indexes].copy()
		all_tables[thisdfname + df_output_name].loc[:,'SequencingID'] = SequencingID
		all_tables[thisdfname + df_output_name].loc[:,'ModelID'] = ModelID
		all_tables[thisdfname + df_output_name].loc[:,'IsDefaultEntryForModel'] = isDefaultEntryModel
		all_tables[thisdfname + df_output_name].loc[:,'ModelConditionID'] = ModelConditionID
		all_tables[thisdfname + df_output_name].loc[:,'IsDefaultEntryForMC'] = isDefaultEntryMC
		all_tables[thisdfname + df_output_name] = all_tables[thisdfname + df_output_name][id_columns + [col for col in all_tables[thisdfname + df_output_name].columns if col not in id_columns]]
		all_tables[thisdfname + df_output_name].to_parquet(thisdfname + df_output_name+".parquet", engine="pyarrow", index=False)
		upload_files.append(UploadedFile(name=thisdfname + df_output_name, local_path=thisdfname + df_output_name+".parquet", format=LocalFormat.PARQUET_TABLE))

tc1 = create_taiga_client_v3()
tc1.update_dataset(permaname=release_date, reason="Add transcript level output files from RNA-seq pipeline", additions=upload_files)
#tc.create_dataset(name="test_all_rna", description="dryrun", files=upload_files)


