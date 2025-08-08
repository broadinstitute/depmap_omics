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
parser.add_argument("--release_permaname", type=str, required=True, help='Permaname to use for release')


args = parser.parse_args()
terra_table = args.terra_table # test with terra_table="/localstuff/terra_data_table_rnaseq_25q2_FINAL.tsv"
sample_metadata = args.sample_metadata # test with sample_metadata="/localstuff/2025-05-01-master-mapping-table_v4-internal-release-date-2025-05-01-master-mapping-table.csv"
release_date = args.release_permaname # Permaname to use for release

#export GOOGLE_APPLICATION_CREDENTIALS=/localstuff/depmap-omics-9764fbdbe040.json
terra_samples = pd.read_table(terra_table)
samples_to_process_all = pd.read_csv(sample_metadata)

samples_to_process = samples_to_process_all.loc[(samples_to_process_all["DataType"] == "wgs")]
samples = pd.merge(terra_samples, samples_to_process, left_on ="entity:sample_id", right_on="SequencingID", how="inner")

tc = create_taiga_client_v3()

hgnc_table = tc.get("hgnc-gene-table-e250.4/hgnc_complete_set")
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


relcn_genes = {}
relcn_column = "cnv_cn_by_gene_weighted_mean"

model_cds_dict = dict(zip(samples["SequencingID"],samples["ModelID"]))
mc_cds_dict = dict(zip(samples["SequencingID"],samples["ModelConditionID"]))
is_default_cds_dict_mc = dict(zip(samples["SequencingID"],samples["IsDefaultEntryForMC"]))
is_default_cds_dict_model = dict(zip(samples["SequencingID"],samples["IsDefaultEntryForModel"]))

df_all_relcns = pd.DataFrame()
id_columns = ['SequencingID','ModelID','IsDefaultEntryForModel','ModelConditionID','IsDefaultEntryForMC']
# Establish the order of genes in output files (by alphabetical order) by anchoring on the first sample
geneorder = pd.read_table(samples.loc[0,relcn_column],header=None, sep="\t", names=["chr","start", "end", "ensid", "genelength", "cnfull", "nbins", "relcn_gene"])
geneorder = geneorder["ensid"].sort_values().reset_index(drop=True)
ens_no_suffix = geneorder.str.split("@", n=1).str[0].str.split(".",  n=1).str[0]
ens_no_suffix.loc[geneorder.str.contains("_PAR_Y")] = ens_no_suffix.loc[geneorder.str.contains("_PAR_Y")] + "_PAR_Y"
gene_biotype = pd.Series(ens_no_suffix.map(ens_gene_biotype))
hgnc_names = pd.Series(ens_no_suffix.map(ens_to_hugo_entrez))
hgnc_names = hgnc_names.fillna(ens_no_suffix)


all_relcn_list = []
sample_labels = []

for sample_id, sample_data in list(samples.iterrows()):
	sample_key = sample_data["entity:sample_id"]
	print(sample_id)
	if pd.notna(sample_data[relcn_column]):
		df = pd.read_table(sample_data[relcn_column],header=None, sep="\t", names=["chr","start", "end", "ensid", "genelength", "cnfull", "nbins", "relcn_gene"]).sort_values(by="ensid").reset_index()
		relcn = np.exp2(df["relcn_gene"].values)  # Use NumPy arrays for speed, convert to linear scale
		allensids = df["ensid"].values
		if (allensids != geneorder).any():
			raise ValueError("Gene order is not the same for " + sample_key)
		sample_labels.append(sample_key)
		# Append data to lists
		all_relcn_list.append(relcn)
	else:
		raise ValueError(sample_key + " does not have CN output")

# Convert lists to DataFrame at the end
df_all_relcns = pd.DataFrame(np.column_stack(all_relcn_list), columns=sample_labels)
df_all_relcns.insert(0, 'GeneName', hgnc_names)

df_all_relcns_cp = df_all_relcns.copy()
upload_files = []


df_all_relcns_cp["hgnc_name"] = hgnc_names
df_all_relcns_cp.set_index("hgnc_name", inplace=True)
df_all_relcns_cp = df_all_relcns_cp.drop(['GeneName'], axis = 1)
df_all_relcns_cp = df_all_relcns_cp.T
SequencingID = df_all_relcns_cp.index.to_series()
ModelConditionID = df_all_relcns_cp.index.map(mc_cds_dict)
ModelID = df_all_relcns_cp.index.map(model_cds_dict)
isDefaultEntryForMC = df_all_relcns_cp.index.map(is_default_cds_dict_mc)
isDefaultEntryForModel = df_all_relcns_cp.index.map(is_default_cds_dict_model)
df_all_relcns_cp.loc[:,'SequencingID'] = SequencingID
df_all_relcns_cp.loc[:,'ModelID'] = ModelID
df_all_relcns_cp.loc[:,'IsDefaultEntryForModel'] = isDefaultEntryForModel
df_all_relcns_cp.loc[:,'ModelConditionID'] = ModelConditionID
df_all_relcns_cp.loc[:,'IsDefaultEntryForMC'] = isDefaultEntryForMC
df_all_relcns_cp = df_all_relcns_cp[id_columns + [col for col in df_all_relcns_cp.columns if col not in id_columns]]
df_all_relcns_cp.to_parquet("OmicsCNGeneWGS.parquet",engine="pyarrow",  index=True)

upload_files.append(UploadedFile(name="OmicsCNGeneMCWGS", local_path="OmicsCNGeneWGS.parquet", format=LocalFormat.PARQUET_TABLE))

tc = create_taiga_client_v3()
tc.update_dataset(permaname=release_date, reason="CN aggregated table", additions=upload_files)
