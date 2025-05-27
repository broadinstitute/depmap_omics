# pip install gcsfs
# pip install fsspec
# pip install biomart
# fissfc entity_tsv -w depmap-omics-rna -p broad-firecloud-ccle -t sample | sed "s/\[\"//g" | sed 's/\"\]//g' > terra_data_table_rnaseq_25q2.tsv

from collections import defaultdict
from taigapy import create_taiga_client_v3
from taigapy.client_v3 import UploadedFile, LocalFormat
import argparse
import pandas as pd
import numpy as np
import signatureanalyzer as sa

parser = argparse.ArgumentParser()
parser.add_argument("--terra_table", type=str, required=True, help='Terra table to use')
parser.add_argument("--sample_metadata", type=str, required=True, help='Sample metadata file to use')
parser.add_argument("--release_permaname", type=str, required=True, help='Permaname to use for release')


args = parser.parse_args()
terra_table = args.terra_table # test with terra_samples_file="/localstuff/terra_table_25Q2_SIGNATURES.tsv"
sample_metadata = args.sample_metadata # test with sample_metadata_file="/localstuff/2025-05-01-master-mapping-table_v4-internal-release-date-2025-05-01-master-mapping-table.csv"
release_date = args.release_permaname # Permaname to use for release
#export GOOGLE_APPLICATION_CREDENTIALS=/localstuff/depmap-omics-9764fbdbe040.json
terra_samples = pd.read_table(terra_table)
samples_to_process_all = pd.read_csv(sample_metadata)

samples_to_process = samples_to_process_all.loc[(samples_to_process_all["datatype"] == "wgs") & (samples_to_process_all["is_default_entry"] == True)]
samples = pd.merge(terra_samples, samples_to_process, left_on ="entity:sample_id", right_on="sequencing_id", how="inner")

tc = create_taiga_client_v3()

bigsigtable = pd.DataFrame()

for sample_id, sample_data in list(samples.iterrows()):
	sample_key = sample_data["entity:sample_id"]
	print(sample_id)
	if pd.notna(sample_data["mutational_sig_row"]):
		df = pd.read_csv(sample_data["mutational_sig_row"],)
		df['model_id'] = sample_data["model_id"]
		bigsigtable = pd.concat([bigsigtable, df], axis=0)
	else:
		raise ValueError(sample_key + " does not have signature output")

bigsigtable.rename(columns={"Unnamed: 0": "sequencing_id"}, inplace=True)
#bigsigtable.drop(columns=["max","max_id","max_norm", "sequencing_id"], inplace=True)
bigsigtable.columns = [
    col.split("_")[0] if col.startswith("SBS") else col
    for col in bigsigtable.columns
]
bigsigtable = bigsigtable.rename(columns={"model_id": "ModelID"})
#bigsigtable = bigsigtable.groupby(["ModelID"]).mean()
bigsigtable.to_csv("/localstuff/MolecularSignatureMatrix.csv", index=True)
tc = create_taiga_client_v3()
uploadfiles = []
#uploadfiles.append(UploadedFile(name="MolecularSignatureMatrix",local_path="/localstuff/MolecularSignatureMatrix.csv", format=LocalFormat.CSV_MATRIX))
etiologies = sa.context.signature_composite
etiologies_df = pd.DataFrame.from_dict(etiologies, orient='index').reset_index()
etiologies_df.columns = ['Signature_ID', 'label']
etiologies_df.to_csv("/localstuff/MolecularSignatureEtiologies.csv", index=False)
#uploadfiles.append(UploadedFile(name="MolecularSignatureEtiologies", local_path="/localstuff/MolecularSignatureEtiologies.csv", format=LocalFormat.CSV_TABLE))
#tc.update_dataset(permaname=release_date, reason="25Q2 molecular signature matrix", additions=uploadfiles)


