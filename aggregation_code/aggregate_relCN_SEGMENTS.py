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

relcn_column = "cnv_segments"

model_cds_dict = dict(zip(samples["SequencingID"],samples["ModelID"]))
mc_cds_dict = dict(zip(samples["SequencingID"],samples["ModelConditionID"]))
is_default_cds_dict_mc = dict(zip(samples["SequencingID"],samples["IsDefaultEntryForMC"]))

bigsegmentstable = pd.DataFrame()
id_columns = ['SequencingID','ModelID','IsDefaultEntryForModel','ModelConditionID','IsDefaultEntryForMC']
for sample_id, sample_data in list(samples.iterrows()):
	sample_key = sample_data["entity:sample_id"]
	model_id = sample_data["ModelID"]
	model_condition_id = sample_data["ModelConditionID"]
	is_default_entry_model = sample_data["IsDefaultEntryForModel"]
	is_default_entry_mc = sample_data["IsDefaultEntryForMC"]
	print(sample_id)
	if pd.notna(sample_data[relcn_column]):
		df = pd.read_table(sample_data[relcn_column])
		df['SequencingID'] = sample_key
		df['ModelID'] = model_id
		df['IsDefaultEntryForModel'] = is_default_entry_model
		df['ModelConditionID'] = model_condition_id
		df['IsDefaultEntryForMC'] = is_default_entry_mc
		bigsegmentstable = pd.concat([bigsegmentstable, df], ignore_index=True)
	else:
		raise ValueError(sample_key + " does not have CN output")

upload_files = []
bigsegmentstable['SEGMENT_COPY_NUMBER'] = np.exp2(bigsegmentstable['LOG2_COPY_RATIO_POSTERIOR_50'])
bigsegmentstable = bigsegmentstable.drop(columns=['LOG2_COPY_RATIO_POSTERIOR_50']) 
bigsegmentstable = bigsegmentstable[id_columns + [col for col in bigsegmentstable.columns if col not in id_columns]]
bigsegmentstable.to_parquet("OmicsCNSegments_MC_WGS.parquet",engine="pyarrow",  index=True)

upload_files.append(UploadedFile(name="OmicsCNSegments_MC_WGS", local_path="OmicsCNSegments_MC_WGS.parquet", format=LocalFormat.PARQUET_TABLE))
tc = create_taiga_client_v3()
tc.update_dataset(permaname=release_date, reason="CN aggregated table", additions=upload_files)
