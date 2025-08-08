from taigapy import create_taiga_client_v3
from taigapy.client_v3 import UploadedFile, LocalFormat
import gumbo_rest_client 
import pandas as pd
import json
import datetime
import numpy as np
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("--date_on_or_before", type=str, required=True, help='Date for filtering output. Enter in "YYYY-MM-DD" format or "None" to include all data')
args = parser.parse_args()
release_date = args.date_on_or_before

release_date = '2025-09-01'
#sed '/_comment/d' model_to_omics_data_config.json > model_to_omics_data_config.nocomments.json
gc = gumbo_rest_client.Client(
		authed_session=gumbo_rest_client.create_authorized_session(
			use_default_service_account=True
			),
			username="depmap_omics_upload",
			base_url=gumbo_rest_client.const.prod_url
		)

omics_mapping_table = gc.get("omics_mapping")
omics_mapping_table.priority[omics_mapping_table['datatype'] == 'wes'] = 1000
omics_mapping_table['sequence_type'] = omics_mapping_table['datatype']
omics_mapping_table['sequence_type'][(omics_mapping_table.datatype == 'wgs') | (omics_mapping_table.datatype == "wes")] = 'dna'
date_fields = ['internal_release_date', 'public_release_date']
omics_profile = gc.get("omics_profile")[['id','product','shared_to_dbgap']+date_fields].set_index('id')
omics_sequencing = gc.get("omics_sequencing")[['id','sequencing_date','stranded']].set_index('id')
models = gc.get("model")[['id','stripped_cell_line_name','depmap_model_type_id','growth_pattern']].set_index('id')
model_conditions = gc.get("model_condition")[['id','source','cell_format','media_id']].set_index('id')
depmap_codes = gc.get("depmap_model_type")[['id', 'oncotree_code', 'lineage', 'subtype', 'primary_disease']].set_index('id')
omics_mapping_all_merged = omics_mapping_table.merge(omics_profile, left_on='omics_profile_id', right_on='id', how = "inner")
omics_mapping_all_merged = omics_mapping_all_merged.merge(omics_sequencing, left_on='omics_sequencing_id', right_on='id', how = "inner")
omics_mapping_all_merged = omics_mapping_all_merged.merge(models, left_on='model_id', right_on='id', how = "inner")
omics_mapping_all_merged = omics_mapping_all_merged.merge(model_conditions, left_on='model_condition_id', right_on='id', how = "inner")
omics_mapping_all_merged = omics_mapping_all_merged.merge(depmap_codes, left_on='depmap_model_type_id', right_on='id', how = "inner")
output_field_dict = {"model_id" : "ModelID",
"isDefaultEntryForModel" : "IsDefaultEntryForModel",
"stripped_cell_line_name" : "StrippedCellLineName",
"depmap_model_type_id" : "DepMapCode",
"lineage": "Lineage",
"model_condition_id": "ModelConditionID",
"isDefaultEntryForMC" : "IsDefaultEntryForMC",
"source": "SourceModelCondition",
"cell_format": "CellFormat",
"media_id": "GrowthMedia",
"growth_pattern": "GrowthPattern",
"omics_profile_id" : "ProfileID",
"product": "SequencingPlatform",
"omics_sequencing_id": "SequencingID",
"datatype": "DataType",
"stranded": "Stranded",
"shared_to_dbgap": "SharedToDbGaP",
"sequencing_date": "SequencingDate"
}
for field in date_fields:
	omics_mapping_all_merged[field].fillna("2262-04-11", inplace=True)  # Use a far future date to avoid filtering out data
	omics_mapping_all_merged[field] = pd.to_datetime(omics_mapping_all_merged[field]) # Convert to datetime
	omics_mapping = omics_mapping_all_merged.loc[omics_mapping_all_merged[field] <= pd.to_datetime(release_date)] # Filter by date <= release_date
	omics_default_index_for_model_id = omics_mapping.groupby(['model_id','sequence_type'])['priority'].idxmin() # Get index of highest priority entry for each model
	omics_mapping_by_model_and_mc = omics_mapping.copy()
	omics_mapping_by_model_and_mc['isDefaultEntryForModel'] =  'No'
	omics_mapping_by_model_and_mc.loc[omics_default_index_for_model_id,'isDefaultEntryForModel'] = 'Yes'
	# By model condition
	omics_default_index_for_model_condition_id = omics_mapping_by_model_and_mc.groupby(['model_condition_id','sequence_type'])['priority'].idxmin() # Get index of highest priority entry for each model condition
	omics_mapping_by_model_and_mc['isDefaultEntryForMC'] =  'No'
	omics_mapping_by_model_and_mc.loc[omics_default_index_for_model_condition_id,'isDefaultEntryForMC'] = 'Yes'
	omics_mapping_by_model_and_mc = omics_mapping_by_model_and_mc.rename(columns = output_field_dict)
	#omics_mapping_by_model_and_mc.set_index(['ModelConditionID', 'datatype', 'isDefaultEntryForMC'], inplace = True, verify_integrity=True)
	omics_mapping_by_model_and_mc[list(output_field_dict.values())].to_csv('/localstuff/' + field + '_' + str(release_date) + '.OmicsMappingByModelAndModelCondition.csv', index = False)



#tc = create_taiga_client_v3()
#tc.update_dataset(permaname="2025-05-01-master-mapping-table-28c2",reason="Updated for bug fix pertaining to public datasets", additions=upload_files)