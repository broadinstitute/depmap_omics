# pip install gcsfs
# pip install fsspec
# pip install biomart
# fissfc entity_tsv -w depmap-omics-24q4 -p broad-firecloud-ccle -t sample | sed "s/\[\"//g" | sed 's/\"\]//g' > terra_data_table_rnaseq.tsv

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
parser.add_argument("--rules_json", type=str, required=True, help='json config that holds rules for assigning correct dataset to model')
args = parser.parse_args()
release_date = args.date_on_or_before
config_json = args.rules_json

#sed '/_comment/d' model_to_omics_data_config.json > model_to_omics_data_config.nocomments.json
gc = gumbo_rest_client.Client(
		authed_session=gumbo_rest_client.create_authorized_session(
			use_default_service_account=True
			),
			username="depmap_omics_upload",
			base_url=gumbo_rest_client.const.prod_url
		)

configdata = json.load(open(config_json))
print(configdata)
tables_to_filter_names = list(configdata["table_filters"].keys())
filtered_tables_names = [key + "_filtered" for key in tables_to_filter_names]
for tablename in tables_to_filter_names:
	print(tablename)
	globals()[tablename]  = gc.get(tablename)
	thistableglobal = globals()[tablename]
	# Fill in missing values with a default value based on the data type of the column
	# This is to ensure that the filters work correctly
	for db_column in thistableglobal.columns:
		if pd.api.types.is_bool_dtype(thistableglobal[db_column]):
			thistableglobal[db_column].fillna(False, inplace = True)
			thistableglobal[db_column] = thistableglobal[db_column].astype(bool)
		elif pd.api.types.is_integer_dtype(thistableglobal[db_column]):
			thistableglobal[db_column].fillna(1000, inplace = True)
		elif pd.api.types.is_string_dtype(thistableglobal[db_column]):
			thistableglobal[db_column].fillna("None", inplace = True)
		elif pd.api.types.is_datetime64_any_dtype(thistableglobal[db_column]):
			thistableglobal[db_column].fillna("2262-04-11", inplace = True)
			thistableglobal[db_column] = pd.to_datetime(thistableglobal[db_column])
		else:
			thistableglobal[db_column].fillna("NA", inplace = True)

# The following code block is to filter the tables based on the exclusion and inclusion filters. The filters 
# are specified in the config file. The filtered tables are stored in a new variable with the same name as the
# original table with "_filtered" appended to the name. See config file for more details on how to specify filters
# and filter modes.

for thistablename,filteredtablename in zip(tables_to_filter_names, filtered_tables_names):
	thistable = globals()[thistablename]
	thistable_filters = configdata.get("table_filters").get(thistablename)
	exclusion_filters = thistable_filters.get("exclusion_filters", {})
	exclusion_columns = exclusion_filters.get("columns")
	exclusion_mode = exclusion_filters.get("filter_mode")
	inclusion_filters = thistable_filters.get("inclusion_filters", {})
	inclusion_columns = inclusion_filters.get("columns", [])
	inclusion_mode = inclusion_filters.get("filter_mode")
	# This code block builds the filter string to query on, based on the filter mode and the filters specified in the config file
	filter_string_ex = ''
	filter_string_in = ''
	if exclusion_mode == "all" and exclusion_filters:
		filter_string_ex = " & ".join([f"{col} not in {value!r}" for col, value in exclusion_filters.items() if col in exclusion_columns])
	elif exclusion_mode == "any" and exclusion_filters:
		filter_string_ex = " | ".join([f"{col} not in {value!r}" for col, value in exclusion_filters.items() if col in exclusion_columns])
	if inclusion_mode == "all" and inclusion_filters:
		filter_string_in = " & ".join([f"{col} in {value!r}" for col, value in inclusion_filters.items() if col in inclusion_columns])
	elif inclusion_mode == "any"  and inclusion_filters :
		filter_string_in = " | ".join([f"{col} in {value!r}" for col, value in inclusion_filters.items() if col in inclusion_columns])
	print(filter_string_ex)
	print(filter_string_in)
	print(filteredtablename)
	if filter_string_ex != '' or filter_string_in != '':
		if filter_string_in != '' and filter_string_ex != '':
			globals()[filteredtablename] = thistable.query(filter_string_ex).query(filter_string_in)
		elif filter_string_in != '' and filter_string_ex == '':
			globals()[filteredtablename] = thistable.query(filter_string_in)
		elif filter_string_ex != '' and filter_string_in == '':
			globals()[filteredtablename] = thistable.query(filter_string_ex)
		else:
			print("should never come here")



# Merge tables based on merge parameters specified in the config file
merged_table_all = None
for thistablename, thisfilteredtablename in zip(tables_to_filter_names, filtered_tables_names):
	merge_params = configdata["table_filters"][thistablename].get("merge_parameters", None)
	if merged_table_all is None:
		merged_table_all = globals()[thisfilteredtablename]
	if merge_params:
		primary_key = merge_params["primary"][0]
		foreign_key = merge_params["foreign"][0]
		merged_table_all = pd.merge(globals()[thisfilteredtablename], 
						  merged_table_all, 
						  left_on = primary_key, 
						  right_on = foreign_key, 
						  how = 'inner', 
						  suffixes=['_' + thistablename,'_merged'])

merged_all = dict()

if release_date != {'None'}:
	release_date_columns = configdata["release_date_columns"]["columns"]
	for release_date_column in release_date_columns:
		merged_table_all.loc[merged_table_all[release_date_column] == 'None', [release_date_column]] = "2262-04-11"
		merged_table_all[release_date_column] = pd.to_datetime(merged_table_all[release_date_column])
		#Separate data tables by release date
		filtered = merged_table_all.loc[merged_table_all[release_date_column] <= pd.to_datetime(release_date)]
		merged_all[release_date_column] = filtered


#merged_table = merged_table.loc[merged_table["include"] == True]

# Add priority columns to the merged table based on the priority rankings specified in the config file
# The priority columns are added to the merged table with the suffix "_priority"
# The priority ladder is used for a final outcome of selecing the best row for each model, building from the bottom up.

priority_rankings = configdata["priority_rankings"]
for release_date_column, merged_table in merged_all.items():
	merged_table = merged_all[release_date_column].copy()
	for rankingvar, ranking in priority_rankings["rankings"].items():
		if pd.api.types.is_bool_dtype(merged_table[rankingvar]):
			merged_table[rankingvar] = merged_table[rankingvar].fillna(False)
		elif pd.api.types.is_integer_dtype(merged_table[rankingvar]):
			merged_table[rankingvar] = merged_table[rankingvar].fillna(1000)
		elif pd.api.types.is_string_dtype(merged_table[rankingvar]):
			merged_table[rankingvar] = merged_table[rankingvar].fillna("None")
		else:
			merged_table[rankingvar] = merged_table[rankingvar].fillna("NA")
		rankingvalue = merged_table[rankingvar].astype(str).map(ranking)
		rankingvalue = rankingvalue.fillna(len(ranking)+1)
		merged_table[rankingvar + '_priority'] = rankingvalue
		merged_all[release_date_column] = merged_table


# Special case to accommodate for the way gumbo is structured. 
# You can have the same profile generating two sequences of the same type wgs or rna. This is something GP just delivers without asking
# You can also have two different profiles for each CDS - this would have happened when Sam manually gave each CDS ID its own profile_id 
# for clarity. At the moment it is unclear how many fall under each case, so filtering twice with same parameters (version)
# Wish I could have found a way to not hardcode column names these next 3 lines, but I am not sure how to go about that.

selcols = configdata["final_output_columns"]["columns"]
upload_files = list()
for release_date_column, merged_table in merged_all.items():
	merged_table = merged_all[release_date_column].copy()
	idx = []
	idx.append(merged_table.groupby(["sequencing_id"])["version_merged"].idxmax()) # indexes with latest version in seq_table
	idx.append(merged_table.loc[idx[0]].groupby(["profile_id"])["version_merged"].idxmax()) # indexes with latest version in seq_table
	# Profile_ID:Sequencing_ID 1:1 mapping at this stage, if we filter indexes as merged_table.loc[idx[1]]
	all_rankings = list(configdata["priority_rankings"]["rankings"].keys())
	all_rankings_priority_columns = [col + '_priority' for col in all_rankings]
	idxcount = len(idx)
	ops = (configdata["priority_rankings"]["order_of_operations"])
	for operation in ops:
		ranking_order = configdata["priority_rankings"]["ranking_orders"][operation] 
		groupbyvar = ranking_order["groupbyvar"]
		for priorityvar, minmax in zip(ranking_order["priorityvar"], ranking_order["select_index"]):
			zeropriorityrows_index = merged_table.loc[idx[idxcount-1]].loc[merged_table[priorityvar + '_priority'] == 0].index
			nonzeropriorityrows_index = merged_table.loc[idx[idxcount-1]].loc[merged_table[priorityvar + '_priority'] != 0].index
			nonzeropriorityrows = merged_table.loc[nonzeropriorityrows_index]
			grouped = nonzeropriorityrows.groupby(groupbyvar)[priorityvar + '_priority']
			if len(zeropriorityrows_index) > 0:
				all_indexes = np.concatenate([zeropriorityrows_index.tolist(),grouped.apply(minmax).tolist()])
			else:
				all_indexes = grouped.apply(minmax).tolist()
			idx.append(all_indexes)
			idxcount = idxcount + 1
		#merged_table.loc[idx[idxcount-1]].to_csv("" + operation + ".csv")
		# Pick one row per model based on the priority rankings set above. This will be default data row for that model+datatype combo
	merged_table["is_default_entry"] = 'False'
	merged_table.loc[idx[len(idx)-1],"is_default_entry"] = 'True'
	merged_all[release_date_column] = merged_table
	print(release_date_column +"***"+release_date)
	merged_all[release_date_column].loc[:,selcols].to_csv(release_date_column + "." + str(release_date) + ".master_mapping_table.csv", index=False)
	upload_files.append(UploadedFile(name="OmicsProfiles", local_path=release_date_column + "." + str(release_date) + ".master_mapping_table.csv", format=LocalFormat.CSV_TABLE))
	



tc = create_taiga_client_v3()
tc.create_dataset(name="2025-05-01-master-mapping-table-V2",description="Updated for bug fix pertaining to public datasets", files=upload_files)