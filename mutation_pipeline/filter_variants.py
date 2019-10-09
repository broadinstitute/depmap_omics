import argparse
import pandas as pd
import numpy as np
import re
import os

############################################
# Filter variants based on several databases
# Inputs:
#   See args
# Outputs:
# mcadosch 8/17
############################################

#######################################################################################################
# 0. Parse command line inputs
#######################################################################################################
parser = argparse.ArgumentParser(description='Filter variant SNVs called')
parser.add_argument('--sample_type', type=str, required=True,
                    help='Sample type: Normal or Tumor')
parser.add_argument('--sample_id', type=str, required=True,
                    help='Sample ID')
parser.add_argument('--mutect1_oncotated_maf', type=str, required=True,
                    help='Path to Oncotator MAF output from MuTect1 variant calls')
parser.add_argument('--mutect1_callstats', type=str, required=True,
					help='Path to callstats output from MuTect1')
parser.add_argument('--mutect2_oncotated_maf', type=str, required=True,
                    help='Path to Oncotator MAF output from MuTect2 variant calls')
parser.add_argument('--TCGAhotspots_path', type=str, required=True,
                    help='Path to TCGA hotspot mutations table')
parser.add_argument('--match_normal_sample_id', type=str, required=False,
					help='Sample ID of match normal')
parser.add_argument('--match_normal_filtered_variants', type=str, required=False,
					help='Path to table of filtered variants of match normal (only if this sample type is Tumor and has match normal)')

####################################################################################################
# 1. Read data into DF
####################################################################################################
def read_data():
	print("Reading data...")
	# SNV data comes from MuTect1 and Indels come from MuTect2

	# Columns of interest
	cols_of_interest = [  'Chromosome','Start_position', 'End_position', 'Hugo_Symbol', \
						  'Strand', 'Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2', \
						  'cDNA_Change', 'Variant_Classification', \
						  'Genome_Change', 'Protein_Change', \
						  'COSMIC_n_overlapping_mutations', 'COSMIC_overlapping_mutation_descriptions', \
						  'COSMIC_total_alterations_in_gene', 'COSMIC_overlapping_mutations', \
						  'gencode_transcript_name', 'ExAC_AF',  \
						  'Variant_Type', 'Entrez_Gene_Id', \
						  'Tumor_Sample_Barcode', 'Matched_Norm_Sample_Barcode', \
						  't_alt_count', 't_ref_count', \
						  'dbSNP_Val_Status', 'gene_type', \
						  'CGC_Tumor_Types_Germline', 'CGC_Tumor_Types_Somatic', 'CGC_Other_Diseases']



	### SNVs
	# Read MuTect1 Oncotated file
	mutect1_oncotated 		= pd.read_table(args.mutect1_oncotated_maf, comment='#', encoding='latin-1')
	# Keep columns of interest only if they are available
	cols_to_keep_mutect1 	= [c for c in cols_of_interest if c in mutect1_oncotated.columns]
	missing_cols 			= [c for c in cols_of_interest if c not in mutect1_oncotated.columns]
	mutect1_oncotated 		= mutect1_oncotated[cols_to_keep_mutect1]
	# For the columns of interest that were not created in the Oncotated MuTect1 file, create them and set their values to NaN
	print("This sample is missing the columns of interest: %s"%missing_cols)
	for missing_col in missing_cols:
		mutect1_oncotated[missing_col] = np.nan

	# Unfortunately Oncotator does not annotate tumor_f, so we must use callstats.txt data to retrieve it
	# Use callstats data to obtain tumor_f
	mutect1_callstats 		= pd.read_table(args.mutect1_callstats, comment='#', usecols=['position', 'tumor_f'])
	snvs 					= pd.merge(mutect1_oncotated, mutect1_callstats, left_on='Start_position', right_on='position')

	### Indels
	mutect2_oncotated 		= pd.read_table(args.mutect2_oncotated_maf, comment='#', encoding='latin-1')
	cols_of_interest_mutect2 = cols_of_interest + ['tumor_f']
	# Keep columns of interest only if they are available
	cols_to_keep_mutect2 	= [c for c in cols_of_interest_mutect2 if c in mutect2_oncotated.columns]
	mutect2_oncotated 		= mutect2_oncotated[cols_to_keep_mutect2]
	indels 					= mutect2_oncotated[ (mutect2_oncotated['Variant_Type'] == "INS") | (mutect2_oncotated['Variant_Type'] == "DEL") ]

	### Merge SNVs + Indels
	data  			= pd.concat([snvs, indels], axis=0)
	return data

####################################################################################################
# 2. Filter out non-somatic variant calls (germline)
####################################################################################################
def filter_germline_mutations(data):
	print("Filtering Germline Mutations...")
	print("There are %s unfiltered variant calls" % data.shape[0])

	### Filter germline calls found in ExAC
	data['to_remove_exac'] = False
	data.loc[ pd.notnull(data['ExAC_AF']), 'to_remove_exac'] = True
	print("Filtered %s germline variants found in exac" % data['to_remove_exac'].sum() )

	### Remove variants with low tumor fraction
	data['to_remove_low_tf'] = False
	# Set tumor_f vaules to numeric type
	data['tumor_f'] = pd.to_numeric(data['tumor_f'], errors="coerce").fillna(0)
	data.loc[ data['tumor_f'] < 0.03, 'to_remove_low_tf'] = True
	print("Filtered %s germline variants with low tumor fraction" % data['to_remove_low_tf'].sum() )

	### Rescue somatic mutations found in COSMIC
	data['to_rescue_cosmic'] = False
	data.loc[ (data['COSMIC_n_overlapping_mutations'] > 3) & (data['tumor_f'] >= 0.03), 'to_rescue_cosmic'] = True
	print("Rescued %s somatic variants found in cosmic" % data['to_rescue_cosmic'].sum() )

	### For normal samples, remove mutations that may be germline
	data['to_remove_normal_germline'] = False
	if args.sample_type == "Normal":
	    data.loc[ (((data['tumor_f'] < 0.55) & (data['tumor_f'] > 0.4)) | (data['tumor_f'] > 0.95)) , 'to_remove_normal_germline'] = True
	print("Filtered %s germline variants in normal sample" % data['to_remove_normal_germline'].sum() )

	### Rescue TCGA hotspot somatic mutations
	tcga_hotspots = pd.read_table(args.TCGAhotspots_path)
	# TCGA hotspots found in our list of variant calls
	tcga_hotspots_in_data = pd.merge(data, tcga_hotspots, left_on="Start_position", right_on="pos")
	data['to_rescue_tcga_hotspots'] = data.apply(lambda row: row['Start_position'] in tcga_hotspots_in_data['Start_position'].tolist(), axis=1)
	print("Rescued %s somatic variants found in TCGA" % data['to_rescue_tcga_hotspots'].sum() )

	### Remove variants found in match normal
	data['to_remove_match_normal_variants'] = False
	if args.sample_type == "Tumor" and args.match_normal_sample_id != "NA":
		print("Tumor sample with match normal")
		match_normal_variants = pd.read_table(args.match_normal_filtered_variants)
		intersecting_variants = pd.merge(data, match_normal_variants, on="Start_position", how="inner", suffixes=["_tumor", "_normal"])
		data['to_remove_match_normal_variants'] = data['Start_position'].isin(intersecting_variants['Start_position'])
	print("Filtered %s germline variants in found match normal sample" % data['to_remove_match_normal_variants'].sum() )

	print("Adding number of reads")
	# Set non-integer values and NaNs in "t_alt_count" column to zero
	data['t_alt_count'] = pd.to_numeric(data['t_alt_count'], errors="coerce").fillna(0)
	data['t_ref_count'] = pd.to_numeric(data['t_ref_count'], errors="coerce").fillna(0)
	data['num_reads'] = data['t_alt_count'] + data['t_ref_count']

	print("Filtering calls with low read count")
	data['insufficient_reads'] = False
	data.loc[data['num_reads'] <= 20, 'insufficient_reads'] = True
	print("There are %s calls with insufficient number of supporting reads"%data['insufficient_reads'].sum())

	print("Filtering non-applicable variant types")
	data['nonapplicable_variant_type'] = False
	data.loc[data['Variant_Classification'].isin(["RNA", "UTR", "3'UTR" "IGR", "Intron", "Silent", "Flank"]), 'nonapplicable_variant_type'] = True
	print("There are %s calls with non-applicable variant type"%data['nonapplicable_variant_type'].sum())

	print("Filtering calls with unknown Hugo symbol")
	data['hugo_symbol_unknown'] = False
	data.loc[data['Hugo_Symbol']=="Unknown", 'hugo_symbol_unknown'] = True
	print("There are %s calls with unknown Hugo symbol"%data['hugo_symbol_unknown'].sum())

	data.loc[:, 'keep'] = 	~(data['to_remove_normal_germline'] | data['to_remove_exac'] | data['to_remove_low_tf'] | data['to_remove_match_normal_variants'] ) | \
							 (data['to_rescue_cosmic'] | data['to_rescue_tcga_hotspots'])

	data.loc[:, 'keep'] = data.loc[:, 'keep'] & ~(data['insufficient_reads']) & ~(data['nonapplicable_variant_type']) & ~(data['hugo_symbol_unknown'])

	return data

####################################################################################################
# 4. Save filtered variants
####################################################################################################
def save_data(data):
	### Save file indicating if sample has clear SNVs, if it has at least one SNV
	os.system("touch ./clear_snvs.txt")
	out_file = open("clear_snvs.txt", "w")
	### No mutations
	if data.shape[0] == 0:
		print("No SNVs found to filter")
		out_file.write("false")
		out_file.close()
	else:
		snvs_to_keep = data.loc[data['keep'] == True]
		# No mutations to keep
		if snvs_to_keep.shape[0] > 1:
			out_file.write("true")
		else:
			out_file.write("false")
		out_file.close()

	data.to_csv("%s.filtered.variants.tsv"%args.sample_id, sep="\t", index=None)

	return

####################################################################################################
# MAIN
####################################################################################################
args 					= parser.parse_args()
data 					= read_data()

if data.shape[0] == 0:
	save_data(data)
	import sys; sys.exit()

data_filtered_germline 	= filter_germline_mutations(data)
save_data(data_filtered_germline)
