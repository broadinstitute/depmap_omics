#cn.py
from IPython.display import Image, display
from matplotlib import pyplot as plt
import dalmatian as dm
import pandas as pd
from taigapy import TaigaClient
tc = TaigaClient()
import seaborn as sns
import os
from genepy import mutations as mut

def renameColumns(df):
	"""
	rename some of the main columns names from RSEM, GATK.. to more readable column names

	Args:
	-----
		df: the df to rename

	Returns:
	------
		df the renamed df
	"""
	return df.rename(columns={'Sample': 'DepMap_ID', 'CONTIG': 'Chromosome', 'START': 'Start',
														'END': 'End', 'seqnames': 'Chromosome', 'start': 'Start', 'end': 'End'})

def loadFromGATKAggregation(refworkspace,  sortby=["DepMap_ID", 
																										'Chromosome', 
																										"Start", 
																										"End"],
														todrop=[], showPlots=False, colname="combined_seg_file",
														plotColname="modeled_segments_plot_tumor", tempFolder="temp/",
														sampleset="all", colRenaming={'CONTIG': 'Chromosome',
																													'START': 'Start',
																													'END': 'End',
																													'end': 'End',
																													'seqnames': 'Chromosome', 
																													'start': 'Start',
																													'Sample': "DepMap_ID",
																													'NUM_POINTS_COPY_RATIO': 'Num_Probes',
																													'MEAN_LOG2_COPY_RATIO': 'Segment_Mean',
																													'CALL': 'Status'}):
	"""
 	"""
	wm = dm.WorkspaceManager(refworkspace)
	aggregated = wm.get_entities(
		'sample_set').loc[sampleset, colname]
	segments = pd.read_csv(aggregated, sep = '\t').rename(columns=colRenaming)
 
	# TODO: copy allelic calls as well
	segments = segments[~segments.DepMap_ID.isin(todrop)].reset_index(drop=True)
	if "chr" in segments['Chromosome'][0]:
 		segments['Chromosome'] = [i[3:] for i in segments['Chromosome']]
	segments.Segment_Mean = 2**segments.Segment_Mean
	segments.Start = segments.Start.astype(int)
	segments.End = segments.End.astype(int)
	segments.loc[segments[segments.Chromosome.isin(
			['X', 'Y'])].index, 'Segment_Mean'] = segments[segments.Chromosome.isin(['X', 'Y'])]['Segment_Mean']/2
	segments = segments.sort_values(by=sortby)
	
	print("loading "+ str(len(set(segments.DepMap_ID)))+ " rows")
	if showPlots:
		# plotting results of CN calls for this new sample set
		for i, (k, val) in enumerate(wm.get_samples().loc[wm.get_sample_sets().loc[
				sampleset].samples].iterrows()):
			plot = val[plotColname]
			os.system('gsutil cp '+plot+' '+tempFolder)
			print(k)
			print(val['arxspan_id'], val['sex'])
			if i > 30:
					continue
			display(Image(os.path.join(tempFolder,plot.split('/')[-1])))
	return segments


def postProcess(refwm, sampleset='all', sortby=[
            "DepMap_ID", 'Chromosome', "Start", "End"], topdrop=[],
        genechangethresh=0.025, segmentsthresh=2000):
	print('loading WES from Terra')
	segments = loadFromGATKAggregation(refwm, sampleset=sampleset, sortby=sortby, todrop=todrop)
	print('making gene level copy number')
 
	genecn = mut.toGeneMatrix(mut.manageGapsInSegments(segments), gene_mapping)

	# validation step
	print('summary of the gene cn data:')
	print(genecn.values.min(), genecn.values.mean(), genecn.values.max())
	mut.checkGeneChangeAccrossAll(genecn, thresh=genechangethresh)
	wesfailed = mut.checkAmountOfSegments(segments, thresh=segmentsthresh)
	print("failed our QC")
	print(wesfailed)
	
 %store wesfailed
 
	segments = segments[~segments.DepMap_ID.isin(
		set(wesfailed)-set(wes_toprefer.keys()))].reset_index(drop=True)
	genecn = genecn[~genecn.index.isin(set(wesfailed)-set(wes_toprefer.keys()))]

	#resetting the source
	for v in set(segments.DepMap_ID):
			segments.loc[segments[segments.DepMap_ID == v].index,
									'Source'] = ccle_refsamples[ccle_refsamples.index == v].source.values[0]
	segments.Source = segments.Source.replace({'CCLF': 'Broad WES', 'CHORDOMA': 'Chordoma WES', 'SANGER': 'Sanger WES',
																						'IBM': 'Broad WES', np.nan: 'Broad WES', 'DEPMAP': 'Broad WES', 'IBM WES': "Broad WES", 'Broad CCLF': "Broad WES"})

	#saving
	print('saving files')
	segments.to_csv('temp/segments_allWES_withreplicates_' +
									samplesetname+'.csv', index=False)
	genecn.to_csv('temp/gene_cn_allWES_withreplicates_'+samplesetname+".csv")

def CCLEPostProcessing(wesrefwm, sampleset="all", wrongwes=set(), deletedwes=set(),
		):
  postProcess(wesrefwm, sampleset, wrongwes | deletedwes)
