from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
from taigapy import TaigaClient
tc = TaigaClient()

def plot_gene_cn_comparison(number_of_dots=10000, savefig=False,
                            release1={
                                'name': 'internal-20q3-00d0', 'version': 9},
                            release2={
                                'name': 'internal-21q1-4fc4', 'version': 7},
                            names=None,
                            filenames1={'gene_cn': 'CCLE_gene_cn',
                                        'segment_cn': 'CCLE_segment_cn'},
                            filenames2={'gene_cn': 'CCLE_gene_cn',
                                        'segment_cn': 'CCLE_segment_cn'}
                            ):
	if names == None:
		names = [release1['name'], release2['name']]
	CCLE_gene_cn_1 = tc.get(
		name=release1['name'], version=release1['version'], file=filenames1['gene_cn'])
	CCLE_gene_cn_2 = tc.get(
		name=release2['name'], version=release2['version'], file=filenames2['gene_cn'])
	CCLE_gene_cn_1 = CCLE_gene_cn_1.stack()
	CCLE_gene_cn_2 = CCLE_gene_cn_2.stack()
	CCLE_gene_cn_12 = pd.concat([CCLE_gene_cn_1, CCLE_gene_cn_2], axis=1)
	del CCLE_gene_cn_1, CCLE_gene_cn_2

	cols = ['{:s}.{:d}'.format(release1['name'], release1['version']),
					'{:s}.{:d}'.format(release2['name'], release2['version'])]
	CCLE_gene_cn_12.columns = cols
	CCLE_gene_cn_12.reset_index(inplace=True)
	CCLE_gene_cn_12.rename(
		columns={'level_0': 'DepMap_ID', 'level_1': 'gene'}, inplace=True)

	CCLE_segment_cn_1 = tc.get(
		name=release1['name'], version=release1['version'], file=filenames1['segment_cn'])
	CCLE_segment_cn_2 = tc.get(
		name=release2['name'], version=release2['version'], file=filenames2['segment_cn'])
	sources = pd.merge(CCLE_segment_cn_1[['DepMap_ID', 'Source']].drop_duplicates(),
											CCLE_segment_cn_2[['DepMap_ID',
																				'Source']].drop_duplicates(),
											on='DepMap_ID', suffixes=['_'+names[0], '_'+names[1]])
	del CCLE_segment_cn_1, CCLE_segment_cn_2
	sources['source_change'] = sources.apply(
		lambda x: '{:s} -> {:s}'.format(x['Source_'+names[0]], x['Source_'+names[1]]), axis=1)
	sources['source_has_changed'] = (
		sources['Source_'+names[0]] != sources['Source_'+names[1]])

	CCLE_gene_cn_12 = pd.merge(CCLE_gene_cn_12, sources, on='DepMap_ID')

	plt.figure(figsize=(20, 10))
	sns.scatterplot(data=CCLE_gene_cn_12.sample(number_of_dots, random_state=0), x=cols[0], y=cols[1],
									hue='source_change', style='source_has_changed', alpha=0.5, cmap='Tab20')

	if savefig:
		# .replace('/', '.')
		filename = ('{}-vs-{}.png'.format(cols[1], cols[0]))
		print('saving {}'.format(filename))
		plt.savefig(filename, bbox_inches='tight')

	return CCLE_gene_cn_12


def plotCNchanges(newgenecn, prevgenecn, newsegments, prevsegments, depmap_id="DepMap_ID", source="Source", prevname='prev', newname="new"):
  """
  makes a Javad Plot on the gene copy number dataset

  Args:
  -----
    newgenecn: pd dataframe the newg gene copy numbe matrix
    prevgenecn: pd dataframe the prev gene copy numbe matrix
    newsegments: pd datagrame the newg segments which contain a "source" columns
    prevsegments: pd datagrame the prev segments which contain a "source" columns
    depmap_id: str: the sample id colname
    source: str: the source colname
    prevname: str the name of newg on the plot
    newname: str the name of prev on the plot
  """
  grouped = pd.concat(
    [prevgenecn.stack(), newgenecn.stack()], axis=1)
  grouped.columns = [prevname, newname]

  grouped.reset_index(inplace=True)
  grouped.rename(
    columns={'level_0': depmap_id, 'level_1': 'gene'}, inplace=True)
  sources = pd.merge(prevsegments[[depmap_id, source]].drop_duplicates(),
                     newsegments[[depmap_id, source]].drop_duplicates(),
                     on=depmap_id, suffixes=['_'+prevname, '_'+newname])

  sources['source_change'] = sources.apply(
    lambda x: '{:s} -> {:s}'.format(x[source+'_'+prevname], x[source+'_'+newname]), axis=1)
  sources['source_has_changed'] = (
    sources[source+'_'+prevname] != sources[source+"_"+newname])
  grouped = pd.merge(grouped, sources, on=depmap_id)
  plt.figure(figsize=(20, 10))
  sns.scatterplot(data=grouped.sample(1000000, random_state=0), x=prevname, y=newname,
                  hue='source_change', style='source_has_changed', alpha=0.5, cmap='Tab20')
