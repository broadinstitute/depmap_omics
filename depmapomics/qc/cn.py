from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
from taigapy import TaigaClient
tc = TaigaClient()


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
