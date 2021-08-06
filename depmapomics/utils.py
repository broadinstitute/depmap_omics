import io
import os.path
import pandas as pd
from biomart import BiomartServer
from genepy.utils.helper import createFoldersFor

from depmapomics.config import CACHE_PATH, ENSEMBL_SERVER_V

defattr = ['ensembl_gene_id', 'clone_based_ensembl_gene', 'hgnc_symbol', 'gene_biotype',
  'entrezgene_id']

def _fetchFromServer(ensemble_server, attributes):
  server = BiomartServer(ensemble_server)
  ensmbl = server.datasets['hsapiens_gene_ensembl']
  res = pd.read_csv(io.StringIO(ensmbl.search({
    'attributes': defattr+attributes
  }, header=1).content.decode()), sep='\t')
  return res

def generateGeneNames(ensemble_server=ENSEMBL_SERVER_V, useCache=False, cache_folder=CACHE_PATH,
  attributes=[]):
  """generate a genelist dataframe from ensembl's biomart

  Args:
      ensemble_server ([type], optional): [description]. Defaults to ENSEMBL_SERVER_V.
      useCache (bool, optional): [description]. Defaults to False.
      cache_folder ([type], optional): [description]. Defaults to CACHE_PATH.

  Raises:
      ValueError: [description]

  Returns:
      [type]: [description]
  """
  assert cache_folder[-1] == '/'
  cache_folder = os.path.expanduser(cache_folder)
  createFoldersFor(cache_folder)
  assert(os.path.exists(cache_folder))
  cachefile = os.path.join(cache_folder, 'biomart_ensembltohgnc.csv')
  if useCache & os.path.isfile(cachefile):
    print('fetching gene names from biomart cache')
    res = pd.read_csv(cachefile)
  else:
    print('downloading gene names from biomart')
    res = _fetchFromServer(ensemble_server, attributes)
    res.to_csv(cachefile, index=False)

  res.columns = defattr+attributes
  if type(res) is not type(pd.DataFrame()):
    raise ValueError('should be a dataframe')
  res = res[~(res[
    'clone_based_ensembl_gene'].isna() & res['hgnc_symbol'].isna())]
  res.loc[res[res.hgnc_symbol.isna()].index, "hgnc_symbol"] = \
    res[res.hgnc_symbol.isna()
                  ]['clone_based_ensembl_gene']

  return res
