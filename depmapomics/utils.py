import io
import os.path
import pandas as pd
from biomart import BiomartServer
from genepy.utils.helper import createFoldersFor

from depmapomics.config import CACHE_PATH, ENSEMBL_SERVER_V


def generateGeneNames(ensemble_server=ENSEMBL_SERVER_V, useCache=False, cache_folder=CACHE_PATH):
    """
    # TODO: to document
    """
    assert cache_folder[-1] == '/'
    createFoldersFor(cache_folder)
    cachefile = os.path.join(cache_folder, 'biomart_ensembltohgnc.csv')
    cachefile = os.path.expanduser(cachefile)
    if useCache & os.path.isfile(cachefile):
        print('fetching gene names from biomart cache')
        ensembltohgnc = pd.read_csv(cachefile)
    else:
        print('downloading gene names from biomart')
        server = BiomartServer(ensemble_server)
        ensmbl = server.datasets['hsapiens_gene_ensembl']
        ensembltohgnc = pd.read_csv(io.StringIO(ensmbl.search({
            'attributes': ['ensembl_gene_id', 'clone_based_ensembl_gene',
                           'hgnc_symbol', 'gene_biotype', 'entrezgene_id']
        }, header=1).content.decode()), sep='\t')
        ensembltohgnc.to_csv(cachefile, index=False)

    ensembltohgnc.columns = ['ensembl_gene_id', 'clone_based_ensembl_gene',
                             'hgnc_symbol', 'gene_biotype', 'entrezgene_id']
    if type(ensembltohgnc) is not type(pd.DataFrame()):
        raise ValueError('should be a dataframe')
    ensembltohgnc = ensembltohgnc[~(ensembltohgnc[
        'clone_based_ensembl_gene'].isna() & ensembltohgnc['hgnc_symbol'].isna())]
    ensembltohgnc.loc[ensembltohgnc[ensembltohgnc.hgnc_symbol.isna()].index, "hgnc_symbol"] = \
        ensembltohgnc[ensembltohgnc.hgnc_symbol.isna()
                      ]['clone_based_ensembl_gene']

    gene_rename = {i.ensembl_gene_id: i.hgnc_symbol+' ('+i.ensembl_gene_id+')'
                   for _, i in ensembltohgnc.iterrows()}
    protcod_rename = {
        i.ensembl_gene_id: i.hgnc_symbol+' ('+str(int(i.entrezgene_id))+')'
        for _, i in
        ensembltohgnc[(~ensembltohgnc.entrezgene_id.isna()) &
                      (ensembltohgnc.gene_biotype == 'protein_coding')].iterrows()}
    return gene_rename, protcod_rename, ensembltohgnc
