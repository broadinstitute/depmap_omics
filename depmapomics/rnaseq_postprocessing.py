from genepy.google.gcp import cpFiles
from biomart import BiomartServer
import io
import pandas as pd

def download_rnaseq_files_to_tmp(workspace_manager, sample_set_id):
    res = workspace_manager.get_sample_sets().loc[sample_set_id]
    rsem_genes_expected_count = res['rsem_genes_expected_count']
    rsem_genes_tpm = res['rsem_genes_tpm']
    rsem_transcripts_tpm = res['rsem_transcripts_tpm']
    rsem_transcripts_expected_count = res['rsem_transcripts_expected_count']
    cpFiles([rsem_genes_expected_count], "/tmp/rsem_genes_expected_count" )
    cpFiles([rsem_genes_tpm], "/tmp/rsem_genes_tpm")
    cpFiles([rsem_transcripts_tpm], "/tmp/rsem_transcripts_tpm")
    cpFiles([rsem_transcripts_expected_count], "/tmp/rsem_transcripts_expected_count")

def generate_gene_names():
    server = BiomartServer( "http://www.ensembl.org/biomart" )
    ensmbl = server.datasets['hsapiens_gene_ensembl']
    ensembltohgnc = pd.read_csv(io.StringIO(ensmbl.search({
      'attributes': ['ensembl_gene_id','clone_based_ensembl_gene','hgnc_symbol','gene_biotype','entrezgene_id']
    }, header=1).content.decode()), sep='\t')

    ensembltohgnc.columns = ['ensembl_gene_id','clone_based_ensembl_gene','hgnc_symbol','gene_biotype','entrezgene_id']
    ensembltohgnc = ensembltohgnc[~(ensembltohgnc['clone_based_ensembl_gene'].isna() & ensembltohgnc['hgnc_symbol'].isna())]
    ensembltohgnc.loc[ensembltohgnc[ensembltohgnc.hgnc_symbol.isna()].index,"hgnc_symbol"] = ensembltohgnc[ensembltohgnc.hgnc_symbol.isna()]['clone_based_ensembl_gene']

    gene_rename =  {i.ensembl_gene_id: i.hgnc_symbol+' ('+i.ensembl_gene_id+')' for k,i in ensembltohgnc.iterrows()}
    protcod_rename = {i.ensembl_gene_id: i.hgnc_symbol+' ('+str(int(i.entrezgene_id))+')' for _,i in ensembltohgnc[(~ensembltohgnc.entrezgene_id.isna()) & (ensembltohgnc.gene_biotype=='protein_coding')].iterrows()}
    return gene_rename, protcod_rename
