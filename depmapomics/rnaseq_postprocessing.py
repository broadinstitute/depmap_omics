import functools
import io
import os
import os.path

import dalmatian as dm
import pandas as pd
from biomart import BiomartServer
from genepy.google.gcp import cpFiles

from depmapomics.config import CACHE_PATH
from depmapomics.omics_pipeline import OmicsPipeline
from genepy.utils import helper as h

class RNAseq(OmicsPipeline):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.gene_rename, self.protcod_rename, self.ensembltohgnc = self.generate_gene_names()

    def download_rnaseq_files_to_tmp(self, sample_set_id):
        res = self.samplesets.loc[sample_set_id]
        cpFiles([res['rsem_genes_expected_count']], "{}/rsem_genes_expected_count".format(self.tmp_path))
        cpFiles([res['rsem_genes_tpm']], "{}/rsem_genes_tpm".format(self.tmp_path))
        cpFiles([res['rsem_transcripts_tpm']], "{}/rsem_transcripts_tpm".format(self.tmp_path))
        cpFiles([res['rsem_transcripts_expected_count']],
            "{}/rsem_transcripts_expected_count".format(self.tmp_path))

    def generate_gene_names(self, cached=True):
        cachefile = os.path.join(CACHE_PATH, 'biomart_ensembltohgnc.csv')
        cachefile = os.path.expanduser(cachefile)
        if cached & os.path.isfile(cachefile):
            self.printv('fetching gene names from biomart cache')
            ensembltohgnc = pd.read_csv(cachefile)
        else:
            self.printv('downloading gene names from biomart')
            server = BiomartServer("http://www.ensembl.org/biomart")
            ensmbl = server.datasets['hsapiens_gene_ensembl']
            ensembltohgnc = pd.read_csv(io.StringIO(ensmbl.search({
              'attributes': ['ensembl_gene_id','clone_based_ensembl_gene','hgnc_symbol','gene_biotype','entrezgene_id']
            }, header=1).content.decode()), sep='\t')
            ensembltohgnc.to_csv(cachefile, index=False)

        ensembltohgnc.columns = ['ensembl_gene_id','clone_based_ensembl_gene','hgnc_symbol','gene_biotype','entrezgene_id']
        ensembltohgnc = ensembltohgnc[~(ensembltohgnc['clone_based_ensembl_gene'].isna() & ensembltohgnc['hgnc_symbol'].isna())]
        ensembltohgnc.loc[ensembltohgnc[ensembltohgnc.hgnc_symbol.isna()].index,"hgnc_symbol"] = ensembltohgnc[ensembltohgnc.hgnc_symbol.isna()]['clone_based_ensembl_gene']

        gene_rename =  {i.ensembl_gene_id: i.hgnc_symbol+' ('+i.ensembl_gene_id+')' for k,i in ensembltohgnc.iterrows()}
        protcod_rename = {i.ensembl_gene_id: i.hgnc_symbol+' ('+str(int(i.entrezgene_id))+')' for _,i in ensembltohgnc[(~ensembltohgnc.entrezgene_id.isna()) & (ensembltohgnc.gene_biotype=='protein_coding')].iterrows()}
        return gene_rename, protcod_rename, ensembltohgnc

    def cleanup_columns(self):
        self.printv('cleaning up columns')
        rsem_files = [
          "{}/rsem_transcripts_tpm".format(self.tmp_path),
          "{}/rsem_genes_tpm".format(self.tmp_path),
          "{}/rsem_genes_expected_count".format(self.tmp_path),
          "{}/rsem_transcripts_expected_count".format(self.tmp_path)]
        files = {}
        for val in rsem_files:
            file = pd.read_csv(val, compression='gzip', header=0,
                              sep='\t', quotechar='"', error_bad_lines=False)
            file.rename(columns = {'transcript_id(s)': 'transcript_id'}, inplace=True)
            if val in [x for x in rsem_files if 'gene' in x]:
                if ('gene_id' not in file.columns) & ('Unnamed: 0' in file.columns):
                    file.rename(columns={'Unnamed: 0': 'gene_id'}, inplace=True)

                assert {'transcript_id', 'gene_id'} - set(file.columns) == set(), \
                    '{}:\n{}'.format(val, ', '.join(file.columns))
                    # '{}:\n{}'.format(val, file.head(2))
            #     file = file.drop(columns = 'transcript_id').set_index('gene_id', drop=True)
            #     file = file[(file.sum(1) != 0) & (file.var(1) != 0)]
            files[val.split('/')[-1]] = file
        return files

    def subset_gene_columns(self, files):
        self.printv('subsetting gene columns')
        for val in ['rsem_genes_expected_count','rsem_genes_tpm']:
            file = files[val].drop(columns='transcript_id').set_index('gene_id')
            file = file[(file.sum(1) != 0) & (file.var(1) != 0)]
            r = [i.split('.')[0] for i in file.index]
            dup = h.dups(r)
            if len(dup)>0:
                print(dup)
                raise ValueError('duplicate genes')
            file.index = r
            files[val.replace('genes','proteincoding_genes')] = file[file.index.isin(
              set(self.ensembltohgnc[
                self.ensembltohgnc.gene_biotype == 'protein_coding'].ensembl_gene_id))]
            files[val] = file.rename(index=self.gene_rename).T
            files[val.replace('genes','proteincoding_genes')] = files[
              val.replace('genes','proteincoding_genes')].rename(index=self.protcod_rename).T
        return files

    def subset_transcript_columns(self, files):
        self.printv('subsetting transcript columns')
        rename_transcript = {}
        missing = []
        for val in ['rsem_transcripts_expected_count','rsem_transcripts_tpm']:
            file = files[val]
            file = file[(file[file.columns[2:]].sum(1) != 0) & (file[file.columns[2:]].var(1) != 0)]
            r = [i.split('.')[0] for i in file.transcript_id]
            dup = h.dups(r)
            if len(dup)>0:
                print(dup)
                raise ValueError('duplicate genes')
            file.loc[:, 'transcript_id'] = r
            if len(rename_transcript)==0:
                for _,v in file.iterrows():
                    if v.gene_id.split('.')[0] in self.gene_rename:
                        rename_transcript[v.transcript_id] = self.gene_rename[v.gene_id.split('.')[0]].split(' (')[0] + ' (' + v.transcript_id + ')'
                    else:
                        missing.append(v.gene_id.split('.')[0])
                print('missing: '+str(len(missing))+' genes')
            files[val] = file.set_index('transcript_id').drop(columns = 'gene_id').rename(index = rename_transcript).T
        return files

    def format_rsem_results(self):
        files = self.cleanup_columns()
        files = self.subset_gene_columns(files)
        files = self.subset_transcript_columns(files)
        files = self.remove_genes_that_do_not_match(files)
        return files

    def remove_genes_that_do_not_match(self, files):
        self.printv('removing genes that do not match')
        files['rsem_proteincoding_genes_expected_count'] = files[
          'rsem_proteincoding_genes_expected_count'][[
            i for i in files['rsem_proteincoding_genes_expected_count'].columns if ' (' in i]]
        files['rsem_proteincoding_genes_tpm'] = files[
          'rsem_proteincoding_genes_tpm'][[
            i for i in files['rsem_proteincoding_genes_tpm'].columns if ' (' in i]]
        return files

    def save_files(self, files, release):
        self.printv('storing files in {}'.format(self.tmp_path))
        for k,val in files.items():
            val.to_csv(os.path.join(self.tmp_path, k.replace('rsem','expression_'+release)+'.csv'))

    def postprocess(self, sample_set, save_output=True, release=None):
        if release is None:
            release = sample_set
        self.download_rnaseq_files_to_tmp(sample_set)
        files = self.format_rsem_results()
        if save_output:
            self.save_files(files, release=release)
        return files
