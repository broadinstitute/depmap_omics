import io
import os.path
import asyncio

import dalmatian as dm
from numpy.lib.npyio import save
import pandas as pd
import numpy as np
from scipy.stats import zscore 

from biomart import BiomartServer
from genepy.google import gcp
from genepy.utils import helper as h
from genepy.utils.helper import createFoldersFor
from genepy.google.google_sheet import dfToSheet
from genepy import rna, terra

from depmapomics import terra as myterra
from depmapomics.qc import rna as myQC
from depmapomics import tracker

from depmapomics.config import CACHE_PATH, TMP_PATH, ENSEMBL_SERVER_V

from gsheets import Sheets
from taigapy import TaigaClient
tc = TaigaClient()


def solveQC(tracker, failed, save=""):
    """
    # TODO todocument
    """
    newfail = []
    rename = {}
    # finding other replicates to solve failed ones
    for val in failed:
        a = tracker.loc[val].arxspan_id
        res = tracker[(tracker.datatype == 'rna')
                            & (tracker.arxspan_id == a)]
        if len(res) > 1:
            for k in res.index:
                if k not in failed:
                    rename[val] = k
        else:
            newfail.append(val)
    print("samples that failed:")
    print(newfail)
    if save:
        h.listToFile(newfail, save+"_rnafailed.txt")
    return rename
    

def updateTracker(refworkspace, selected, lowqual, tracker, samples, samplesetname, 
                  sheetname='ccle sample tracker', sheetcreds='../.credentials.json', onlycol=[
                    "star_bam_file", 'star_bam_index'], newgs="gs://cclebams/rnasq_hg38/",
                  dry_run=False, keeppath=False, qcname="star_logs", match=".Log.final.out"):
    """
    # TODO: to document
    """
    
    starlogs = myterra.getQC(workspace=refworkspace, only=samples,
                             qcname=qcname, match=match)
    for k, v in starlogs.items():
        if k == 'nan':
            continue
        if tracker.loc[k, 'bam_qc'] != v[0]:
            tracker.loc[k, 'bam_qc'] = v[0]
    
    ## copy star bam file to our cclebams/rnasq_hg38/ bucket
    renamed, _ = terra.changeGSlocation(workspacefrom=refworkspace, newgs=newgs, 
                                        onlysamples=samples, onlycol=onlycol, 
                                        entity="sample", keeppath=keeppath, dry_run=dry_run)

    tracker.loc[samples, ['legacy_size', 'legacy_crc32c_hash']
                        ] = tracker.loc[samples][['size', 'crc32c_hash']].values
    tracker.loc[samples, ['internal_bam_filepath', "internal_bai_filepath"]] = renamed[[
        'star_bam_file', 'star_bam_index']].values
    tracker.loc[samples, 'size'] = [gcp.extractSize(
        i)[1] for i in gcp.lsFiles(renamed['star_bam_file'].tolist(), '-l')]
    tracker.loc[samples, 'crc32c_hash'] = [gcp.extractHash(
        i) for i in gcp.lsFiles(renamed['star_bam_file'].tolist(), '-L')]
    tracker.loc[samples, 'md5_hash'] = [gcp.extractHash(
        i, "md5") for i in gcp.lsFiles(renamed['star_bam_file'].tolist(), '-L')]
    
    tracker.loc[selected, samplesetname]=1
    tracker.loc[tracker[tracker.datatype=='rna'].index,'low_quality']=0
    tracker.loc[lowqual,'low_quality']=1
    dfToSheet(tracker, sheetname, secret=sheetcreds)


def loadFromRSEMaggregate(refwm, todrop=[], filenames=["transcripts_tpm", "genes_tpm", 
                                                        "genes_expected_count", 
                                                        "transcripts_expected_count"],
                          sampleset="all", renamingFunc=None):
    """
    #TODO: to document
    Args:
    ----
  
    """
    files = {}
    renaming = {}
    refwm = dm.WorkspaceManager(refwm)
    res = refwm.get_sample_sets().loc[sampleset]
    for val in filenames:
        file = pd.read_csv(res[val], compression='gzip', header=0,
                        sep='\t', quotechar='"', error_bad_lines=False)
        if renamingFunc is not None:
            # removing failed version
            renaming = renamingFunc(file.columns[2:], todrop)
        else:
            renaming.update({i:i for i in file.columns[2:] if i not in todrop})
        renaming.update({'transcript_id(s)': 'transcript'})
        # we remove the failed samples where we did not found anything else to replace them with
        files[val] = file[file.columns[:2].tolist()+[i for i in file.columns[2:]
                            if i in set(renaming.keys())]].rename(columns=renaming)
    return files, renaming


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
        ensembltohgnc[ensembltohgnc.hgnc_symbol.isna()]['clone_based_ensembl_gene']

    gene_rename = {i.ensembl_gene_id: i.hgnc_symbol+' ('+i.ensembl_gene_id+')'
                   for _, i in ensembltohgnc.iterrows()}
    protcod_rename = {
        i.ensembl_gene_id: i.hgnc_symbol+' ('+str(int(i.entrezgene_id))+')'
        for _, i in
        ensembltohgnc[(~ensembltohgnc.entrezgene_id.isna()) &
                      (ensembltohgnc.gene_biotype == 'protein_coding')].iterrows()}
    return gene_rename, protcod_rename, ensembltohgnc


def subsetGenes(files, gene_rename, filenames=['rsem_transcripts_expected_count',
                                                'rsem_transcripts_tpm'],
                              drop=[], index="transcript_id"):
    """
    # TODO: to document
    """
    print('subsetting transcript columns')
    rename_transcript = {}
    missing = []
    for val in filenames:
        file = files[val].drop(columns=drop).set_index(index)
        file = file[(file.sum(1) != 0) & (file.var(1) != 0)]
        r = [i.split('.')[0] for i in file.index]
        dup = h.dups(r)
        if len(dup) > 0:
            print(dup)
            raise ValueError('duplicate genes')
        file.index = r
        if len(rename_transcript) == 0 and index == "transcript_id":
            for _, v in files[val].iterrows():
                if v.gene_id.split('.')[0] in gene_rename:
                    rename_transcript[v.transcript_id] = gene_rename[
                        v.gene_id.split('.')[0]].split(' (')[0] + ' (' + v.transcript_id + ')'
                else:
                    missing.append(v.gene_id.split('.')[0])
            print('missing: '+str(len(missing))+' genes')
        files[val].rename(index=rename_transcript if len(
            rename_transcript) != 0 else gene_rename).T
    return files


def extractProtCod(files, ensembltohgnc, protcod_rename,
                    filenames=['rsem_transcripts_expected_count', 'rsem_transcripts_tpm'],
                    rep=('genes', 'proteincoding_genes')):
    """
    # TODO: to document
    """
    for val in filenames:
        name = val.replace(rep[0], rep[1])
        files[name] = files[name][files[name].columns.isin(set(ensembltohgnc.ensembl_gene_id))]
        files[name] = files[name].rename(columns=protcod_rename)
        # removing genes that did not match.. pretty unfortunate
        files[name] = files[name][[i for i in files[name].columns if ' (' in i]]
        # check: do we have any duplicates?
        # if we do, managing duplicates
        if len(set(h.dups(files[name].columns.tolist()))) > 0:
            print("we have duplicate gene names!!")
            for dup in h.dups(files[name].columns):
                a = files[name][dup].sum()
                files[name].drop(columns=dup)
                files[name][dup] = a
            
    return files


async def ssGSEA(tpm_genes, pathtogenepy="../",
                 geneset_file="data/genesets/msigdb.v7.2.symbols.gmt"):
    """
    #TODO: todocument
    """
    tpm_genes = tpm_genes.copy()
    tpm_genes.columns = [i.split(' (')[0] for i in tpm_genes.columns]

    # summing the different exons/duplicates
    for i in h.dups(tpm_genes.columns):
        val = tpm_genes[i].sum(1)
        tpm_genes = tpm_genes.drop(columns=i)
        tpm_genes[i] = val

    # total size of data
    print("total size of data")
    print(len(set([val for val in tpm_genes.columns if '.' not in val])))
    tpm_genes = pd.DataFrame(data=zscore(np.log2(
        tpm_genes+1), nan_policy="omit"), columns=tpm_genes.columns, index=tpm_genes.index)

    # MAYBE NOT NEEDED
    #### merging splicing variants into the same gene
    #counts_genes_merged, _, _= h.mergeSplicingVariants(counts_genes.T, defined='.')

    enrichments = (await rna.gsva(tpm_genes.T, pathtogenepy=pathtogenepy,
                                  geneset_file=geneset_file, method='ssgsea')).T
    enrichments.index = [i.replace('.', '-') for i in enrichments.index]
    return enrichments


def saveFiles(files, release, folder=TMP_PATH):
    """
    # TODO: to document
    """
    print('storing files in {}'.format(folder))
    for k, val in files.items():
        val.to_csv(os.path.join(folder, k.replace(
            'rsem', 'expression_'+release)+'.csv'))
        if 'tpm' in k:
            val.apply(lambda x: np.log2(x + 1)).to_csv(os.path.join(folder,
                                                   k.replace('rsem', 'expression_' + release) + 
                                                   '_logp1.csv'))


def postprocessExpression(refworkspace, samplesetname,
                          save_output="", doCleanup=False,
                          colstoclean=[], ensemblserver=ENSEMBL_SERVER_V,
                          previousQCfail=[], samplesetToLoad="all",
                          geneLevelCols=[ "genes_tpm", "genes_expected_count"],
                          trancriptLevelCols=["transcripts_tpm", "transcripts_expected_count"],
                          ssGSEAcol="gene_tpm", renamingFunc=None, useCache=False
                          ):
    """
    # TODO: to document
    """
    if not samplesetToLoad:
        samplesetToLoad = samplesetname
    refwm = dm.WorkspaceManager(refworkspace)
    print("load QC and generate QC report")
    samplesinset = [i['entityName'] for i in refwm.get_entities(
        'sample_set').loc[samplesetname].samples]

    _, lowqual, failed = myQC.plot_rnaseqc_results( refworkspace, samplesinset,
                                                      save=bool(save_output),  
                                                      output_path=save_output+"/"+
                                                      samplesetname+"_rna_qcs/")

    failed = failed.index.tolist()
    print('you want to copy that up top, to save it for next time', failed)
    failed.extend(previousQCfail)
    
    if doCleanup:
        print("cleaninp up data")
        res = refwm.get_samples()
        for val in colstoclean:
            refwm.disable_hound().delete_entity_attributes(
                'sample', res[val], delete_files=True)
    
    print("generating gene names")
    gene_rename, protcod_rename, ensembltohgnc = generateGeneNames(
        ensemble_server=ensemblserver, useCache=useCache)
    
    print("loading files")
    files, renaming = loadFromRSEMaggregate(refwm, todrop=failed, filenames=trancriptLevelCols+geneLevelCols, 
                                  sampleset=samplesetToLoad, renamingFunc=renamingFunc)
    if save_output:
        h.dictToFile(renaming, save_output+"_sample_renaming.json")
        lowqual.to_csv(save_output+"_lowqual_samples.csv")
        h.listToFile(failed, save_output+"_failed_samples.txt")
    print("renaming files")
    # gene level
    if len(geneLevelCols) > 0:
        files = subsetGenes(files, gene_rename, filenames=geneLevelCols, 
                                    drop="transcript_id", index="gene_id")

        files = extractProtCod(files, ensembltohgnc[ensembltohgnc.gene_biotype == 'protein_coding'], 
                                protcod_rename, filenames=geneLevelCols)
    if len(trancriptLevelCols) > 0:
        files = subsetGenes(
            files, gene_rename, filenames=trancriptLevelCols, drop="gene_id", index="transcript_id")
    
    print("doing ssGSEA")
    enrichments = asyncio.run(ssGSEA(files[ssGSEAcol]))
    print("saving files")
    enrichments.to_csv('temp/gene_sets_'+samplesetname+'_all.csv')
    saveFiles(files, samplesetname)
    return files, enrichments, failed, samplesinset, renaming
    

def CCLEPostProcessing(refworkspace, samplesetname, refsheet_url="https://docs.google.com/spreadsheets/d/1Pgb5fIClGnErEqzxpU7qqX6ULpGTDjvzWwDN8XUJKIY",
                       colstoclean=['fastq1', 'fastq2', 'recalibrated_bam', 'recalibrated_bam_index'], 
                       ensemblserver=ENSEMBL_SERVER_V, doCleanup=True,
                       previousQCfail=[], my_id='~/.client_secret.json',
                       mystorage_id="~/.storage.json", samplesetToLoad="all", 
                       tocompare={"genes_expected_count": "CCLE_RNAseq_reads",
                                  "genes_tpm": "CCLE_expression_full",
                                  "proteincoding_genes_tpm": "CCLE_expression"},
                       sheetname='ccle sample tracker', sheetcreds='../.credentials.json',
                       **kwargs):

    
    sheets = Sheets.from_files(my_id, mystorage_id)
    ccle_refsamples = sheets.get(refsheet_url).sheets[0].to_frame(index_col=0)
    
    def rn(r, todrop):
        renaming = tracker.removeOlderVersions(
            names=r, refsamples=ccle_refsamples[ccle_refsamples.datatype=="rna"], 
            arxspan_id="arxspan_id", version="version")
        # if we have a replaceable failed version in our dataset
        rename = solveQC(ccle_refsamples, todrop)
        for k, _ in renaming.items():
            if k in rename:
                renaming[rename[k]] = renaming.pop(k)
        return renaming

    folder = os.path.join("data", samplesetname, "rna_")
    files, enrichments, _, samplesinset, _ = postprocessExpression(refworkspace, samplesetname,
                                                        save_output=folder,
                            doCleanup=doCleanup,
                            colstoclean=colstoclean, ensemblserver=ensemblserver,
                            previousQCfail=previousQCfail, samplesetToLoad=samplesetToLoad,
                            geneLevelCols=["genes_tpm", "genes_expected_count"],
                            trancriptLevelCols=["transcripts_tpm",
                                                "transcripts_expected_count"],
                            ssGSEAcol="gene_tpm", renamingFunc=rn)
    
    print("doing validation")
    prevcounts = tc.get(name='depmap-a0ab', file='CCLE_RNAseq_reads')
    nonoverlap = set(prevcounts.columns) ^ set(
        files.get('genes_expected_count').columns)
    print("number of non overlaping genes:")
    print(len(nonoverlap))
    # have we lost any samples compared to last release?
    lost = set(prevcounts.index) - set(files.get('genes_expected_count').index)
    print("of which, lost genes:")
    print(lost)
    # do we have samples that are missanotated compared to previous releases (replicate level)
    #notindataset, missannotated, unmatched = findMissAnnotatedReplicates(replevel, prevcounts, renaming)
    #for k,v in unmatched.items():
    #    if ccle_refsamples.loc[k].arxspan_id!=v:
    #        print(k,v)
    # do we have samples that are missanotated compared to previous releases (sample level)
    unmatched = rna.getDifferencesFromCorrelations(
        files.get('genes_expected_count'), prevcounts, minsimi=0.95)
    print("differences in correlations against the previous release")
    print(unmatched)
    # Is it because of  duplicate version?
    print('do we see it as a duplicate in the tracker?')
    rnasamples = ccle_refsamples[ccle_refsamples.datatype == 'rna']
    for i, _ in unmatched:
        print(len(rnasamples[rnasamples.arxspan_id == i]))
    
    #CCLE_expression, CCLE_expression_full, ,
    print("comparing to previous release")
    #h.compareDfs(files["rsem_transcripts_tpm"], tc.get(name='depmap-a0ab', file='CCLE_RNAseq_transcripts'))
    #h.compareDfs(files["rsem_transcripts_expected_count"], tc.get(name='depmap-a0ab', file='CCLE_expression_transcripts_expected_count'))
    # h.compareDfs(enrichments, tc.get(name='depmap-a0ab', file='CCLE_fusions_unfiltered'))
    for key, val in tocompare.items():
        _, omissmatchCols, _, omissmatchInds, newNAs, new0s = h.compareDfs(
            files[key], tc.get(name='depmap-a0ab', file=val))
        print(key)
        assert omissmatchCols == 0
        assert omissmatchInds == 0
        assert newNAs == 0
        assert new0s == 0
    
    print("updating the tracker")
    renaming = h.fileToDict(folder+"_sample_renaming.json")
    lowqual = pd.read_csv(folder+"_lowqual_samples.csv", index_col=0)
    updateTracker(refworkspace, set(renaming.keys()) - set(['transcript_id(s)']),
                lowqual[lowqual.sum(1) > 3].index.tolist(),
                ccle_refsamples, samplesinset, samplesetname, sheetname=sheetname, sheetcreds=sheetcreds)
    print("done")
