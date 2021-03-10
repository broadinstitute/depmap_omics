import dalmatian as dm
import pandas as pd
from genepy.google.gcp import cpFiles
from src.CCLE_postp_function import filterFusions
import re

def standardize_leftgene_rightgene(fusions):
    '''
    converts [GENE_NAME]^[ENSG] --> [GENE_NAME] ([ENSG])
    Example: "SMAD4^ENSG00000141646.14" --> "SMAD4 (ENSG00000141646.14)"
    '''
    fusions[['LeftGene', 'RightGene']] = fusions[['LeftGene', 'RightGene']]\
        .applymap(lambda x: '{} ({})'.format(*x.split(r'^')))
    return fusions

def postprocess_fusions(refworkspace, sample_id_name='DepMap_ID', sample_set_name = 'all',
                        output_fusion_file='/tmp/fusions.csv',
                        output_fusion_filtered_file='/tmp/filtered_fusions.csv'):
    refwm = dm.WorkspaceManager(refworkspace).disable_hound()
    aggregated = refwm.get_sample_sets().loc[sample_set_name]['fusions_star']
    cpFiles([aggregated], "/tmp/expression.fusion.tsv")
    fusions = pd.read_csv("/tmp/expression.fusion.tsv",
                          names=[sample_id_name, 'FusionName', 'JunctionReadCount',
                                 'SpanningFragCount', 'SpliceType', 'LeftGene', 'LeftBreakpoint',
                                 'RightGene', 'RightBreakpoint', 'LargeAnchorSupport', 'FFPM',
                                 'LeftBreakDinuc', 'LeftBreakEntropy', 'RightBreakDinuc',
                                 'RightBreakEntropy', 'annots'], skiprows=1, sep='\t')

    CCLE_count = fusions[['LeftBreakpoint', 'RightBreakpoint']]\
        .value_counts()\
        .to_frame(name='CCLE_count')
    fusions = pd.merge(fusions, CCLE_count, on=['LeftBreakpoint', 'RightBreakpoint'])

    fusions_filtered = filterFusions(fusions, maxfreq=0.1, sample_id_name=sample_id_name)

    fusions = standardize_leftgene_rightgene(fusions)
    fusions_filtered = standardize_leftgene_rightgene(fusions_filtered)

    fusions.to_csv(output_fusion_file, index=False)
    fusions_filtered.to_csv(output_fusion_filtered_file, index=False)

    return fusions, fusions_filtered
