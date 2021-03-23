import sys

import pandas as pd
from genepy import rna
from src.CCLE_postp_function import getQC
import igv
from genepy.utils.helper import parrun

from depmapomics.config import (RNASEQC_THRESHOLDS_FAILED,
                                RNASEQC_THRESHOLDS_LOWQUAL)

def plot_rnaseqc_results(workspace, samplelist, output_path='data/rna_qc_plots/'):
    rnaqc = getQC(workspace=workspace, only=samplelist, qcname="rnaseqc2_metrics")
    assert pd.Series(rnaqc).map(lambda x: x[0]).notnull().all()

    import numpy as np
    qcs = pd.DataFrame()
    for k, val in rnaqc.items():
        if val[0] is not np.nan:
            qcs = pd.concat([qcs, pd.read_csv(val[0], sep='\t', index_col=0)], axis=1)
    qcs = qcs[~((qcs.mean(1) == 1.0) | (qcs.mean(1) == 0.0))]

    # pandas implementation
    # qcs = pd.Series(rnaqc).map(lambda x: x[0]).dropna()
    # qcs = qcs.apply(lambda x: pd.read_csv(x, sep='\t', index_col=0))
    # qcs = pd.concat(qcs.tolist(), axis=1)
    # qcs = qcs[~((qcs.mean(1)==1.0) | (qcs.mean(1)==0.0))]

    print('Low quality samples')

    sys.stdout.flush()

    lowqual = rna.filterRNAfromQC(qcs, thresholds=RNASEQC_THRESHOLDS_LOWQUAL,
                                  folder='{}/lowqual/'.format(output_path), plot=True,
                                  qant1=0.1, qant3=0.9)
    print('Failed QC samples')
    sys.stdout.flush()
    failed = rna.filterRNAfromQC(qcs, thresholds=RNASEQC_THRESHOLDS_FAILED,
                                 folder='{}/lowqual/'.format(output_path), plot=True,
                                 qant1=0.07, qant3=0.93)

    return qcs, lowqual, failed


def get_gcloud_auth_token():
    parrun(["echo $(gcloud auth application-default print-access-token) > /tmp/gcloud_token"], cores=1)

    with open('/tmp/gcloud_token', 'r') as f:
        token = f.readline().strip()
    return token

from gsheets import Sheets

def load_sample_tracker():
    SAMPLE_TRACKER_URL = 'https://docs.google.com/spreadsheets/d/1Pgb5fIClGnErEqzxpU7qqX6ULpGTDjvzWwDN8XUJKIY'
    sheets_obj = Sheets.from_files('~/.client_secret.json', '~/.storage.json')

    sheets = sheets_obj.get(SAMPLE_TRACKER_URL).sheets
    depmap_samples = sheets[0].to_frame(header=0, index_col=0)

    depmap_samples.loc[depmap_samples['datatype'].isin(['wes', 'rna', 'wgs']),
                       'hg19_bam'] = depmap_samples['legacy_bam_filepath']
    depmap_samples.loc[depmap_samples['datatype'].isin(['wes', 'rna', 'wgs']),
                       'hg19_bai'] = depmap_samples['legacy_bai_filepath']
    depmap_samples.loc[depmap_samples['datatype'].isin(['RRBS', 'raindance', 'hybrid_capture']),
                       'hg19_bam'] = depmap_samples['internal_bam_filepath']
    depmap_samples.loc[depmap_samples['datatype'].isin(['RRBS', 'raindance', 'hybrid_capture']),
                       'hg19_bai'] = depmap_samples['internal_bai_filepath']

    depmap_samples.loc[depmap_samples['datatype'].isin(['wes', 'rna', 'wgs']),
                       'hg38_bam'] = depmap_samples['internal_bam_filepath']
    depmap_samples.loc[depmap_samples['datatype'].isin(['wes', 'rna', 'wgs']),
                       'hg38_bai'] = depmap_samples['internal_bai_filepath']


    return depmap_samples


def subset_depmap_samples(depmap_samples, arxspan_ids, genome, datatypes=None):
    depmap_samples_subset = depmap_samples[depmap_samples['arxspan_id'].isin(arxspan_ids)]
    depmap_samples_subset = depmap_samples_subset.dropna(subset=['{}_bam'.format(genome), '{}_bai'.format(genome)])
    if datatypes is not None:
        depmap_samples_subset = depmap_samples_subset[depmap_samples_subset['datatype'].isin(datatypes)]
    return depmap_samples_subset

class depmap_igv(igv.Browser):
    '''
    this class can be used to plot IGV tracks for DepMap cell lines

    example code:
        b = depmap_igv({"genome": "hg19", 'locus':'11:62296943'})
        b.get_the_tracks(['ACH-000090', 'ACH-002184'], ['wes', 'wgs', 'RRBS', 'rna'])
        b.show()
    '''
    def __init__(self, config):
        super().__init__({**config, "oauthToken": get_gcloud_auth_token()})
        self.depmap_samples = load_sample_tracker()

    def get_the_track(self, row, fileformat='bam'):
        sample_info = row[1]
        self.load_track(
        {
            "name": '{} ({}): {}'.format(sample_info['arxspan_id'],
                                         sample_info['stripped_cell_line_name'],
                                         sample_info['datatype']),
            "url": sample_info['{}_bam'.format(self.config['genome'])],
            'indexURL': sample_info['{}_bai'.format(self.config['genome'])],
            "format": fileformat
        })

    def get_the_tracks_from_dataframe(self, depmap_samples):
        for row in depmap_samples.sort_values(['arxspan_id', 'datatype']).iterrows():
            self.get_the_track(row)

    def get_the_tracks(self, arxspan_ids, datatypes=None):
        '''
        create an IGV plot with tracks corresponding to DepMap data
        arxspan_ids: list of arxspan IDs
        datatypes: list of datatypes to show. If None will show all (default is None)
        '''
        depmap_samples_subset = subset_depmap_samples(self.depmap_samples, arxspan_ids,
                                                      self.config['genome'], datatypes=datatypes)
        self.get_the_tracks_from_dataframe(depmap_samples_subset)
