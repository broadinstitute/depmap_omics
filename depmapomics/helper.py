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

class depmap_igv(igv.Browser):
    def __init__(self, config):
        super().__init__({**config, "oauthToken": get_gcloud_auth_token()})

    def get_the_track(self, depmap_sample, datatype, iloc_value=0, column_type='legacy', fileformat='bam'):
        sample_info = depmap_sample[depmap_sample['datatype'] == datatype].iloc[iloc_value]
        self.load_track(
        {
            "name": '{} ({}): {}'.format(sample_info['arxspan_id'],
                                         sample_info['stripped_cell_line_name'],
                                         sample_info['datatype']),
            "url": sample_info['{}_bam_filepath'.format(column_type)],
            'indexURL': sample_info['{}_bai_filepath'.format(column_type)],
            "format": fileformat
        })

    def get_the_tracks(self, depmap_sample, schemas):
        for schema in schemas:
            self.get_the_track(depmap_sample, **schema)
