from genepy import rna
from depmapomics.config import RNASEQC_THRESHOLDS_LOWQUAL, RNASEQC_THRESHOLDS_FAILED

import pandas as pd

def plot_rnaseqc_results(workspace, samplelist, output_path = 'data/rna_qc_plots/'):
    rnaqc = getQC(workspace=workspace, only=samplelist, qcname="rnaseqc2_metrics")
    assert pd.Series(rnaqc).map(lambda x: x[0]).notnull().all()

    import numpy as np
    qcs = pd.DataFrame()
    for k,val in rnaqc.items():
        if val[0] is not np.nan:
            qcs = pd.concat([qcs, pd.read_csv(val[0],sep='\t', index_col=0)],axis=1)
    qcs = qcs[~((qcs.mean(1)==1.0) | (qcs.mean(1)==0.0))]

    # qcs = pd.Series(rnaqc).map(lambda x: x[0]).dropna()
    # qcs = qcs.apply(lambda x: pd.read_csv(x, sep='\t', index_col=0))
    # qcs = pd.concat(qcs.tolist(), axis=1)
    # qcs = qcs[~((qcs.mean(1)==1.0) | (qcs.mean(1)==0.0))]

    print('Low quality samples')
    import sys
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
