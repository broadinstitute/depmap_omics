from depmapomics import constants
import sys

import pandas as pd
from mgenepy import rna
from depmapomics import terra
import numpy as np


def export_qc(
    workspace, selected_samples=[], colname="rnaseqc2_metrics", allow_missing=False
):
    """
    export qcs and aggregate them into one df
    """
    rnaqc = terra.getQC(workspace=workspace, only=selected_samples, qcname=colname)
    if allow_missing:
        missing_qc = [k for k, v in rnaqc.items() if len(v) == 0]
        print("The following samples don't have RNAseq QC in the workspace: ")
        print(missing_qc)
    else:
        assert (
            pd.Series(rnaqc).map(lambda x: x[0]).notnull().all()
        ), "Some samples have no QC data"
    qcs = pd.DataFrame()
    for _, val in rnaqc.items():
        if val[0] is not np.nan:
            qcs = pd.concat([qcs, pd.read_csv(val[0], sep="\t", index_col=0)], axis=1)
    return qcs


def plot_rnaseqc_results(
    workspace,
    samplelist,
    output_path="data/rna_qcs/",
    qcname="rnaseqc2_metrics",
    rnaqc={},
    save=True,
):
    """
    TODO: to document
    """
    qcs = export_qc(workspace, selected_samples=samplelist, colname=qcname)

    qcs = qcs[~((qcs.mean(1) == 1.0) | (qcs.mean(1) == 0.0))]

    # pandas implementation
    # qcs = pd.Series(rnaqc).map(lambda x: x[0]).dropna()
    # qcs = qcs.apply(lambda x: pd.read_csv(x, sep='\t', index_col=0))
    # qcs = pd.concat(qcs.tolist(), axis=1)
    # qcs = qcs[~((qcs.mean(1)==1.0) | (qcs.mean(1)==0.0))]

    print("Low quality samples")

    sys.stdout.flush()

    lowqual = rna.filterRNAfromQC(
        qcs,
        thresholds=constants.RNASEQC_THRESHOLDS_LOWQUAL,
        folder="{}/lowqual/".format(output_path),
        plot=save,
        qant1=0.1,
        qant3=0.9,
    )
    print("Failed QC samples")
    sys.stdout.flush()
    failed = rna.filterRNAfromQC(
        qcs,
        thresholds=constants.RNASEQC_THRESHOLDS_FAILED,
        folder="{}/failed/".format(output_path),
        plot=save,
        qant1=0.07,
        qant3=0.93,
    )

    return qcs, lowqual, failed
