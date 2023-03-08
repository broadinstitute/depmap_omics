import os
import pandas as pd
import subprocess
from mgenepy.utils import helper as h
from matplotlib import pyplot as plt
import numpy as np
import math
import seaborn as sns


async def gsva(data, geneset_file, method="ssgsea", recompute=True):
    print("you need to have R installed with GSVA and GSEABase library installed")
    pathtogenepy = os.path.dirname(os.path.abspath(__file__))
    if (
        not recompute
        and os.path.exists("/tmp/data_genepyhelper_gsva.csv")
        and os.path.exists("/tmp/res_genepy_ssGSEA.tsv")
    ):
        print("trying to bypass computing...")
        v = pd.read_csv("/tmp/data_genepyhelper_gsva.csv", index_col=0)
        if v.shape[0] != data.shape[0] or v.shape[1] != data.shape[1]:
            print("WARNING: recompute to false but not the same df for sure")
        return pd.read_csv("/tmp/res_genepy_ssGSEA.tsv", sep="\t")
    data.to_csv("/tmp/data_genepyhelper_gsva.csv")
    cmd = (
        "Rscript "
        + pathtogenepy
        + "/ssGSEA.R /tmp/data_genepyhelper_gsva.csv "
        + geneset_file
        + " "
        + method
    )
    res = subprocess.run(cmd, shell=True, capture_output=True)
    if res.returncode != 0:
        raise ValueError("issue with the command: " + str(res))
    print(res)
    res = pd.read_csv("/tmp/res_genepy_ssGSEA.tsv", sep="\t")
    return res


def filterRNAfromQC(
    rnaqc,
    folder="tempRNAQCplot/",
    plot=True,
    qant1=0.07,
    qant3=0.93,
    thresholds={},
    num_cols=10,
    figsize=(10, 0.2),
):
    thresh = {
        "minmapping": 0.8,  # Mapping Rate
        "minendmapping": 0.75,
        "minefficiency": 0.6,  # Expression Profiling Efficiency
        "maxendmismatch": 0.025,  # Base Mismatch end wise
        "maxmismatch": 0.02,  # Base Mismatch
        "minhighqual": 0.6,  # High Quality Rate
        "minexon": 0.6,  # Exonic Rate
        "maxambiguous": 0.2,  # Ambiguous Alignment Rate
        "maxsplits": 0.1,  # Avg. Splits per Read
        "maxalt": 0.65,  # Alternative Alignments rate
        "maxchim": 0.3,  # Chimeric Alignment Rate
        "minreads": 20000000,
        "minlength": 80,  # Read Length
        "maxgenes": 35000,
        "mingenes": 10000,
    }
    thresh.update(thresholds)

    qcs = rnaqc.T
    tot = []
    a = qcs[
        (qcs["Mapping Rate"] < thresh["minmapping"])
        | (qcs["Base Mismatch"] > thresh["maxmismatch"])
        | (qcs["End 1 Mapping Rate"] < thresh["minendmapping"])
        | (qcs["End 2 Mapping Rate"] < thresh["minendmapping"])
        | (qcs["End 1 Mismatch Rate"] > thresh["maxendmismatch"])
        | (qcs["End 2 Mismatch Rate"] > thresh["maxendmismatch"])
        | (qcs["Expression Profiling Efficiency"] < thresh["minefficiency"])
        | (qcs["High Quality Rate"] < thresh["minhighqual"])
        | (qcs["Exonic Rate"] < thresh["minexon"])
        | (qcs["Ambiguous Alignment Rate"] > thresh["maxambiguous"])
        | (qcs["Avg. Splits per Read"] < thresh["maxsplits"])
        | (qcs["Alternative Alignments"] > thresh["maxalt"] * qcs["Total Reads"])
        | (qcs["Chimeric Alignment Rate"] > thresh["maxchim"])
        | (qcs["Total Reads"] < thresh["minreads"])
        | (qcs["Read Length"] < thresh["minlength"])
        | (thresh["maxgenes"] < qcs["Genes Detected"])
        | (qcs["Genes Detected"] < thresh["mingenes"])
    ].index.tolist()

    tot.append(
        [
            1
            if i in qcs[(qcs["Mapping Rate"] < thresh["minmapping"])].index.tolist()
            else 0
            for i in a
        ]
    )
    tot.append(
        [
            1
            if i in qcs[(qcs["Base Mismatch"] > thresh["maxmismatch"])].index.tolist()
            else 0
            for i in a
        ]
    )
    tot.append(
        [
            1
            if i
            in qcs[(qcs["End 1 Mapping Rate"] < thresh["minendmapping"])].index.tolist()
            else 0
            for i in a
        ]
    )
    tot.append(
        [
            1
            if i
            in qcs[(qcs["End 2 Mapping Rate"] < thresh["minendmapping"])].index.tolist()
            else 0
            for i in a
        ]
    )
    tot.append(
        [
            1
            if i
            in qcs[
                (qcs["End 1 Mismatch Rate"] > thresh["maxendmismatch"])
            ].index.tolist()
            else 0
            for i in a
        ]
    )
    tot.append(
        [
            1
            if i
            in qcs[
                (qcs["End 2 Mismatch Rate"] > thresh["maxendmismatch"])
            ].index.tolist()
            else 0
            for i in a
        ]
    )
    tot.append(
        [
            1
            if i
            in qcs[
                (qcs["Expression Profiling Efficiency"] < thresh["minefficiency"])
            ].index.tolist()
            else 0
            for i in a
        ]
    )
    tot.append(
        [
            1
            if i
            in qcs[(qcs["High Quality Rate"] < thresh["minhighqual"])].index.tolist()
            else 0
            for i in a
        ]
    )
    tot.append(
        [
            1
            if i in qcs[(qcs["Exonic Rate"] < thresh["minexon"])].index.tolist()
            else 0
            for i in a
        ]
    )
    tot.append(
        [
            1
            if i
            in qcs[
                (qcs["Ambiguous Alignment Rate"] > thresh["maxambiguous"])
            ].index.tolist()
            else 0
            for i in a
        ]
    )
    tot.append(
        [
            1
            if i
            in qcs[(qcs["Avg. Splits per Read"] < thresh["maxsplits"])].index.tolist()
            else 0
            for i in a
        ]
    )
    tot.append(
        [
            1
            if i
            in qcs[
                (qcs["Alternative Alignments"] > thresh["maxalt"] * qcs["Total Reads"])
            ].index.tolist()
            else 0
            for i in a
        ]
    )
    tot.append(
        [
            1
            if i
            in qcs[(qcs["Chimeric Alignment Rate"] > thresh["maxchim"])].index.tolist()
            else 0
            for i in a
        ]
    )
    tot.append(
        [
            1
            if i in qcs[(qcs["Total Reads"] < thresh["minreads"])].index.tolist()
            else 0
            for i in a
        ]
    )
    tot.append(
        [
            1
            if i in qcs[(qcs["Read Length"] < thresh["minlength"])].index.tolist()
            else 0
            for i in a
        ]
    )
    tot.append(
        [
            1
            if i in qcs[(thresh["maxgenes"] < qcs["Genes Detected"])].index.tolist()
            else 0
            for i in a
        ]
    )
    tot.append(
        [
            1
            if i in qcs[(qcs["Genes Detected"] < thresh["mingenes"])].index.tolist()
            else 0
            for i in a
        ]
    )

    res = pd.DataFrame(
        index=a,
        columns=[
            "Mapping Rate",
            "Base Mismatch",
            "End 1 Mapping Rate",
            "End 2 Mapping Rate",
            "End 1 Mismatch Rate",
            "End 2 Mismatch Rate",
            "Expression Profiling Efficiency",
            "High Quality Rate",
            "Exonic Rate",
            "Ambiguous Alignment Efficiency",
            "Avg. Splits per Read",
            "Alternative Alignments",
            "Chimeric Alignment Rate",
            "Total Reads",
            "Read Length",
            "Min Genes Detected",
            "Max Genes Detected",
        ],
        data=np.array(tot).astype(bool).T,
    )

    print(a)
    if len(res) > 0:
        h.createFoldersFor(folder)
        res.to_csv(folder + "_qc_results.csv")
        if plot:
            _, ax = plt.subplots(figsize=(figsize[0], math.ceil(len(res) * figsize[1])))
            plot = sns.heatmap(res, xticklabels=True, yticklabels=True, cbar=False)
            plt.yticks(rotation=0)
            plt.show()
            plot.get_figure().savefig(folder + "failed_qc.pdf")

            num_rows = math.ceil(len(rnaqc) / num_cols)
            _, axes = plt.subplots(num_rows, num_cols, figsize=(20, num_rows * 2))
            for val_idx, val in enumerate(rnaqc.index):
                ax = axes.flatten()[val_idx]
                qc = rnaqc.loc[val]
                sns.violinplot(y=qc, ax=ax)
                q1 = qc.quantile(qant1)
                q3 = qc.quantile(qant3)
                outlier_top_lim = q3 + 1.5 * (q3 - q1)
                outlier_bottom_lim = q1 - 1.5 * (q3 - q1)
                for k, v in qc[
                    (qc < outlier_bottom_lim) | (qc > outlier_top_lim)
                ].iteritems():
                    ax.text(
                        0.05,
                        v,
                        k,
                        ha="left",
                        va="center",
                        color="red" if k in a else "black",
                    )
            plt.tight_layout()
            plt.show()
            plt.savefig("{}/qc_metrics.pdf".format(folder), bbox_inches="tight")
    return res
