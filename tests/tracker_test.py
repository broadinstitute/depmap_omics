import pandas as pd

from depmapomics.ccle import expressionRenaming
from depmapomics.config import *
from depmapomics import tracker as track
from depmapomics import expressions
from depmapomics import ccle


def test_expression_outputs():
    trackerobj = track.initTracker()
    testset_name = "test_postprocessing"
    folder = "temp/dryrun"

    ccle_refsamples = trackerobj.read_tracker()

    import asyncio

    files, failed, _, renaming, lowqual, _ = asyncio.run(
        expressions.postProcess(
            RNAWORKSPACE,
            testset_name,
            save_output=folder,
            samplesetToLoad=testset_name,
            geneLevelCols=RSEMFILENAME_GENE,
            trancriptLevelCols=RSEMFILENAME_TRANSCRIPTS,
            compute_enrichment=False,
            ssGSEAcol="genes_tpm",
            renamingFunc=ccle.expressionRenaming,
            dry_run=True,
        )
    )

    # check output file formatting
    proteincoding_genes_tpm_logp1 = pd.read_csv(
        "temp/dryrun/proteincoding_genes_tpm_logp1.csv", index_col=0
    )

    assert all(
        [x[:4] == "ACH-" for x in proteincoding_genes_tpm_logp1.index]
    ), "proteincoding_genes_tpm_logp1: indices are not arxspans in"

    assert (
        ~proteincoding_genes_tpm_logp1.isnull().values.any()
    ), "proteincoding_genes_tpm_logp1: NANs found"

    # updated_tracker = expressions._CCLEPostProcessing(
    #     samplesetname=testset_name,
    #     trackerobj=trackerobj,
    #     samplesetToLoad=testset_name,
    #     dry_run=True,
    #     recompute_ssgsea=False,
    #     compute_enrichment=False,
    # )

    return True
