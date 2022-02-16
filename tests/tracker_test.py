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
    samplesintestset = [
        "CDS-NdhRa4",
        "CDS-p8j9qw",
        "CDS-27EI5P",
        "CDS-AVK2Vf",
        "CDS-NVVPN2",
    ]
    d = {
        "genes_tpm": "gs://fc-secure-8759d74a-3dc5-43df-9753-561478fd087d/06641010-77c3-430b-b706-8a626e035cce/RNA_aggregate/fabfecc3-b7e2-4448-b64f-2ef6d5fa63b2/call-rsem_aggregate_results/test_postprocessing.rsem_genes_tpm.txt.gz",
        "genes_expected_count": "gs://fc-secure-8759d74a-3dc5-43df-9753-561478fd087d/06641010-77c3-430b-b706-8a626e035cce/RNA_aggregate/fabfecc3-b7e2-4448-b64f-2ef6d5fa63b2/call-rsem_aggregate_results/test_postprocessing.rsem_genes_expected_count.txt.gz",
        "transcripts_tpm": "gs://fc-secure-8759d74a-3dc5-43df-9753-561478fd087d/06641010-77c3-430b-b706-8a626e035cce/RNA_aggregate/fabfecc3-b7e2-4448-b64f-2ef6d5fa63b2/call-rsem_aggregate_results/test_postprocessing.rsem_transcripts_tpm.txt.gz",
        "transcripts_expected_count": "gs://fc-secure-8759d74a-3dc5-43df-9753-561478fd087d/06641010-77c3-430b-b706-8a626e035cce/RNA_aggregate/fabfecc3-b7e2-4448-b64f-2ef6d5fa63b2/call-rsem_aggregate_results/test_postprocessing.rsem_transcripts_expected_count.txt.gz",
    }
    rsemfilelocs = pd.Series(
        data=d,
        index=[
            "genes_tpm",
            "genes_expected_count",
            "transcripts_tpm",
            "transcripts_expected_count",
        ],
    )
    rnaqclocs = {
        "CDS-27EI5P": [
            "gs://fc-secure-639c94ba-2b0d-4960-92fc-9cd50046a968/ed62eed5-989b-405c-ad3c-4b84a014527e/RNA_pipeline/c55cb752-5e51-4b9d-a005-ced21a135db5/call-rnaseqc2/CDS-27EI5P.metrics.tsv"
        ],
        "CDS-AVK2Vf": [
            "gs://fc-secure-639c94ba-2b0d-4960-92fc-9cd50046a968/ed62eed5-989b-405c-ad3c-4b84a014527e/RNA_pipeline/4981c4d0-d3cc-4c77-be3a-767e5b29a803/call-rnaseqc2/CDS-AVK2Vf.metrics.tsv"
        ],
        "CDS-NdhRa4": [
            "gs://fc-secure-639c94ba-2b0d-4960-92fc-9cd50046a968/ed62eed5-989b-405c-ad3c-4b84a014527e/RNA_pipeline/72a26499-02d7-4661-ae37-1e18499a10f2/call-rnaseqc2/CDS-NdhRa4.metrics.tsv"
        ],
        "CDS-NVVPN2": [
            "gs://fc-secure-639c94ba-2b0d-4960-92fc-9cd50046a968/ed62eed5-989b-405c-ad3c-4b84a014527e/RNA_pipeline/dde5a72b-84ef-4e48-ad3a-170c4addd93e/call-rnaseqc2/CDS-NVVPN2.metrics.tsv"
        ],
        "CDS-p8j9qw": [
            "gs://fc-secure-639c94ba-2b0d-4960-92fc-9cd50046a968/ed62eed5-989b-405c-ad3c-4b84a014527e/RNA_pipeline/a1e4d474-e59c-4f65-9400-d937d701a18e/call-rnaseqc2/CDS-p8j9qw.metrics.tsv"
        ],
    }

    ccle_refsamples = trackerobj.read_tracker()

    import asyncio

    files, failed, _, renaming, lowqual, _ = asyncio.run(
        expressions.postProcess(
            None,
            testset_name,
            save_output=folder,
            samplesetToLoad=testset_name,
            geneLevelCols=RSEMFILENAME_GENE,
            trancriptLevelCols=RSEMFILENAME_TRANSCRIPTS,
            compute_enrichment=False,
            ssGSEAcol="genes_tpm",
            renamingFunc=ccle.expressionRenaming,
            dry_run=True,
            samplesinset=samplesintestset,
            rsemfilelocs=rsemfilelocs,
            rnaqclocs=rnaqclocs,
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
