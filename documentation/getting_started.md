This document summarizes how to set up Terra and get access to services that are required to run the pipeline.

## For Internal Users <a name="internal-users"></a>

### Getting Terra Access

1. You will need to request access to the following terra workspaces:
  - [RNASeq](https://app.terra.bio/#workspaces/broad-firecloud-ccle/DepMap_hg38_RNAseq)
  - [WES](https://app.terra.bio/#workspaces/broad-firecloud-ccle/DepMap_WES_CN_hg38)
  - [WGS (including the mutation pipeline)](https://app.terra.bio/#workspaces/broad-firecloud-ccle/DepMap_WGS_CN)

The current owners of these workspaces should give you access.

2. For the mutation pipeline you will also need to request dbGaP access (required for TCGA workflows). *See CCLE/new hiree section on Asana for details*.
3. Acquire access to required billing projects (e.g. broad-firecloud-ccle). *See CCLE/new hiree section on Asana for details*.
4. Get access to the following Terra groups:
  - depmap_ccle_data
  - depmap-pipelines
  - ccle-pipeline
5. If you need to get access to the data delivered by Genomics Platform (GP) at the Broad, use the following links:
  - __WES__ [IBM](https://app.terra.bio/#workspaces/terra-broad-cancer-prod/Getz_IBM_CellLines_Exomes)
  - __WES__ [Broad](https://app.terra.bio/#workspaces/terra-broad-cancer-prod/CCLE_DepMap_WES)
  - __RNA__ [IBM](https://app.terra.bio/#workspaces/terra-broad-cancer-prod/Getz_IBM_CellLines_RNASeqData)
  - __RNA__ [Broad](https://app.terra.bio/#workspaces/terra-broad-cancer-prod/CCLE_DepMap_RNAseq)
  - __WGS__ [IBM](https://app.terra.bio/#workspaces/terra-broad-cancer-prod/Getz_IBM_CellLines_WGS)
  - __WGS__ [Broad](https://app.terra.bio/#workspaces/terra-broad-cancer-prod/DepMap_WGS)
  - __WGS__ [Broad hg38 cram](https://app.terra.bio/#workspaces/broad-genomics-data/DepMap_WGS)
6. Request access to the data bucket `gs://cclebams/`


### Additional python dependencies:
- [Taiga](https://cds.team/taiga) is a platform that allows the Cancer Data Science team to store and share data. In order to access and upload data, you will need to login to [taiga](https://cds.team/taiga) with your broad google account and [set up your token](https://github.com/broadinstitute/taigapy#prerequisites) for the python client.
- We are currently using a relational database, Gumbo, to track our cell lines' metadata and release status. In order to interact with Gumbo through python, follow the instruction and install the Gumbo client [here](https://github.com/broadinstitute/gumbo_client).
- In order to use internal-only functions involved in the loading, depmap data post-processing, and uploading, you need to install the [depmap_omics_upload repo](https://github.com/broadinstitute/depmap_omics_upload)


## For External Users <a name="external-users"></a>

In order to run the processing pipelines written in [WDL](https://software.broadinstitute.org/wdl/documentation/), you will need to set up workspaces on Terra:

### Creating your Terra Workspaces:

1. You first need a [Terra](https://app.terra.bio/#) account with correct billing setup. See [here](https://support.terra.bio/hc/en-us/articles/360034677651-Account-setup-and-exploring-Terra) for a tutorial on getting started.
2. If you haven't already, create a workspace under the billing project of your choice. If you need to process both RNA and WGS data, we recommend creating one workspace for each. 
3. Import the WDL scripts by following the links to dockstore and clicking on *launch with terra* (note: you'll need both *_pipeline and *_aggregate for each data type):
    - __WGS__ [Pipeline](https://dockstore.org/workflows/github.com/broadinstitute/depmap_omics/WGS_pipeline)
    - __WGS__ [Aggregate](https://dockstore.org/workflows/github.com/broadinstitute/depmap_omics/WGS_aggregate)
    - __RNA__ [Pipeline](https://dockstore.org/workflows/github.com/broadinstitute/depmap_omics/RNA_pipeline)
    - __RNA__ [Aggregate](https://dockstore.org/workflows/github.com/broadinstitute/depmap_omics/RNA_aggregate)
4. DepMap's workspace configurations are saved after each data release under `data/`. We recommend using configurations from the latest quarter. For example, if the latest release is `21Q4`, you should be able to find the configurations in `data/21Q4/RNAconfig/all_configs.json` and `data/21Q4/RNAconfig/all_configs.json` for RNA and WGS, respectively.
5. Set up the right inputs and outpus for your workflows according to `inputs_[WORKFLOW_NAME].json` and `outputs_[WORKFLOW_NAME].json` files, which are under the same directory as `all_configs.json`.
6. Load your samples to the sample table so that their bam and bam index google storage filepaths get listed in the right data column to WGS_pipeline and RNA pipeline (e.g. internal_bam_filepath contains hg38 aligned bam files whereas hg19_bam_filepath contains hg19 aligned bam files). 
7. Create a sample set with the set of samples you want to analyse. Make sure the name of this sample set on terra is the same as `SAMPLESETNAME` in `config_global.py`.

Once this is done, you can launch your jupyter notebook server and run the `*_CCLE` jupyter notebooks corresponding to our RNA pipeline and WGS pipeline (older versions for WES (CN and mutations are available in a previous commit labelled 20Q4)).

Remark:
  1. you will need to use the `postProcesssing()` functions for post processing instead of the CCLE ones in the `dm_omics.py` module.
  2. you will need to change some of the variables in the `config_global.py` and `config_prod.py`.
  3. you won't be able to run the function conditional on the `isCCLE` boolean. You can however reimplement them to create your own pipeline.
