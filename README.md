# depmap_omics

__have a look at [DepMap](https://www.depmap.org)__

![](documentation/depmap-logo_white.png)

This repository contains code that processes data for the biannual DepMap data release. 

[Here](https://docs.google.com/presentation/d/1i0HI31dBejTYmzI9Cp6Ij7--t6eSR2r9vcC22TxSnNI/edit#slide=id.g525fd14bef_0_116) is an overview of the pipeline and the data it produces.

## Table of Contents
- [Getting Started](#quickstart)
  - [Installation](#installation)
  - [For Internal Users](#internal-users)
  - [For External Users](#external-users)
- [Repository File Structure](#file-structure)
- [Running the Pipeline](#running-pipeline)
  - [Uploading and Preprocessing](#upload-preprocess)
  - [Running Terra Pipelines](#running-terra-pipelines)
  - [Downloading and Postprocessing](#downloading-postprocessing)
  - [QC, Groupding and Uploading](#qc-grouping-uploading)
- [Auxiliary Data for the Pipeline](#data)

## Getting Started <a name="quickstart"></a>

The processing pipeline rely on the following tools:

- [python](https://www.learnpython.org/)
- [R](https://www.codecademy.com/learn/learn-r)
- [jupyter](https://jupyter.org/index.html)
- [WDL](https://software.broadinstitute.org/wdl/documentation/)
- [gcp](https://cloud.google.com/sdk/docs/quickstart-macos)
- [docker](https://docs.docker.com/get-started/)
- [Terra](https://software.broadinstitute.org/firecloud/documentation/)
- [dalmatian](https://github.com/broadinstitute/dalmatian)
- [Terra and gcp](https://docs.google.com/document/d/1zTtaN-Px64f8JvgydZNdBbzBpFWyZzEpshSNxQh43Oc/edit#heading=h.dz5wh0l4bu9g)


### Installatiion <a name="installation"></a>

`git clone http://github.com/BroadInstitute/depmap_omics.git && cd depmap_omics`

`pip install -e .`

### :warning: this repository needs other repos

Some important data and code from the [genepy Library](https://github.com/broadinstitute/genepy).

Use the instructions in the genepy page to install the package.

### :warning: you would need the following R python packages

1. You will need to install jupyter notetbooks and google cloud sdk
  - install [Google Cloud SDK](https://cloud.google.com/sdk/docs/downloads-interactive).
  - authenticate my SDK account by running `gcloud auth application-default login` in the terminal, and follow the instrucion to log in.

2. and GSVA for ssGSEA in R `R` run `R -e 'if(!requireNamespace("BiocManager", quietly = TRUE)){install.packages("BiocManager")};BiocManager::install(c("GSEABase", "erccdashboard", "GSVA", "DESeq2"));'`

3. For Python use the requirements.txt file `pip install -r requirements.txt` 

## For Internal Users <a name="internal-users"></a>

> To learn about the tools we use in the pipeline, see [this section](#running-pipeline) for a detailed walkthrough

Before running the pipeline, make sure to set up the following:

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
6. Request access to the data bucket `gs://cclebams/`


### Additional python dependencies:
- [Taiga](https://cds.team/taiga) is a platform that allows the Cancer Data Science team to store and share data. In order to access and upload data, you will need to login to [taiga](https://cds.team/taiga) with your broad google account and [set up your token](https://github.com/broadinstitute/taigapy#prerequisites) for the python client.
- We are currently using a relational database, Gumbo, to track our cell lines' metadata and release status. In order to interact with Gumbo through python, follow the instruction and install the Gumbo client [here](https://github.com/broadinstitute/gumbo_client).


## For External Users <a name="external-users"></a>

> To learn about the tools we use in the pipeline, see [this section](#running-pipeline) for a detailed walkthrough

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


## Repository File Structure <a name="file-structure"></a>

For our 2 computation pipelines for depmap omics:
- Expression and Fusion (RNA)
- Copy number and Mutations (WGS)

Each:
- can be run in a jupyter notebook file,
- gets data from Terra workspace's gcp buckets managed by Broad's Genomics Platform + DevOps (internal only), 
- updates the sample TSVs in the processing workspace on terra with paths to the files, 
- computes the results for each samples by running workflows, 
- downloads the results and postprocesses them with additional local functions,
- performs QC and uploads them to taiga (internal only).

__ccle_tasks/__ Contains a notebook for each of the different additional processing that the CCLE team has to perform as well as one-off tasks run by the omics team

__data/__ Contains important information used for processing, including terra workspace configurations from past quarters

__depmapomics/__ Contains the core python code used in the pipeline and called by the processing jupyter notebooks

__\*\_pipeline/__ Contains some of the pipeline's workflows' wdl files and script files used by these workflows 

__temp/__ Contains the temp file that can get removed after processing (should be empty)

__documentation/__ Contains some additional files and diagrams for documenting the pipelines

__tests/__ Contains automated pytest functions used internally for development

## Pipeline Walkthrough <a name="running-pipeline"></a>

refer to [DepMap_processing_pipeline](DepMap_processing_pipeline.md)

## Auxiliary Data for the Pipeline <a name="data"></a>


### PONS

for CN pons are made from a set of ~400 normals from the GTEx project as they were sequenced in the same fashion as CCLE samples with the same set of baits. you can see the ID of the samples in `data/samples_for_pons.csv`.
We have created pons for each bait set and using XY only.
We have used workflow from the pipeline:
`gatk/CNV_Somatic_Panel_Workflow`


### targets

The data we are presenting comes from different WES targets/baits/intervals.

We are currently using Illumina ICE intervals and Agilent intervals. you can find their related PON files and interval files as parameters in our workspace files in `data/xQx.json`

__additional auxilliary data is used and listed in some of our workflow, like the CGA pipeline. You will be able to find them by looking at the wdl scripts of each pipelines and looking into the data/xQx.json` files for workspace data files.__


@[jkobject](https://www.jkobject.com)
@gkugener
@gmiller
@5im1z
@__[BroadInsitute](https://www.broadinstitute.org)

If you have any feedback or run into any issues, feel free to post an issue on the github repo.
