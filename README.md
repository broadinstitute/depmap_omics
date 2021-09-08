# depmap_omics

__have a look at [DepMap](https://www.depmap.org)__

![](https://github.com/broadinstitute/ccle_processing/blob/editREADME/documentation/depmap-logo_white.png)

What you need to process the Quarterly DepMap-Omics releases from Terra.

[Here](https://docs.google.com/presentation/d/1i0HI31dBejTYmzI9Cp6Ij7--t6eSR2r9vcC22TxSnNI/edit#slide=id.g525fd14bef_0_116) is a presentation of the pipeline.

## Table of Contents
- [Getting Started](#quickstart)
  - [Installation](#installation)
  - [For Internal Users](#internal-users)
  - [For External Users](#external-users)
- [Running the Pipeline](#running-pipeline)
  - [Uploading and Preprocessing](#upload-preprocess)
  - [Running Terra Pipelines](#running-terra-pipelines)
  - [Downloading and Postprocessing](#downloading-postprocessing)
  - [QC, Groupding and Uploading](#qc-grouping-uploading)
- [Repository File Structure](#file-structure)
- [Auxiliary Data for the Pipeline](#data)

## Getting Started <a name="quickstart"></a>

If you are not familiar with these notions, we will first recommend you get more knowledge into each:

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

`pip install .`

### :warning: this repository needs other repos

Some important data and code from the [genepy Library](https://github.com/broadinstitute/genepy).

Use the instructions in the genepy page to install the package.

### :warning: you would need the approriate R packages and python packages

1. You will need to install jupyter notetbooks and google cloud sdk
  - install [Google Cloud SDK](https://cloud.google.com/sdk/docs/downloads-interactive).
  - authenticate my SDK account by running `gcloud auth application-default login` in the terminal.

2. and GSVA for ssGSEA in R `R` run `R -e 'if(!requireNamespace("BiocManager", quietly = TRUE)){install.packages("BiocManager")};BiocManager::install(c("GSEABase", "erccdashboard", "GSVA", "DESeq2"));`

3. For Python use the requirements.txt file `pip install -r requirements.txt` 

## For Internal Users <a name="internal-users"></a>

> To learn about the tools we use in the pipeline, see [here](#running-pipeline) for a detailed walkthrough

### Getting Terra Access

1. You will need to request access to the following terra workspaces:
  - https://app.terra.bio/#workspaces/broad-firecloud-ccle/DepMap_Mutation_Calling_CGA_pipeline
  - https://app.terra.bio/#workspaces/broad-firecloud-ccle/DepMap_hg38_RNAseq
  - https://app.terra.bio/#workspaces/broad-firecloud-ccle/DepMap_WES_CN_hg38
  - https://app.terra.bio/#workspaces/broad-firecloud-ccle/DepMap_WGS_CN
The current owners of these workspaces should give you access.
2. For the mutation pipeline you will also need to request dbGaP access (required for TCGA workflows). See CCLE/new hiree section on Asana for details.
3. Aquire access to required billing projects (e.g. broad-firecloud-ccle). See CCLE/new hiree section on Asana for details.
4. Get access to the following Terra groups:
  - DEPMAP_CCLE_DATA
  - DEPMAP-PIPELINES
  - CCLE2-DATA
  - CCLE-PIPELINE
5. If you need to get access to the depmap data, use the following links:
  - __Exome__ [IBM](https://app.terra.bio/#workspaces/terra-broad-cancer-prod/Getz_IBM_CellLines_Exomes)
  - __Exome__ [Broad](https://app.terra.bio/#workspaces/terra-broad-cancer-prod/CCLE_DepMap_WES)
  - __RNA__ [IBM](https://app.terra.bio/#workspaces/terra-broad-cancer-prod/Getz_IBM_CellLines_RNASeqData)
  - __RNA__ [Broad](https://app.terra.bio/#workspaces/terra-broad-cancer-prod/CCLE_DepMap_RNAseq)
  - __WGS__ [IBM](https://app.terra.bio/#workspaces/terra-broad-cancer-prod/Getz_IBM_CellLines_WGS)
  - __WGS__ [Broad](https://app.terra.bio/#workspaces/terra-broad-cancer-prod/DepMap_WGS)
6. Request access to the data bucket `gs://cclebams/`
7. You will need also access to the billing project `broad-firecloud-ccle`


### Additional Logins:
- In order to access and upload data, you will need to login to [taiga](https://cds.team/taiga) with your broad google account and [set up your token](https://github.com/broadinstitute/taigapy#prerequisites).
- In order to run the imports [gsheets](https://pypi.org/project/gsheets/), you need to make sure your broad google account has access to the cell line info sheets. In addition, you will need to obtain the following google API credential files:
  - Go to [console](https://console.developers.google.com/), acquire access if you don't have it already. Under credentials -> create credentials -> OAuth Client ID -> application type = Desktop app -> create. Download `client_secreat.json` and save as `~/.client_secret.json` (Refer to quickstart [here](https://gsheets.readthedocs.io/en/stable/)).
  - After calling `gsheets.from_files` for the first time using `~/.client_secret.json`, log in via google following the prompt in browser, `storage.json` will be created automatically. Save it as `~/.storage.json`.
  - Follow instruction [here](https://cloud.google.com/docs/authentication/production?authuser=1#create_service_account), create service account if needed, and save key file as `../.credentials.json`.

*Remember to share relevant gsheets with the service account (`client_email` in `../.credentials.json`).

### Boot up

- You first need to go to [taiga](https://cds.team/taiga/dataset) and create some new datasets for the virtual release

We are instantiating all the parameters needed for this pipeline to run

#### Adding new data

We are looking for new samples in a range of workspaces.

They are quite messy and might contain duplicates and/or broken file paths...

- We are thus looking at the bam files one by one and comparing them with our own bams. 
- We remove broken files, duplicates and add new version of a cell line's bam if we find some.

#### Check that we have all the cell lines we expect for this release

This involves comparing to the list in the Google sheet "Cell Line Profiling Status."

_As the list cannot be parsed, we are not comparing it for now_


## For External Users <a name="external-users"></a>

> To learn about the tools we use in the pipeline, see [here](#running-pipeline) for a detailed walkthrough

### Creating your Terra Workspaces:

1. You first need a terra account with correct billing setup
2. use the different *_pipeline folders, they contain everything to setup the workspace correctly
  1. use the WDL scripts available in:
    - https://dockstore.org/my-workflows/github.com/broadinstitute/depmap_omics/WGS_pipeline
    - https://dockstore.org/my-workflows/github.com/broadinstitute/depmap_omics/WGS_aggregate
    - https://dockstore.org/my-workflows/github.com/broadinstitute/depmap_omics/RNA_pipeline
    - https://dockstore.org/my-workflows/github.com/broadinstitute/depmap_omics/RNA_aggregate
  2. export them to your workspace
  3. add he required workspace data listed in all_configs.json (under GENERAL.workspace.attributes)
  4. setup the right inputs and outpus for your workflows using the inputs_[WORKFLOW_NAME].json and outputs_[WORKFLOW_NAME].json files
3. load your samples so that their bam and bam index google storage filepath get listed in the right data column to WGS_pipeline and RNA pipeline.
4. create a sample set with the set of samples you want to analyse.

Once this is done, you can run your jupyter notebook server and open the `CCLE_\*` jupyter files corresponding to our RNA pipeline and WGS pipeline (older versions for WES (CN and mutations are available in a previous commit labelled 20Q4)).

Remark:
  1. you cannot run the `_...()` functions listed in the notebooks.
  2. you will need to use the `postProcesssing()` functions for post processing instead of the CCLE ones.
  3. you will need to change some of the variables in the prod config file.
  4. you wont be able to run the function conditional on the CCLE boolean. You can however reimplement them to create your own pipeline.

The notebook architectures are as follows:

1. UpLoading and preprocessing
2. Running the Terra Pipelines
3. DownLoading and postProcessing the samples
4. QC, grouping and uploading to the DepMap portal

## Pipeline Walkthrough <a name="running-pipeline"></a>

![schema](https://github.com/broadinstitute/ccle_processing/blob/editREADME/documentation/architecture_diagram_white.png)

To run the CCLE pipeline we follow the installation process above and then boot up a GCP instance to run the notebooks from it.

You can find more documentation about the range of workspaces that have been created: [here](https://github.com/broadinstitute/ccle_processing/blob/master/documentation/firecloud_documentation.md)

We are using a set of key tools to process the sequencing output:
- __star__:
  - [https://www.ncbi.nlm.nih.gov/pubmed/23104886](https://www.ncbi.nlm.nih.gov/pubmed/23104886)
  - [https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)
- __rsem__: 
  - [https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-323](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-323)
- __star fusion__: 
  - [https://github.com/STAR-Fusion/STAR-Fusion/wiki](https://github.com/STAR-Fusion/STAR-Fusion/wiki)
  - [http://biorxiv.org/content/early/2017/03/24/120295](http://biorxiv.org/content/early/2017/03/24/120295)
- __mutect__: 
  - [https://software.broadinstitute.org/cancer/cga/mutect](https://software.broadinstitute.org/cancer/cga/mutect)
  - [https://youtu.be/rN-cLrb5aGs](https://youtu.be/rN-cLrb5aGs)
  - [https://gatk.broadinstitute.org/hc/en-us/articles/360036490432-Mutect2](https://gatk.broadinstitute.org/hc/en-us/articles/360036490432-Mutect2)
  - [https://www.nature.com/articles/nbt.2514](https://www.nature.com/articles/nbt.2514)
- __gatk cnv__:
  - [https://software.broadinstitute.org/gatk/documentation/article?id=11682](https://software.broadinstitute.org/gatk/documentation/article?id=11682)
- __strelka__:
  - [https://www.nature.com/articles/s41592-018-0051-x](https://www.nature.com/articles/s41592-018-0051-x)
  - [https://github.com/Illumina/strelka](https://github.com/Illumina/strelka)

The following flowchart provides another good overview of what the pipeline is doing.

![](https://github.com/broadinstitute/ccle_processing/blob/editREADME/documentation/updated-flowchart.png)

_What is explained below comes from the notebook's documentations and might be better understood by reading them directly on the notebooks_


### 1. Uploading and Preprocessing <a name="upload-preprocess"></a>

The first phase really is about getting samples generated at the broad and located into different places. Looking for duplicates and finding/adding the metadata we will need in order to have coherent and complete sample information. __This is not something that you would need to run. you can skip directly to part2__.

**Remarks:** 
- in the initialization you might want to remove any import related to `taiga` and `gsheet` to not cause any errors.
- feel free to reuse `createDatasetWithNewCellLines`, `GetNewCellLinesFromWorkspaces` or any other function for your own needs.

### 2. Running Terra Pipelines <a name="running-terra-pipelines"></a>

We are using Dalmatian to send request to Terra, so before running this part, you need to make sure that your dalmatian `WorkspaceManager` object is initialized with the right workspace you created and that the functions take as input you workflow names. You also need to make sure that you created your sample set with all your samples and that you initialized the `sampleset` string with its name
You can then run this part for the pipeline to run on your samples. It should take around a day.

#### Copy Numbers

We are running the following workflows in this order to generate the copy number dataset:

[WGS_pipeline](https://dockstore.org/workflows/github.com/broadinstitute/depmap_omics/WGS_pipeline:master?tab=info) imports and runs several sub-processes to generate copy number segments.

[WGS_aggregate](https://dockstore.org/workflows/github.com/broadinstitute/depmap_omics/WGS_aggregate:master?tab=info) aggregates CN segment files into one.

This output file for download will be saved under the sample set under the combined_seg_file attribute.

#### Mutation

We are running a set of 6 functions/workflows To generate the mutation dataset:

*   For new samples in DepMap, run the ICE version of this task. CCLE2 samples used Agilent targets, so this pipeline should be used instead. The pipelines are identical in terms of their outputs, but the proper targets, baits, and pseudo normal should be used based on how the samples were sequenced.



*NOTE:* Some further filtering occurs afterwards within the Jupyter Notebook, including hard filtering of allele frequency and coverage (see the section on local data analysis).

This outputs to be downloaded will be saved in the sample set that was run. The output we use for the release is:


*   **passedCGA_filteredMAF_aggregated** 

There are several other tasks in this workspace. In brief:

*   **CGA_Production_Analysis_Pipeline_Cell_Lines** (lelagina/CGA_Production_Analysis_Pipeline_Cell_LinesSnapshot ID: 12). This task is the same as the ICE and AGILENT prefixed version above, except that it relied on pulling the baits and targets to use from the metadata stored for the samples. Having AGILENT and ICE versions specified made the uploading and running process easier.
*   **SANGER_CGA_Production_Analysis_Pipeline_Cell_Lines** (cclf/CGA_Production_Analysis_Pipeline_Cell_Lines_debuggingSnapshot ID: 22). This task was trying to run the CGA pipeline on the Sanger WES data, using a Sanger pseudo normal. In its current implementation, this task fails to complete for the samples.
*   **UNFILTERED_aggregateMAFs_selectFields** (ccle_mg/aggregateMAFs_selectFieldsSnapshot ID: 1). Aggregates the MAF outputted by the CGA cell line pipeline prior to the common variant filter and germline filtering tasks. This can give us insight to which mutations are getting filtered out when. We may want to potentially include this MAF in the release so people can see why certain mutations of interest may be getting filtered out.
*   WES_DM_Mutation_Calling_Pipeline_(standard |expensive) (gkugener/WES_DM_Mutation_Calling_PipelineSnapshot ID: 2). This was a previous mutation calling pipeline implemented for CCLE. We do not use this pipeline any more as the CGA pipeline looks better.
*   aggregate_filterMAF_CGA (CCLE/aggregate_filterMAF_CGASnapshot ID: 1). An aggregation MAF task that we used in the past. We do not use this task anymore.
*   calculate_mutational_burden (breardon/calculate_mutational_burdenSnapshot ID: 21). This task can be used to calculate the mutational rate of the samples. We do not make use of this data in the release although it could be of interest.
*   summarizeWigFile (breardon/summarizeWigFileSnapshot ID: 5). CCLF ran this task (might be necessary for the mutational burden task). For our workflow, we do not run it.


#### RNAseq

We are generating both expression and fusion datasets with RNAseq data. We use the [GTEx pipeline](https://github.com/broadinstitute/gtex-pipeline/blob/master/TOPMed_RNAseq_pipeline.md) to generate the expression dataset, and [STAR-Fusion](https://github.com/STAR-Fusion/STAR-Fusion/wiki) to generate gene fusion calls. This task also contains a flag that lets you specify if you want to delete the intermediates (fastqs) that can be large and might cost a lot to store. Run the following tasks on all samples that you need, in this order:

[RNA_pipeline](https://dockstore.org/workflows/github.com/broadinstitute/depmap_omics/RNA_pipeline:master?tab=info) imports and runs several sub-processes to generate RNA expression and fusion data matrices.

[RNA_aggregate](https://dockstore.org/workflows/github.com/broadinstitute/depmap_omics/RNA_aggregate:master?tab=info) aggregates expression and fusion data files into their respective aggregated file.

The outputs to be downloaded will be saved under the sample set that you ran. The outputs we use for the release are:

*   rsem_genes_expected_count
*   rsem_genes_tpm
*   rsem_transcripts_tpm
*   fusions_star

### On local

__We then save the workflow configurations used__
__,delete unmapped bams generated during the process__
__and move the hg38 aligned bams to our own datastorage bucket__
__We download and reprocess removing the appended version and keeping only the newest versions__


**Remarks:**
- for the copy number pipeline we have parametrized both an XX version and an XY version, we recommend using the XY version as it covers the entire genome
- for the mutation pipeline we are working on Tumor-Normal pairs which explain some of the back and forth between the two workspace data table. (workflows works as well with or without matched normals.)
- for the expression pipeline, we have an additional set of workflows to call mutations from RNAseq, this might not be relevant to your need.

### 3. Downloading and Postprocessing (often called **2.2-4 on local** in the notebooks) <a name="downloading-postprocessing"></a>

This step will do a set of tasks:
- clean some of the workspace for large useless files.
- retrieve from the workspace interesting QC results.
- copy realigned bam files to some bucket.
- download the results.
- remove all duplicate samples from our downloaded file (keeping only thee latest version of each samples).
- saving the current pipeline configuration.

_You would only be interested here at minima in the result downloading_
 
...and post processing tasks. The main postprocessing functions (and their wrappers for internal use) can be found in the following files:

`copynumbers.py`: `postProcess()` is responsible for CN postprocessing (segmentes and gene-level CN files)

`expressions.py`: `postProcess()` is responsible for postprocessing expression data 

`fusions.py`: `postProcess()` is responsible for postprocessing aggregated fusion files

`mutations.py`: `postProcess()` is responsible for postprocessing aggregated MAF files

Please refer to the files and functions above for more details.

The post processing happens in R using guillaume's and Allie's functions, in brief:

#### CN

- processSegments
- filterForCCLE
- interpolateGapsInSegmented
- extendEndsOfSegments

##### creating gene copy number

- generateEntrezGenes
- generateGeneLevelMatrixFromSegments

##### Validation step

Once the files are saved, we load them back in python and do some validations, in brief:

- mean,max,var...
- to previous version: same mean,max,var...
- checkAmountOfSegments: flag any samples with a very high number of segments
- checkGeneChangeAccrossAll: flag any genes which stay at a similar value across all samples

#### Mutation

Here, rather than rerunning the entire analysis, because we know we are adding only WES samples, we can download the previous release's MAF, add the samples, update any annotations, and perform any global filters at the end.

First we need to do an additional step of filtering on coverage and number 

- readMutations
- createSNPs
- addToMainMutation
- filterAllelicFraction
- filterMinCoverage
- mergeAnnotations
- addAnnotation
- maf_add_variant_annotations
- mutation_maf_to_binary_matrix (x3)

##### validation

__Compare to previous release (broad only)__

I would run some checks here comparing the results to the previous releases MAF. Namely:

- Count the total number of mutations per cell line, split by type (SNP, INS, DEL)
- Count the total number of mutations observed by position (group by chromosome, start position, end position and count the number of mutations)
- Look at specific differences between the two MAFs (join on DepMap_ID, Chromosome, Start position, End position, Variant_Type). I would do this for WES only

__check if important mutations are present or not__

#### RNA

Here we get all data and remove the duplicates directly with the function `removeDuplicates`

We then run:

- readTranscripts
- readCounts
- readTPM
- renameFunction

- Allie's gene renaming / filtering and log transform

##### Validation

- mean,max,var...
- to previous version: same mean,max,var...
- we QC on the amount of genes with 0 counts for each samples

##### Generate filtered fusion table

Release: `r release`

We want to apply filters to the fusion table to reduce the number of artifacts in the dataset. Specifically, we filter the following:

* Remove fusions involving mitochondrial chromosomes, or HLA genes, or immunoglobulin genes
* Remove red herring fusions (from STAR-Fusion annotations column)
* Remove recurrent in CCLE (>= 25 samples)
* Remove fusion with (SpliceType=" INCL_NON_REF_SPLICE" and LargeAnchorSupport="No" and FFPM < 0.1)
* Remove fusions with FFPM < 0.05 (STAR-Fusion suggests using 0.1, but looking at the translocation data, this looks like it might be too aggressive)

With the functions:
- readFusions
- filterFusions
- prepare_depmap_fusion_data_for_taiga


**Remarks:**
- in the RNAseq pipeline we have an additional sub-pipeline at the end of the notebook to process the fusion calls from starFusion
- to get the exact same results as in CCLE, be sure to run `genecn = genecn.apply(lambda x: np.log2(1+x))` to the genecn dataframe in the CNV pipeline (present at the end of the validation steps).
- we do not yet have integrated our germline calling in the mutation pipeline but you can still find the haplotypeCaller\|DeepVariant workflows and their parameters


### 4. QC, Grouping and Uploading to the Portal <a name="qc-grouping-uploading"></a>

These tasks should not be very interesting for any outside user as they revolve around manual checks of the data, comparison to previous releases, etc.

We are also preparing the data to be released to different groups, removing the samples per access category: Blacklist\|Internal\|DepMapConsortium\|Public.

We are then uploading the data to a server called taiga where it will be used in the depmap portal


## Repository File Structure <a name="file-structure"></a>

There is for now 3 computation pipeline for depmap omics:
- Expression
- Mutations
- Copy number

Each:
- is contained in an jupyter notebook file
- gets data from Terra workspace's gcp buckets managed by Broad's Genomics Platform + DevOps, 
- updates the sample TSVs on Terra with path to the files, 
- compute the results for each samples by running workflows, 
- download the results, post process them with additional local functions and QC them.
- uploads them to taiga.

__data/__ Contains important information used for processing

__src__ Contains the location of function files

__\*\_pipeline__ Contains some of the pipeline's workflows' wdl files and script files used by these workflows 

__ccle_tasks__ Contains a notebook for each of the different additional processing that the CCLE team has to perform

__legacy__ Contains the previous R markdown files that were used as templates for the previous pipeline's post-processing

__readmes__ Contains some of the depmap readmes 

__temp__ Contains the temp file that can get removed after processing (should be empty)

__documentation__ Contains some additional files for documenting the pipelines


## Auxiliary Data for the Pipeline <a name="data"></a>


### PONS

for CN pons are made from a set of ~400 normals from the GTEx project as they were sequenced in the same fashion as CCLE samples with the same set of baits. you can see the ID of the samples in `data/samples_for_pons.csv`.
We have created pons for each bait set and using XY only.
We have used workflow from the pipeline:
`gatk/CNV_Somatic_Panel_Workflow`


### targets

The data we are presenting comes from different WES targets/baits/intervals.

We are currently using Illumina ICE intervals and Agilent intervals. you can find their related PON files and interval files as parameters in our workspace files in `data/xQx.json`




@[jkobject](https://www.jkobject.com)
@gkugener
@gmiller
@__[BroadInsitute](https://www.broadinstitute.org)

For any questions/issues, send us an email at ccle-help@broadinstitute.org or write down an issue on the github repo
