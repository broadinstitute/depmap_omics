# DepMap Omics: How to

## Introduction
*Give some introduction here ... Something like this (feel free to change/add):*

DepMap analyzes raw omics data on a quarterly basis... The data consists of RNAseq, etc ... 

*Put a figure pointing to the pipeline itself (I made some tweaks to the old figures. See slides 3 and 6 [here](https://docs.google.com/presentation/d/11xwUTIp2f16VLynk50mFXYVM8q3tJ7kSDTyplkt_QlY/edit?usp=sharing))*

For a more detailed understanding of the pipelines refer to [this presentation](https://docs.google.com/presentation/d/1i0HI31dBejTYmzI9Cp6Ij7--t6eSR2r9vcC22TxSnNI/edit#slide=id.g525fd14bef_0_116).

![](/assets/images/depmap-logo.png)

The following tools are used in our pipelines:
*this can be just a list with links converted into hyperlinks*


[https://www.ncbi.nlm.nih.gov/pubmed/23104886](https://www.ncbi.nlm.nih.gov/pubmed/23104886, https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)

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
  - [https://software.broadinstitute.org/gatk/documentation/tooldocs/4.beta.4/org_broadinstitute_hellbender_tools_walkers_mutect_Mutect2.php](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.beta.4/org_broadinstitute_hellbender_tools_walkers_mutect_Mutect2.php)
  - [https://www.nature.com/articles/nbt.2514](https://www.nature.com/articles/nbt.2514)
- __gatk cnv__:
  - [https://software.broadinstitute.org/gatk/documentation/article?id=11682](https://software.broadinstitute.org/gatk/documentation/article?id=11682)
- __strelka__:
  - [https://www.nature.com/articles/s41592-018-0051-x](https://www.nature.com/articles/s41592-018-0051-x)
  - [https://github.com/Illumina/strelka](https://github.com/Illumina/strelka)

Furthermore the pipelines make use of several software development tools. In particular some familiarity with the following tools are recommneded:


- [python](https://www.learnpython.org/)
- [R](https://www.codecademy.com/learn/learn-r)
- [jupyter](https://jupyter.org/index.html)
- [WDL](https://software.broadinstitute.org/wdl/documentation/)
- [gcp](https://cloud.google.com/sdk/docs/quickstart-macos)
- [docker](https://docs.docker.com/get-started/)
- [Terra](https://software.broadinstitute.org/firecloud/documentation/)
- [dalmatian](https://github.com/broadinstitute/dalmatian)
- [Terra and gcp](https://docs.google.com/document/d/1zTtaN-Px64f8JvgydZNdBbzBpFWyZzEpshSNxQh43Oc/edit#heading=h.dz5wh0l4bu9g)


The next sections are a detailed walkthrough to run DempMap omics pipelines on the [Terra](https://app.terra.bio/) platform. 

## Installation 

### clone the required repositories

This repo uses some important data and code from the [JKBio Library](https://www.github.com/jkobject/JKBio) and [gkugener](https://github.com/broadinstitute/gkugener) custom repositories. This repository and the other two dependecy repositories should be cloned into the same path using the following command:

```bash
git clone https://github.com/broadinstitute/ccle_processing.git
git clone https://github.com/jkobject/JKBio.git
git clone https://github.com/broadinstitute/gkugener.git
```


### :warning: you would need the approriate R packages and python packages -- what are the appropriate de

1. You will need to install [Jupyter Notebook](https://jupyter.org/install) and Google Cloud SDK
  - install [Google Cloud SDK](https://cloud.google.com/sdk/docs/downloads-interactive).
  - authenticate your Google account by running `gcloud auth application-default login` in the terminal.
2. For R packages, a loading function contains all required ones (in [here](https://github.com/broadinstitute/gkugener/blob/master/RScripts/load_libraries_and_annotations.R)) *how would people use this?*
install the following R packages using the provided commands:
  - install `cdsomics`: `cd src/cdsomics && R CMD INSTALL . && cd -`.
  - inst `taigr`: `cd ../JKBio/taigr && R CMD INSTALL . && cd -`.
  - And `Celllinemapr`: `cd ../JKBio/cell_line_mapping-master/celllinemapr && R CMD INSTALL . && cd -`.
3. Install the Python dependencies: `pip install -r requirements.txt`.


### Creating your Terra Workspaces:

use the `data/xQx/.json` which lists the parameters used for each workflows of each off the 3 workspaces in our pipeline (the CSV lists the workflows with their correct name):
- import the workflows, with their parameters listed in here.
- import the workspace parameters/data listed in the `GENERAL` field.

Once you have set up your workspace with the corresponding workflows, workspace data and workflow input/output parameters, you can import your data to be processed. __Tools available in `TerraFunction` in the `JKBio` package as well as the `dalmatian` package can be used to automate this process.__

Once this is done, you can run your jupyter notebook server and open one of the 3 `CCLE_\*` jupyter files corresponding to our 3 pipelines.

This notebook architecture is as follows:

1. UpLoading and preprocessing
2. Running Terra Pipelines
3. DownLoading and postProcessing
4. QC, grouping and uploading to the portal

## Running the pipeline

### 1. UpLoading and preprocessing 

The first phase really is about getting samples generated at the BroadInstitute and located into different places. Looking for duplicates and finding/adding the metadata we will need in order to have coherent and complete sample information. __This is not something that you would need to run. you can skip directly to part2__.

**Remarks:** 
- in the initialization you might want to remove any import related to `taiga` and `gsheet` to not cause any errors.
- feel free to reuse `createDatasetWithNewCellLines`, `GetNewCellLinesFromWorkspaces` or any other function for your own needs.

### 2. Running Terra Pipelines

Before running this part, you need to make sure that your dalmatian `workspacemanager` object is initialized with the right workspace you created and that the `submission`functions take as input your workflows' names. You also need to make sure that you created your sample set with all your samples and that you initialized the `sampleset` string with its name.

You can then run this part for the pipeline to run on your samples. It should take around a day.

**Remarks:**
- for the copy number pipeline we have parametrized both an XX version and an XY version, we recommend using the XY version as it covers the entire genome
- for the mutation pipeline we are working on Tumor-Normal pairs which explain some of the back and forth between the two workspace data table. (workflows works as well with or without matched normals.)
- for the expression pipeline, we have an additional set of workflows to call mutations from RNAseq, this might not be relevant to your need.

### 3. DownLoading and postProcessing (often called **2.2-4 on local** in the notebooks)

This step will do a set of tasks:
- clean some of the workspace for large useless files.
- retrieve from the workspace interesting QC results.
- copy realigned BAM files to some bucket.
- download the results.
- remove all duplicate samples from our downloaded file (keeping only thee latest version of each samples).
- saving the current pipeline configuration.

_You would only be interested here at minima in the result downloading_
 
...and post processing tasks.

> Unfortunately for now the postProcessing tasks are not all made to be easily run outside of the CCLE pipeline. Most of them are in R and are run with the Rpy2 tool.

So amongst these functions, some of them might be of a lesser interest to an external user. The most important ones for each pipelines are:

- `processSegments`
- `interpolateGapsInSegmented`
- `extendEndsOfSegments`
- `generateEntrezGenes`
- `generateGeneLevelMatrixFromSegments`

- `readMutations`
- `createSNPs`
- `filterAllelicFraction`
- `filterMinCoverage` *you would need to rewrite this function as it now runs on DepMap columns. In one line it requires that the total coverage of that site (aggregated across methods) > 8 and that there are at least 4 alternate alleles for each mutations*
- `maf_add_variant_annotations`
- `mutation_maf_to_binary_matrix`

- `readTranscripts`
- `readCounts`
- `readTPM`
- `readFusions`
- `filterFusions`
- the `step 2.2.5` where we remove samples with more then 39k transcripts with 0 readcounts.
- `prepare_depmap_\*\_for_taiga`

**Remarks:**
- in the RNAseq pipeline we have an additional sub-pipeline at the end of the notebook to process the fusion calls from starFusion.
- to get the exact same results as in CCLE, be sure to run `genecn = genecn.apply(lambda x: np.log2(1+x))` to the `genecn` dataframe in the CNV pipeline (present at the end of the validation steps).
- we do not yet have integrated our germline calling in the mutation pipeline but you can still find the haplotypeCaller\|DeepVariant workflows and their parameters.


### 4. QC, grouping and uploading to the portal

These tasks should not be very interesting for any outside user as they revolve around manual checks of the data, comparison to previous releases, etc.

We are also preparing the data to be released to different groups, removing the samples per access category: Blacklist\|Internal\|DepMapConsortium\|Public.

We are then uploading the data to a server called taiga where it will be used in the DepMap portal.
