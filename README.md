# ccle_processing
What you need to process the Quarterly DepMap-Omics releases from Terra
Here is a presentation of the pipeline: https://docs.google.com/presentation/d/1i0HI31dBejTYmzI9Cp6Ij7--t6eSR2r9vcC22TxSnNI/edit#slide=id.g525fd14bef_0_116

## Instalation

if you are not familiar with these notions, we will first recommand you get more knowledge into each:
- python https://www.learnpython.org/
- R https://www.codecademy.com/learn/learn-r
- jupyter https://jupyter.org/index.html
- WDL https://software.broadinstitute.org/wdl/documentation/
- gcp https://cloud.google.com/sdk/docs/quickstart-macos
- docker https://docs.docker.com/get-started/
- Terra https://software.broadinstitute.org/firecloud/documentation/
- dalmatian https://github.com/broadinstitute/dalmatian
- Terra and gcp https://docs.google.com/document/d/1zTtaN-Px64f8JvgydZNdBbzBpFWyZzEpshSNxQh43Oc/edit#heading=h.dz5wh0l4bu9g

- star:
  - https://www.ncbi.nlm.nih.gov/pubmed/23104886
  - https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf
- rsem: 
  - https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-323
- star fusion: 
  - https://github.com/STAR-Fusion/STAR-Fusion/wiki
  - http://biorxiv.org/content/early/2017/03/24/120295
- mutect: 
  - https://software.broadinstitute.org/cancer/cga/mutect
  - https://software.broadinstitute.org/gatk/documentation/tooldocs/4.beta.4/org_broadinstitute_hellbender_tools_walkers_mutect_Mutect2.php
  - https://www.nature.com/articles/nbt.2514
- gatk cnv 
  - https://software.broadinstitute.org/gatk/documentation/article?id=11682
- strelka:
  - https://www.nature.com/articles/s41592-018-0051-x
  - https://github.com/Illumina/strelka




### /!\ this repository needs other repos 
some important data and code from the [JKBio Library](https://www.github.com/jkobject/JKBio) and [gkugener](https://github.com/broadinstitute/gkugener)
Go to the repos and pull them to the same parent folder as ccle_processing.

### /!\ you would need the approriate R packages and python packages:
- for R packages, a loading function contains all required ones (in [here](https://github.com/broadinstitute/gkugener/blob/master/RScripts/load_libraries_and_annotations.R))
- for Python use the requirements.txt file `pip install -r requirements.txt`


### Get Terra Access

1. you will need to request access to the following terra workspaces:
  - https://app.terra.bio/#workspaces/broad-firecloud-ccle/DepMap_Mutation_Calling_CGA_pipeline
  - https://app.terra.bio/#workspaces/broad-firecloud-ccle/DepMap_hg38_RNAseq
  - https://app.terra.bio/#workspaces/broad-firecloud-ccle/DepMap_WES_CN_hg38
2. For the mutation pipeline you will also need to request DBGap access (required for TCGA workflows): https://docs.google.com/document/d/1hwsjUypqUpse7IeMyBLKEXmdlXUzfBa4e4p9teCVtaw/edit?usp=sharing
3. Ask Sarah Young for access to required billing projects (e.g. broad-firecloud-ccle)
4. get access to the following Terra groups:
  - DEPMAP_CCLE_DATA
  - DEPMAP-PIPELINES
  - CCLE2-DATA
  - CCLE-PIPELINE
5. if you need to get access to the depmap data:
  - __Exome__ [IBM](https://firecloud.terra.bio/#workspaces/broad-genomics-delivery/Getz_IBM_CellLines_Exomes/data)
  - __Exome__ [Broad](https://firecloud.terra.bio/#workspaces/broad-firecloud-ccle/CCLE_DepMap_WES)
  - __Exome__ [Broad](https://firecloud.terra.bio/#workspaces/broad-genomics-delivery/CCLE_DepMap_WES)
  - __rna__ [IBM](https://firecloud.terra.bio/#workspaces/broad-genomics-delivery/Getz_IBM_CellLines_RNASeqData/data)
  - __rna__ [Broad](https://firecloud.terra.bio/#workspaces/broad-firecloud-ccle/CCLE_DepMap_RNAseq)
  - __rna__ [Broad](https://firecloud.terra.bio/#workspaces/broad-genomics-delivery/CCLE_DepMap_RNAseq)
6. request access to the metadata bucket `gs://ccle_default_params/`

*more information on the firecloud workspaces and what they might contain is available on this [document](https://www.github.com/broadinstitute/ccle_processing/firecloud_documentation.html)

## File structure

there is for now 3 computation pipeline for depmap omics:
- expression
- mutations
- copy number

each:
- is contained in an jupyter notebook file
- gets data from buckets, 
- updates the TSVs on Terra, 
- compute the results for each, 
- QC them and do some more filtering,
- Uploads them to taiga.

__data/__ contains important information used for processing

__src__ contains the location of function files

__\*\_pipeline__ contains some of the pipeline wdl files and script files 

__ccle_tasks__ contains a notebook for each of the different additional processing that the CCLE team has to perform

__legacy__ contains the previous R markdown files that were used as templates for the previous pipeline's postprocessing

__readmes__ contains the readmes that are updated at each depmap releases 

__temp__ contains the temp file that can get removed after processing (should be empty)
