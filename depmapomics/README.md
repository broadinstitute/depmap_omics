
## About the files:

> Files are organized according to the part of the pipeline or type of pipeline they belong to.

### config

The config files contains mosts parameters needed to run the pipeline. They can be overrided in the python functions.

### loading

The loading file contains all functions needed for loading new samples from Broad's genomic platform's Terra workspaces.

It uses a lot of heuristics about file possible issues the file and workspaces can contain. it tries to infer metadata and uses some google spreadsheet used by broad ops team to track samples.

### terra

The terra file contains a couple functions that are uniquely related to interacting with the Terra platform.

### fingerprinting

We do SNP fingerprinting of all bam files we receive, this file contains the necessary functions to input and get data from the SNP fingerprinting terra workspace as well as to compare it with the latest fingerprinting data we possess.

### tracker



### copy numbers



### expressions



### fusions



### mutations



### QC



### release

