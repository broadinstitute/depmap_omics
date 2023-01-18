Copyright Broad Institute, 2017

This WDL reverts a SAM or BAM file to uBAMs, one per readgroup 

Requirements/expectations :
- Pair-end sequencing data in SAM or BAM format
- One or more read groups

Outputs :
- Set of unmapped BAMs, one per read group, with reads sorted by queryname

Cromwell version support 
- Successfully tested on v27
- Does not work on versions < v23 due to output syntax

Runtime parameters are optimized for Broad's Google Cloud Platform implementation. 
For program versions, see docker containers. 

LICENSING : 
This script is released under the WDL source code license (BSD-3) (see LICENSE in 
https://github.com/broadinstitute/wdl). Note however that the programs it calls may 
be subject to different licenses. Users are responsible for checking that they are
authorized to run all programs before running this script. Please see the docker 
page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
licensing information pertaining to the included programs.

## BamToUnmappedRGBamsWf

### Inputs

#### Required

  * `input_bam` (File, **required**)
  * `picard_docker` (String, **required**)
  * `picard_path` (String, **required**)
  * `preemptible_tries` (Int, **required**)
  * `ref_fasta` (File, **required**)
  * `ref_fasta_index` (File, **required**)
  * `RevertBamToUnmappedRGBams.disk_size` (Int, **required**)
  * `RevertBamToUnmappedRGBams.java_opt` (String, **required**)
  * `RevertBamToUnmappedRGBams.mem_size` (String, **required**)
  * `RevertBamToUnmappedRGBams.output_dir` (String, **required**)
  * `SortBamByQueryname.disk_size` (Int, **required**)
  * `SortBamByQueryname.java_opt` (String, **required**)
  * `SortBamByQueryname.mem_size` (String, **required**)
  * `ValidateSamFile.disk_size` (Int, **required**)
  * `ValidateSamFile.java_opt` (String, **required**)
  * `ValidateSamFile.mem_size` (String, **required**)

### Outputs

  * `sortsam_out` (Array[File])
  * `validatesam_out` (Array[File])

## RevertBamToUnmappedRGBams

### Inputs

#### Required

  * `disk_size` (Int, **required**)
  * `docker_image` (String, **required**)
  * `input_bam` (File, **required**)
  * `java_opt` (String, **required**)
  * `mem_size` (String, **required**)
  * `output_dir` (String, **required**)
  * `picard_path` (String, **required**)

### Outputs

  * `unmapped_bams` (Array[File])

## SortBamByQueryname

### Inputs

#### Required

  * `disk_size` (Int, **required**)
  * `docker_image` (String, **required**)
  * `input_bam` (File, **required**)
  * `java_opt` (String, **required**)
  * `mem_size` (String, **required**)
  * `picard_path` (String, **required**)
  * `preemptible_tries` (Int, **required**)
  * `sorted_bam_name` (String, **required**)

### Outputs

  * `sorted_bam` (File)

## ValidateSamFile

### Inputs

#### Required

  * `disk_size` (Int, **required**)
  * `docker_image` (String, **required**)
  * `input_bam` (File, **required**)
  * `java_opt` (String, **required**)
  * `mem_size` (String, **required**)
  * `picard_path` (String, **required**)
  * `preemptible_tries` (Int, **required**)
  * `report_filename` (String, **required**)

### Outputs

  * `report` (File)
