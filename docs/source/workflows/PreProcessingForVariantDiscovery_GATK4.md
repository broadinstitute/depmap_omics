Copyright Broad Institute, 2017

This WDL pipeline implements data pre-processing according to the GATK Best Practices 
(June 2016). Example JSONs are provided for the WGS use case but the workflow can be 
applied to Exomes and Targeted Panels.

Requirements/expectations :
- Pair-end sequencing data in unmapped BAM (uBAM) format
- One or more read groups, one per uBAM file, all belonging to a single sample (SM)
- Input uBAM files must additionally comply with the following requirements:
- - filenames all have the same suffix (we use ".unmapped.bam")
- - files must pass validation by ValidateSamFile 
- - reads are provided in query-sorted order
- - all reads must have an RG tag

Output :
- A clean BAM file and its index, suitable for variant discovery analyses.

Software version requirements (see recommended dockers in inputs JSON)
- GATK 4.beta.3 or later
- Picard 2.x
- Samtools (see gotc docker)
- Python 2.7

Cromwell version support 
- Successfully tested on v28
- Does not work on versions < v23 due to output syntax

Runtime parameters are optimized for Broad's Google Cloud Platform implementation.

LICENSING : 
This script is released under the WDL source code license (BSD-3) (see LICENSE in 
https://github.com/broadinstitute/wdl). Note however that the programs it calls may 
be subject to different licenses. Users are responsible for checking that they are
authorized to run all programs before running this script. Please see the dockers
for detailed licensing information pertaining to the included programs.

## PreProcessingForVariantDiscovery_GATK4

### Inputs

#### Required

  * `agg_large_disk` (Int, **required**)
  * `agg_medium_disk` (Int, **required**)
  * `agg_preemptible_tries` (Int, **required**)
  * `agg_small_disk` (Int, **required**)
  * `bwa_commandline` (String, **required**)
  * `compression_level` (Int, **required**)
  * `dbSNP_vcf` (File, **required**)
  * `dbSNP_vcf_index` (File, **required**)
  * `flowcell_medium_disk` (Int, **required**)
  * `flowcell_small_disk` (Int, **required**)
  * `flowcell_unmapped_bams_list` (File, **required**)
  * `gatk_docker` (String, **required**)
  * `gatk_launch_path` (String, **required**)
  * `gotc_docker` (String, **required**)
  * `gotc_path` (String, **required**)
  * `known_indels_sites_VCFs` (Array[File], **required**)
  * `known_indels_sites_indices` (Array[File], **required**)
  * `picard_docker` (String, **required**)
  * `picard_path` (String, **required**)
  * `preemptible_tries` (Int, **required**)
  * `python_docker` (String, **required**)
  * `ref_dict` (File, **required**)
  * `ref_fasta` (File, **required**)
  * `ref_fasta_index` (File, **required**)
  * `ref_name` (String, **required**)
  * `sample_name` (String, **required**)
  * `unmapped_bam_suffix` (String, **required**)
  * `ApplyBQSR.java_opt` (String, **required**)
  * `ApplyBQSR.mem_size` (String, **required**)
  * `BaseRecalibrator.java_opt` (String, **required**)
  * `BaseRecalibrator.mem_size` (String, **required**)
  * `CreateSequenceGroupingTSV.mem_size` (String, **required**)
  * `GatherBamFiles.java_opt` (String, **required**)
  * `GatherBamFiles.mem_size` (String, **required**)
  * `GatherBqsrReports.java_opt` (String, **required**)
  * `GatherBqsrReports.mem_size` (String, **required**)
  * `GetBwaVersion.mem_size` (String, **required**)
  * `MarkDuplicates.java_opt` (String, **required**)
  * `MarkDuplicates.mem_size` (String, **required**)
  * `MergeBamAlignment.java_opt` (String, **required**)
  * `MergeBamAlignment.mem_size` (String, **required**)
  * `SamToFastqAndBwaMem.java_opt` (String, **required**)
  * `SamToFastqAndBwaMem.mem_size` (String, **required**)
  * `SamToFastqAndBwaMem.num_cpu` (String, **required**)
  * `SamToFastqAndBwaMem.ref_amb` (File, **required**)
  * `SamToFastqAndBwaMem.ref_ann` (File, **required**)
  * `SamToFastqAndBwaMem.ref_bwt` (File, **required**)
  * `SamToFastqAndBwaMem.ref_pac` (File, **required**)
  * `SamToFastqAndBwaMem.ref_sa` (File, **required**)
  * `SortAndFixTags.java_opt_fix` (String, **required**)
  * `SortAndFixTags.java_opt_sort` (String, **required**)
  * `SortAndFixTags.mem_size` (String, **required**)

#### Optional

  * `SamToFastqAndBwaMem.ref_alt` (File?)

#### Defaults

  * `base_file_name` (String, default=sample_name + "." + ref_name)
  * `flowcell_unmapped_bams` (Array[File], default=read_lines(flowcell_unmapped_bams_list))

### Outputs

  * `duplication_metrics` (File)
  * `bqsr_report` (File)
  * `analysis_ready_bam` (File)
  * `analysis_ready_bam_index` (File)
  * `analysis_ready_bam_md5` (File)

## GetBwaVersion

### Inputs

#### Required

  * `bwa_path` (String, **required**)
  * `docker_image` (String, **required**)
  * `mem_size` (String, **required**)
  * `preemptible_tries` (Int, **required**)

### Outputs

  * `version` (String)

## SamToFastqAndBwaMem

### Inputs

#### Required

  * `bwa_commandline` (String, **required**)
  * `bwa_path` (String, **required**)
  * `compression_level` (Int, **required**)
  * `disk_size` (Int, **required**)
  * `docker_image` (String, **required**)
  * `input_bam` (File, **required**)
  * `java_opt` (String, **required**)
  * `mem_size` (String, **required**)
  * `num_cpu` (String, **required**)
  * `output_bam_basename` (String, **required**)
  * `picard_path` (String, **required**)
  * `preemptible_tries` (Int, **required**)
  * `ref_amb` (File, **required**)
  * `ref_ann` (File, **required**)
  * `ref_bwt` (File, **required**)
  * `ref_dict` (File, **required**)
  * `ref_fasta` (File, **required**)
  * `ref_fasta_index` (File, **required**)
  * `ref_pac` (File, **required**)
  * `ref_sa` (File, **required**)

#### Optional

  * `ref_alt` (File?)

### Outputs

  * `output_bam` (File)
  * `bwa_stderr_log` (File)

## MergeBamAlignment

### Inputs

#### Required

  * `aligned_bam` (File, **required**)
  * `bwa_commandline` (String, **required**)
  * `bwa_version` (String, **required**)
  * `compression_level` (Int, **required**)
  * `disk_size` (Int, **required**)
  * `docker_image` (String, **required**)
  * `java_opt` (String, **required**)
  * `mem_size` (String, **required**)
  * `output_bam_basename` (String, **required**)
  * `picard_path` (String, **required**)
  * `preemptible_tries` (Int, **required**)
  * `ref_dict` (File, **required**)
  * `ref_fasta` (File, **required**)
  * `ref_fasta_index` (File, **required**)
  * `unmapped_bam` (File, **required**)

### Outputs

  * `output_bam` (File)

## SortAndFixTags

### Inputs

#### Required

  * `compression_level` (Int, **required**)
  * `disk_size` (Int, **required**)
  * `docker_image` (String, **required**)
  * `input_bam` (File, **required**)
  * `java_opt_fix` (String, **required**)
  * `java_opt_sort` (String, **required**)
  * `mem_size` (String, **required**)
  * `output_bam_basename` (String, **required**)
  * `picard_path` (String, **required**)
  * `preemptible_tries` (Int, **required**)
  * `ref_dict` (File, **required**)
  * `ref_fasta` (File, **required**)
  * `ref_fasta_index` (File, **required**)

### Outputs

  * `output_bam` (File)
  * `output_bam_index` (File)
  * `output_bam_md5` (File)

## MarkDuplicates

### Inputs

#### Required

  * `compression_level` (Int, **required**)
  * `disk_size` (Int, **required**)
  * `docker_image` (String, **required**)
  * `input_bams` (Array[File], **required**)
  * `java_opt` (String, **required**)
  * `mem_size` (String, **required**)
  * `metrics_filename` (String, **required**)
  * `output_bam_basename` (String, **required**)
  * `picard_path` (String, **required**)
  * `preemptible_tries` (Int, **required**)

### Outputs

  * `output_bam` (File)
  * `duplicate_metrics` (File)

## CreateSequenceGroupingTSV

### Inputs

#### Required

  * `docker_image` (String, **required**)
  * `mem_size` (String, **required**)
  * `preemptible_tries` (Int, **required**)
  * `ref_dict` (File, **required**)

### Outputs

  * `sequence_grouping` (Array[Array[String]])
  * `sequence_grouping_with_unmapped` (Array[Array[String]])

## BaseRecalibrator

### Inputs

#### Required

  * `dbSNP_vcf` (File, **required**)
  * `dbSNP_vcf_index` (File, **required**)
  * `disk_size` (Int, **required**)
  * `docker_image` (String, **required**)
  * `gatk_launch_path` (String, **required**)
  * `input_bam` (File, **required**)
  * `input_bam_index` (File, **required**)
  * `java_opt` (String, **required**)
  * `known_indels_sites_VCFs` (Array[File], **required**)
  * `known_indels_sites_indices` (Array[File], **required**)
  * `mem_size` (String, **required**)
  * `preemptible_tries` (Int, **required**)
  * `recalibration_report_filename` (String, **required**)
  * `ref_dict` (File, **required**)
  * `ref_fasta` (File, **required**)
  * `ref_fasta_index` (File, **required**)
  * `sequence_group_interval` (Array[String], **required**)

### Outputs

  * `recalibration_report` (File)

## GatherBqsrReports

### Inputs

#### Required

  * `disk_size` (Int, **required**)
  * `docker_image` (String, **required**)
  * `gatk_launch_path` (String, **required**)
  * `input_bqsr_reports` (Array[File], **required**)
  * `java_opt` (String, **required**)
  * `mem_size` (String, **required**)
  * `output_report_filename` (String, **required**)
  * `preemptible_tries` (Int, **required**)

### Outputs

  * `output_bqsr_report` (File)

## ApplyBQSR

### Inputs

#### Required

  * `disk_size` (Int, **required**)
  * `docker_image` (String, **required**)
  * `gatk_launch_path` (String, **required**)
  * `input_bam` (File, **required**)
  * `input_bam_index` (File, **required**)
  * `java_opt` (String, **required**)
  * `mem_size` (String, **required**)
  * `output_bam_basename` (String, **required**)
  * `preemptible_tries` (Int, **required**)
  * `recalibration_report` (File, **required**)
  * `ref_dict` (File, **required**)
  * `ref_fasta` (File, **required**)
  * `ref_fasta_index` (File, **required**)
  * `sequence_group_interval` (Array[String], **required**)

### Outputs

  * `recalibrated_bam` (File)

## GatherBamFiles

### Inputs

#### Required

  * `compression_level` (Int, **required**)
  * `disk_size` (Int, **required**)
  * `docker_image` (String, **required**)
  * `input_bams` (Array[File], **required**)
  * `java_opt` (String, **required**)
  * `mem_size` (String, **required**)
  * `output_bam_basename` (String, **required**)
  * `picard_path` (String, **required**)
  * `preemptible_tries` (Int, **required**)

### Outputs

  * `output_bam` (File)
  * `output_bam_index` (File)
  * `output_bam_md5` (File)
