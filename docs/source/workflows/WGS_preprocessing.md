
## WGS_preprocessing

### Inputs

#### Required

  * `dbSNP_vcf` (File, **required**)
  * `dbSNP_vcf_index` (File, **required**)
  * `input_bam` (File, **required**)
  * `input_bam_index` (File, **required**)
  * `known_indels_sites_VCFs` (Array[File], **required**)
  * `known_indels_sites_indices` (Array[File], **required**)
  * `preemptible_tries` (Int, **required**)
  * `ref_amb` (File, **required**)
  * `ref_ann` (File, **required**)
  * `ref_bwt` (File, **required**)
  * `ref_dict` (File, **required**)
  * `ref_fasta` (File, **required**)
  * `ref_fasta_index` (File, **required**)
  * `ref_name` (String, **required**)
  * `ref_pac` (File, **required**)
  * `ref_sa` (File, **required**)
  * `sample_name` (String, **required**)

#### Optional

  * `ref_alt` (File?)

#### Defaults

  * `unmapped_bam_suffix` (String, default=".bam")
  * `BamToUnmappedBams.additional_disk_size` (Int, default=20)
  * `BamToUnmappedBams.gatk_docker` (String, default="broadinstitute/gatk:latest")
  * `BamToUnmappedBams.gatk_path` (String, default="/gatk/gatk")
  * `PreProcessingForVariantDiscovery_GATK4.agg_large_disk` (Int, default=400)
  * `PreProcessingForVariantDiscovery_GATK4.agg_medium_disk` (Int, default=300)
  * `PreProcessingForVariantDiscovery_GATK4.agg_small_disk` (Int, default=200)
  * `PreProcessingForVariantDiscovery_GATK4.bwa_commandline` (String, default="bwa mem -K 100000000 -p -v 3 -t 16 -Y $bash_ref_fasta")
  * `PreProcessingForVariantDiscovery_GATK4.compression_level` (Int, default=5)
  * `PreProcessingForVariantDiscovery_GATK4.flowcell_medium_disk` (Int, default=200)
  * `PreProcessingForVariantDiscovery_GATK4.flowcell_small_disk` (Int, default=100)
  * `PreProcessingForVariantDiscovery_GATK4.gatk_docker` (String, default="us.gcr.io/broad-gatk/gatk:4.2.0.0")
  * `PreProcessingForVariantDiscovery_GATK4.gatk_path` (String, default="/gatk/gatk")
  * `PreProcessingForVariantDiscovery_GATK4.gotc_docker` (String, default="us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.7-1603303710")
  * `PreProcessingForVariantDiscovery_GATK4.gotc_path` (String, default="/usr/gitc/")
  * `PreProcessingForVariantDiscovery_GATK4.python_docker` (String, default="python:2.7")
  * `BamToUnmappedBams.RevertSam.machine_mem_gb` (Int, default=2)
  * `BamToUnmappedBams.RevertSam.preemptible_attempts` (Int, default=3)
  * `BamToUnmappedBams.SortSam.machine_mem_gb` (Int, default=4)
  * `BamToUnmappedBams.SortSam.preemptible_attempts` (Int, default=3)
  * `PreProcessingForVariantDiscovery_GATK4.ApplyBQSR.mem_size_gb` (Float, default=4)
  * `PreProcessingForVariantDiscovery_GATK4.BaseRecalibrator.mem_size_gb` (Float, default=6)
  * `PreProcessingForVariantDiscovery_GATK4.CreateSequenceGroupingTSV.mem_size_gb` (Float, default=2)
  * `PreProcessingForVariantDiscovery_GATK4.GatherBamFiles.mem_size_gb` (Float, default=3)
  * `PreProcessingForVariantDiscovery_GATK4.GatherBqsrReports.mem_size_gb` (Float, default=4)
  * `PreProcessingForVariantDiscovery_GATK4.GetBwaVersion.mem_size_gb` (Float, default=1)
  * `PreProcessingForVariantDiscovery_GATK4.MarkDuplicates.mem_size_gb` (Float, default=7.5)
  * `PreProcessingForVariantDiscovery_GATK4.MergeBamAlignment.mem_size_gb` (Float, default=4)
  * `PreProcessingForVariantDiscovery_GATK4.SamToFastqAndBwaMem.mem_size_gb` (Float, default=14)
  * `PreProcessingForVariantDiscovery_GATK4.SamToFastqAndBwaMem.num_cpu` (String, default=16)
  * `PreProcessingForVariantDiscovery_GATK4.SortAndFixTags.mem_size_gb` (Float, default=10)

### Outputs

  * `duplication_metrics` (File)
  * `bqsr_report` (File)
  * `analysis_ready_bam` (File)
  * `analysis_ready_bam_index` (File)
  * `analysis_ready_bam_md5` (File)
