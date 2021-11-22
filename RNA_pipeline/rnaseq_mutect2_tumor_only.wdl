version 1.0
## Copyright Broad Institute, 2018
##
## Workflows for processing RNA data for germline short variant discovery with GATK (v3+v4) and related tools
##
## Requirements/expectations :
## - BAM
##
## Output :
##
## Runtime parameters are optimized for Broad's Google Cloud Platform implementation.
## For program versions, see docker containers.
##
## LICENSING :
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.

struct Runtime {
  String gatk_docker
  File? gatk_override
  Int max_retries
  Int preemptible
  Int cpu
  Int machine_mem
  Int command_mem
  Int disk
  Int boot_disk_size
}

workflow RNAseq_mutect2 {
  input {
    # Mutect2 inputs
    File intervals
    File ref_fasta
    File ref_fai
    File ref_dict
    File tumor_bam
    File tumor_bai
    File? normal_reads
    File? normal_reads_index
    File? pon
    File? pon_idx
    Int scatter_count
    File? gnomad
    File? gnomad_idx
    File? variants_for_contamination
    File? variants_for_contamination_idx
    File? realignment_index_bundle
    String? realignment_extra_args
    String? m2_extra_args
    String? m2_extra_filtering_args
    String? getpileupsummaries_extra_args
    String? split_intervals_extra_args
    Boolean? make_bamout
    Boolean? compress_vcfs
    File? gga_vcf
    File? gga_vcf_idx
    String? gcs_project_for_requester_pays

    # Funcotator inputs
    Boolean? run_funcotator
    String? sequencing_center
    String? sequence_source
    String? funco_reference_version
    String? funco_output_format
    Boolean? funco_compress
    Boolean? funco_use_gnomad_AF
    File? funco_data_sources_tar_gz
    String? funco_transcript_selection_mode
    File? funco_transcript_selection_list
    Array[String]? funco_annotation_defaults
    Array[String]? funco_annotation_overrides
    Array[String]? funcotator_excluded_fields
    Boolean? funco_filter_funcotations
    String? funcotator_extra_args

    String funco_default_output_format = "MAF"

    # runtime
    String gatk_docker
    File? gatk_override
    String basic_bash_docker = "ubuntu:16.04"
    Boolean? filter_funcotations

    Int? preemptible
    Int? max_retries
    Int small_task_cpu = 2
    Int small_task_mem = 4
    Int small_task_disk = 100
    Int boot_disk_size = 12
    Int learn_read_orientation_mem = 8000
    Int filter_alignment_artifacts_mem = 9000

    # Use as a last resort to increase the disk given to every task in case of ill behaving data
    Int? emergency_extra_disk

    # These are multipliers to multipler inputs by to make sure we have enough disk to accommodate for possible output sizes
    # Large is for Bams/WGS vcfs
    # Small is for metrics/other vcfs
    Float large_input_to_output_multiplier = 2.25
    Float small_input_to_output_multiplier = 2.0
    Float cram_to_bam_multiplier = 6.0

    String? gitc_docker_override

    Array[File] knownVcfs
    Array[File] knownVcfsIndices

    File dbSnpVcf
    File dbSnpVcfIndex

    String? SO
    String? RGLB
    String? RGPL
    String? RGPU
    String? RGSM
    String? VALIDATION_STRINGENCY

    String? gatk_path_override
  }

  Int preemptible_or_default = select_first([preemptible, 4])
  Int max_retries_or_default = select_first([max_retries, 2])

  Boolean compress = select_first([compress_vcfs, false])
  Boolean make_bamout_or_default = select_first([make_bamout, false])
  Boolean run_funcotator_or_default = select_first([run_funcotator, false])
  Boolean filter_funcotations_or_default = select_first([filter_funcotations, true])

  # Disk sizes used for dynamic sizing
  Int ref_size = ceil(size(ref_fasta, "GB") + size(ref_dict, "GB") + size(ref_fai, "GB"))
  Int tumor_reads_size = ceil(size(tumor_bam, "GB") + size(tumor_bai, "GB"))
  Int gnomad_vcf_size = if defined(gnomad) then ceil(size(gnomad, "GB")) else 0
  Int normal_reads_size = if defined(normal_reads) then ceil(size(normal_reads, "GB") + size(normal_reads_index, "GB")) else 0

  # If no tar is provided, the task downloads one from broads ftp server
  Int funco_tar_size = if defined(funco_data_sources_tar_gz) then ceil(size(funco_data_sources_tar_gz, "GB") * 3) else 100
  Int gatk_override_size = if defined(gatk_override) then ceil(size(gatk_override, "GB")) else 0

  # This is added to every task as padding, should increase if systematically you need more disk for every call
  Int disk_pad = 10 + gatk_override_size + select_first([emergency_extra_disk,0])

  # logic about output file names -- these are the names *without* .vcf extensions
  String output_basename = basename(tumor_bam, ".bam")
  String unfiltered_name = output_basename + "-unfiltered"
  String filtered_name = output_basename + "-filtered"
  String funcotated_name = output_basename + "-funcotated"

  String output_vcf_name = output_basename + ".vcf"

  Int tumor_cram_to_bam_disk = ceil(tumor_reads_size * cram_to_bam_multiplier)
  Int normal_cram_to_bam_disk = ceil(normal_reads_size * cram_to_bam_multiplier)

  String gitc_docker = select_first([gitc_docker_override, "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"])

  ## Optional user optimizations

  Runtime standard_runtime = {"gatk_docker": gatk_docker, "gatk_override": gatk_override,
    "max_retries": max_retries_or_default, "preemptible": preemptible_or_default, "cpu": small_task_cpu,
    "machine_mem": small_task_mem * 1000, "command_mem": small_task_mem * 1000 - 500,
    "disk": small_task_disk + disk_pad, "boot_disk_size": boot_disk_size}

  #call ReorderSam_GATK #optional
  Int tumor_bam_size = ceil(size(tumor_bam, "GB") + size(tumor_bai, "GB"))
  Int disksize = (tumor_bam_size*12) + ref_size + disk_pad

  call picard_CleanAfterStar {
    input:
      input_bam = tumor_bam,
      input_bai = tumor_bai,
      ref_fasta = ref_fasta,
      ref_dict = ref_dict,
      ref_fasta_index = ref_fai,
      SO = SO,
      base_name = output_basename,
      RGLB = RGLB,
      RGPL = RGPL,
      RGPU = RGPU,
      RGSM = RGSM,
      VALIDATION_STRINGENCY = VALIDATION_STRINGENCY,
      preemptible_count = preemptible_or_default,
      docker = gatk_docker,
      disk = disksize,
  }

  call SplitNCigarReads_GATK4 {
    input:
      input_bam = picard_CleanAfterStar.output_bam,
      input_bam_index = picard_CleanAfterStar.output_bam_index,
      base_name = output_basename + ".split",
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fai,
      ref_dict = ref_dict,
      interval_list = intervals,
      preemptible_count = preemptible_or_default,
      docker = gatk_docker,
      disk = disksize,
  }

  call BaseRecalibrator {
    input:
      input_bam = SplitNCigarReads_GATK4.output_bam,
      input_bam_index = SplitNCigarReads_GATK4.output_bam_index,
      recal_output_file = output_basename + ".recal_data.csv",
      dbSNP_vcf = dbSnpVcf,
      dbSNP_vcf_index = dbSnpVcfIndex,
      known_indels_sites_VCFs = knownVcfs,
      known_indels_sites_indices = knownVcfsIndices,
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fai,
      preemptible_count = preemptible_or_default,
      docker = gatk_docker,
  }

  call ApplyBQSR {
    input:
      input_bam =  SplitNCigarReads_GATK4.output_bam,
      input_bam_index = SplitNCigarReads_GATK4.output_bam_index,
      base_name = output_basename + ".aligned.duplicates_marked.recalibrated",
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fai,
      ref_dict = ref_dict,
      recalibration_report = BaseRecalibrator.recalibration_report,
      preemptible_count = preemptible_or_default,
      docker = gatk_docker,
  }

  call ScatterIntervalList_GATK4 {
    input:
      interval_list = intervals,
      scatter_count = scatter_count,
      preemptible_count = preemptible_or_default,
      docker = gatk_docker,
  }

  Int m2_output_size = tumor_bam_size / scatter_count
  Int m2_per_scatter_size = tumor_bam_size + ref_size + gnomad_vcf_size + m2_output_size + disk_pad

  scatter (subintervals in ScatterIntervalList_GATK4.out) {
    call M2 {
      input:
          intervals = subintervals,
          ref_fasta = ref_fasta,
          ref_fai = ref_fai,
          ref_dict = ref_dict,
          tumor_bam = ApplyBQSR.output_bam,
          tumor_bai = ApplyBQSR.output_bam_index,
          pon = pon,
          pon_idx = pon_idx,
          gnomad = gnomad,
          gnomad_idx = gnomad_idx,
          preemptible = preemptible,
          max_retries = max_retries,
          m2_extra_args = m2_extra_args,
          getpileupsummaries_extra_args = getpileupsummaries_extra_args,
          variants_for_contamination = variants_for_contamination,
          variants_for_contamination_idx = variants_for_contamination_idx,
          make_bamout = make_bamout_or_default,
          compress = compress,
          gga_vcf = gga_vcf,
          gga_vcf_idx = gga_vcf_idx,
          gatk_override = gatk_override,
          gatk_docker = gatk_docker,
          disk_space = m2_per_scatter_size,
          gcs_project_for_requester_pays = gcs_project_for_requester_pays
    }
  }

  Int merged_vcf_size = ceil(size(M2.unfiltered_vcf, "GB"))
  Int merged_bamout_size = ceil(size(M2.output_bamOut, "GB"))

  call MergeVCFs {
    input:
      input_vcfs = M2.unfiltered_vcf,
      input_vcfs_indexes =  M2.unfiltered_vcf_idx,
      output_vcf_name = output_basename + ".g.vcf.gz",
      preemptible_count = preemptible_or_default,
      docker = gatk_docker,
  }

  call MergeStats { input: stats = M2.stats, runtime_params = standard_runtime }

  call Filter {
    input:
        ref_fasta = ref_fasta,
        ref_fai = ref_fai,
        ref_dict = ref_dict,
        intervals = intervals,
        unfiltered_vcf = MergeVCFs.output_vcf,
        unfiltered_vcf_idx = MergeVCFs.output_vcf_index,
        output_name = filtered_name,
        compress = compress,
        mutect_stats = MergeStats.merged_stats,
        m2_extra_filtering_args = m2_extra_filtering_args,
        runtime_params = standard_runtime,
        disk_space = ceil(size(MergeVCFs.output_vcf, "GB") * small_input_to_output_multiplier) + disk_pad
  }

  call VariantFiltration {
    input:
      input_vcf = Filter.filtered_vcf,
      input_vcf_index = Filter.filtered_vcf_idx,
      base_name = output_basename + ".variant_filtered.vcf.gz",
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fai,
      ref_dict = ref_dict,
      preemptible_count = preemptible_or_default,
      docker = gatk_docker,
  }

  if (run_funcotator_or_default) {
    call Funcotate {
      input:
        ref_fasta = ref_fasta,
        ref_fai = ref_fai,
        ref_dict = ref_dict,
        input_vcf = VariantFiltration.output_vcf,
        input_vcf_idx = VariantFiltration.output_vcf_index,
        reference_version = select_first([funco_reference_version, "hg19"]),
        output_file_base_name = basename(VariantFiltration.output_vcf, ".vcf") + ".annotated",
        output_format = if defined(funco_output_format) then "" + funco_output_format else funco_default_output_format,
        compress = if defined(funco_compress) then select_first([funco_compress]) else false,
        use_gnomad = if defined(funco_use_gnomad_AF) then select_first([funco_use_gnomad_AF]) else false,
        data_sources_tar_gz = funco_data_sources_tar_gz,
        case_id = M2.tumor_sample[0],
        control_id = M2.normal_sample[0],
        sequencing_center = sequencing_center,
        sequence_source = sequence_source,
        transcript_selection_mode = funco_transcript_selection_mode,
        transcript_selection_list = funco_transcript_selection_list,
        annotation_defaults = funco_annotation_defaults,
        annotation_overrides = funco_annotation_overrides,
        funcotator_excluded_fields = funcotator_excluded_fields,
        filter_funcotations = filter_funcotations_or_default,
        extra_args = funcotator_extra_args,
        runtime_params = standard_runtime,
        disk_space = ceil(size(VariantFiltration.output_vcf, "GB") * large_input_to_output_multiplier)  + funco_tar_size + disk_pad
    }
  }

  output {
    File merged_vcf = MergeVCFs.output_vcf
    File merged_vcf_index = MergeVCFs.output_vcf_index
    File variant_filtered_vcf = VariantFiltration.output_vcf
    File variant_filtered_vcf_index = VariantFiltration.output_vcf_index
    File filtering_stats = Filter.filtering_stats
    File mutect_stats = MergeStats.merged_stats

    File? funcotated_file = Funcotate.funcotated_output_file
    File? funcotated_file_index = Funcotate.funcotated_output_file_index
  }
}

task SplitNCigarReads_GATK4 {
  input {
    File input_bam
    File input_bam_index
    String base_name
    File interval_list

    File ref_fasta
    File ref_fasta_index
    File ref_dict

    String? gatk_override
    String docker
    Int? preemptible_count
    Int? memory
    Int? disk
  }

  command <<<
    export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

    gatk \
      --java-options "-Xmx6G" \
      SplitNCigarReads \
      -R ~{ref_fasta} \
      -I ~{input_bam} \
      -O ~{base_name}.bam
  >>>

  output {
    File output_bam = "~{base_name}.bam"
    File output_bam_index = "~{base_name}.bai"
  }

  runtime {
    disks: "local-disk "+select_first([disk,200])+" HDD"
    docker: docker
    memory: select_first([memory,16])+" GB"
    preemptible: select_first([preemptible_count, 10])
  }
}

task BaseRecalibrator {
  input {
    File input_bam
    File input_bam_index
    String recal_output_file

    File dbSNP_vcf
    File dbSNP_vcf_index
    Array[File] known_indels_sites_VCFs
    Array[File] known_indels_sites_indices

    File ref_dict
    File ref_fasta
    File ref_fasta_index

    String? gatk_override

    String docker
    Int? preemptible_count
    Int? memory
    Int? disk
  }

  command <<<
    export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

    gatk --java-options "-XX:GCTimeLimit=100 -XX:GCHeapFreeLimit=5 -XX:+PrintFlagsFinal \
      -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails \
      -Xloggc:gc_log.log -Xms8192m" \
      BaseRecalibrator \
      -R ~{ref_fasta} \
      -I ~{input_bam} \
      --use-original-qualities \
      -O ~{recal_output_file} \
      -known-sites ~{dbSNP_vcf} \
      -known-sites ~{sep=" --known-sites " known_indels_sites_VCFs}
  >>>

  output {
    File recalibration_report = recal_output_file
  }

  runtime {
    memory: select_first([memory,16])+" GB"
    disks: "local-disk "+ select_first([disk,200])+" HDD"
    docker: docker
    preemptible: select_first([preemptible_count, 10])
  }
}


task ApplyBQSR {
  input {
    File input_bam
    File input_bam_index
    String base_name
    File recalibration_report

    File ref_dict
    File ref_fasta
    File ref_fasta_index

    String? gatk_override

    String docker
    Int? preemptible_count
    Int? memory
  }

  command <<<
    export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

    gatk \
    --java-options "-Xms8196m" \
    ApplyBQSR \
    --add-output-sam-program-record \
    -R ~{ref_fasta} \
    -I ~{input_bam} \
    --use-original-qualities \
    -O ~{base_name}.bam \
    --bqsr-recal-file ~{recalibration_report}
  >>>

  output {
    File output_bam = "~{base_name}.bam"
    File output_bam_index = "~{base_name}.bai"
  }

  runtime {
    memory: select_first([memory,16])+" GB"
    disks: "local-disk 300 HDD"
    preemptible: select_first([preemptible_count, 10])
    docker: docker
  }
}


task MergeVCFs {
  input {
    Array[File] input_vcfs
    Array[File] input_vcfs_indexes
    String output_vcf_name

    Int? disk_size = 10
    Int? memory

    String? gatk_override

    String docker
    Int? preemptible_count
  }

  # Using MergeVcfs instead of GatherVcfs so we can create indices
  # See https://github.com/broadinstitute/picard/issues/789 for relevant GatherVcfs ticket
  command <<<
    export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

    gatk --java-options "-Xms4096m"  \
      MergeVcfs \
      --INPUT ~{sep=' --INPUT ' input_vcfs} \
      --OUTPUT ~{output_vcf_name}
  >>>

  output {
    File output_vcf = output_vcf_name
    File output_vcf_index = "~{output_vcf_name}.tbi"
  }

  runtime {
    memory: select_first([memory,6])+" GB"
    disks: "local-disk " + disk_size + " HDD"
    docker: docker
    preemptible: select_first([preemptible_count, 10])
  }
}

task ScatterIntervalList_GATK4 {
  input {
    File interval_list
    Int scatter_count
    String? gatk_override
    String docker
    Int? preemptible_count
    Int? memory
  }

  command <<<
    set -e
    export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}
    mkdir out
    gatk --java-options "-Xms2g" \
      IntervalListTools \
      --SCATTER_COUNT ~{scatter_count} \
      --SUBDIVISION_MODE BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \
      --UNIQUE true \
      --SORT true \
      --INPUT ~{interval_list} \
      --OUTPUT out
    python <<CODE
    import glob, os
    # Works around a JES limitation where multiples files with the same name overwrite each other when globbed
    intervals = sorted(glob.glob("out/*/*.interval_list"))
    for i, interval in enumerate(intervals):
      (directory, filename) = os.path.split(interval)
      newName = os.path.join(directory, str(i + 1) + filename)
      os.rename(interval, newName)
    print(len(intervals))
    f = open("interval_count.txt", "w+")
    f.write(str(len(intervals)))
    f.close()
    CODE
  >>>

  output {
    Array[File] out = glob("out/*/*.interval_list")
    Int interval_count = read_int("interval_count.txt")
  }

  runtime {
    disks: "local-disk 1 HDD"
    memory: select_first([memory,4])+" GB"
    docker: docker
    preemptible: select_first([preemptible_count, 10])
  }
}

task RevertSam {
  input {
    File input_bam
    String base_name
    String sort_order

    String? gatk_override

    String docker
    Int? preemptible_count

    Int? disk
  }

  command <<<
    export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

    gatk \
      RevertSam \
      --INPUT ~{input_bam} \
      --OUTPUT ~{base_name}.bam \
      --VALIDATION_STRINGENCY SILENT \
      --ATTRIBUTE_TO_CLEAR FT \
      --ATTRIBUTE_TO_CLEAR CO \
      --SORT_ORDER ~{sort_order}
  >>>

  output {
    File output_bam = "~{base_name}.bam"
  }

  runtime {
    docker: docker
    disks: "local-disk "+select_first([disk,200])+" HDD"
    memory: "8 GB"
    preemptible: select_first([preemptible_count, 10])
  }
}

task M2 {
  input {
    File? intervals
    File ref_fasta
    File ref_fai
    File ref_dict
    File tumor_bam
    File tumor_bai
    File? normal_bam
    File? normal_bai
    File? pon
    File? pon_idx
    File? gnomad
    File? gnomad_idx
    String? m2_extra_args
    String? getpileupsummaries_extra_args
    Boolean? make_bamout
    Boolean compress
    File? gga_vcf
    File? gga_vcf_idx
    File? variants_for_contamination
    File? variants_for_contamination_idx

    File? gatk_override

    String? gcs_project_for_requester_pays

    # runtime
    String gatk_docker
    Int? mem
    Int? preemptible
    Int? max_retries
    Int? disk_space
    Int? cpu
    Boolean use_ssd = false
  }

  String output_vcf = "output" + if compress then ".vcf.gz" else ".vcf"
  String output_vcf_idx = output_vcf + if compress then ".tbi" else ".idx"

  String output_stats = output_vcf + ".stats"

  # Mem is in units of GB but our command and memory runtime values are in MB
  Int machine_mem = if defined(mem) then mem * 1000 else 3500
  Int command_mem = machine_mem - 500

  parameter_meta {
    intervals: {localization_optional: true}
    ref_fasta: {localization_optional: true}
    ref_fai: {localization_optional: true}
    ref_dict: {localization_optional: true}
    tumor_bam: {localization_optional: true}
    tumor_bai: {localization_optional: true}
    normal_bam: {localization_optional: true}
    normal_bai: {localization_optional: true}
    pon: {localization_optional: true}
    pon_idx: {localization_optional: true}
    gnomad: {localization_optional: true}
    gnomad_idx: {localization_optional: true}
    gga_vcf: {localization_optional: true}
    gga_vcf_idx: {localization_optional: true}
    variants_for_contamination: {localization_optional: true}
    variants_for_contamination_idx: {localization_optional: true}
  }

  command <<<
      set -e

      export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

      # We need to create these files regardless, even if they stay empty
      touch bamout.bam
      touch f1r2.tar.gz
      echo "" > normal_name.txt

      gatk --java-options "-Xmx~{command_mem}m" GetSampleName -R ~{ref_fasta} -I ~{tumor_bam} -O tumor_name.txt -encode \
      ~{"--gcs-project-for-requester-pays " + gcs_project_for_requester_pays}
      tumor_command_line="-I ~{tumor_bam} -tumor `cat tumor_name.txt`"

      if [[ ! -z "~{normal_bam}" ]]; then
          gatk --java-options "-Xmx~{command_mem}m" GetSampleName -R ~{ref_fasta} -I ~{normal_bam} -O normal_name.txt -encode \
          ~{"--gcs-project-for-requester-pays " + gcs_project_for_requester_pays}
          normal_command_line="-I ~{normal_bam} -normal `cat normal_name.txt`"
      fi

      gatk --java-options "-Xmx~{command_mem}m" Mutect2 \
          -R ~{ref_fasta} \
          $tumor_command_line \
          $normal_command_line \
          ~{"--germline-resource " + gnomad} \
          ~{"-pon " + pon} \
          ~{"-L " + intervals} \
          ~{"--alleles " + gga_vcf} \
          -O "~{output_vcf}" \
          ~{true='--bam-output bamout.bam' false='' make_bamout} \
          ~{m2_extra_args} \
          ~{"--gcs-project-for-requester-pays " + gcs_project_for_requester_pays}

      m2_exit_code=$?

      ### GetPileupSummaries

      # If the variants for contamination and the intervals for this scatter don't intersect, GetPileupSummaries
      # throws an error.  However, there is nothing wrong with an empty intersection for our purposes; it simply doesn't
      # contribute to the merged pileup summaries that we create downstream.  We implement this by with array outputs.
      # If the tool errors, no table is created and the glob yields an empty array.
      set +e

      if [[ ! -z "~{variants_for_contamination}" ]]; then
          gatk --java-options "-Xmx~{command_mem}m" GetPileupSummaries -R ~{ref_fasta} -I ~{tumor_bam} ~{"--interval-set-rule INTERSECTION -L " + intervals} \
              -V ~{variants_for_contamination} -L ~{variants_for_contamination} -O tumor-pileups.table ~{getpileupsummaries_extra_args} \
              ~{"--gcs-project-for-requester-pays " + gcs_project_for_requester_pays}


          if [[ ! -z "~{normal_bam}" ]]; then
              gatk --java-options "-Xmx~{command_mem}m" GetPileupSummaries -R ~{ref_fasta} -I ~{normal_bam} ~{"--interval-set-rule INTERSECTION -L " + intervals} \
                  -V ~{variants_for_contamination} -L ~{variants_for_contamination} -O normal-pileups.table ~{getpileupsummaries_extra_args} \
                  ~{"--gcs-project-for-requester-pays " + gcs_project_for_requester_pays}
          fi
      fi

      # the script only fails if Mutect2 itself fails
      exit $m2_exit_code
  >>>

  runtime {
    docker: gatk_docker
    bootDiskSizeGb: 12
    memory: machine_mem + " MB"
    disks: "local-disk " + select_first([disk_space, 100]) + if use_ssd then " SSD" else " HDD"
    preemptible: select_first([preemptible, 10])
    maxRetries: select_first([max_retries, 0])
    cpu: select_first([cpu, 1])
  }

  output {
    File unfiltered_vcf = "~{output_vcf}"
    File unfiltered_vcf_idx = "~{output_vcf_idx}"
    File output_bamOut = "bamout.bam"
    String tumor_sample = read_string("tumor_name.txt")
    String normal_sample = read_string("normal_name.txt")
    File stats = "~{output_stats}"
    File f1r2_counts = "f1r2.tar.gz"
    Array[File] tumor_pileups = glob("*tumor-pileups.table")
    Array[File] normal_pileups = glob("*normal-pileups.table")
  }
}

task MergeStats {
  input {
    Array[File]+ stats
    Runtime runtime_params
  }

  command <<<
    set -e
    export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.gatk_override}

    gatk --java-options "-Xmx~{runtime_params.command_mem}m" MergeMutectStats \
      -stats ~{sep=" -stats " stats} -O merged.stats
  >>>

  runtime {
    docker: runtime_params.gatk_docker
    bootDiskSizeGb: runtime_params.boot_disk_size
    memory: runtime_params.machine_mem + " MB"
    disks: "local-disk " + runtime_params.disk + " HDD"
    preemptible: runtime_params.preemptible
    maxRetries: runtime_params.max_retries
    cpu: runtime_params.cpu
  }

  output {
    File merged_stats = "merged.stats"
  }
}

task Filter {
  input {
    File? intervals
    File ref_fasta
    File ref_fai
    File ref_dict
    File unfiltered_vcf
    File unfiltered_vcf_idx
    String output_name
    Boolean compress
    File? mutect_stats
    File? artifact_priors_tar_gz
    File? contamination_table
    File? maf_segments
    String? m2_extra_filtering_args

    Runtime runtime_params
    Int? disk_space
  }

  String output_vcf = output_name + if compress then ".vcf.gz" else ".vcf"
  String output_vcf_idx = output_vcf + if compress then ".tbi" else ".idx"

  parameter_meta{
    ref_fasta: {localization_optional: true}
    ref_fai: {localization_optional: true}
    ref_dict: {localization_optional: true}
  }

  command <<<
    set -e

    export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.gatk_override}

    gatk --java-options "-Xmx~{runtime_params.command_mem}m" FilterMutectCalls -V ~{unfiltered_vcf} \
      -R ~{ref_fasta} \
      -O ~{output_vcf} \
      ~{"--contamination-table " + contamination_table} \
      ~{"--tumor-segmentation " + maf_segments} \
      ~{"--ob-priors " + artifact_priors_tar_gz} \
      ~{"-stats " + mutect_stats} \
      --filtering-stats filtering.stats \
      ~{m2_extra_filtering_args}
  >>>

  runtime {
    docker: runtime_params.gatk_docker
    bootDiskSizeGb: runtime_params.boot_disk_size
    memory: runtime_params.machine_mem + " MB"
    disks: "local-disk " + select_first([disk_space, runtime_params.disk]) + " HDD"
    preemptible: runtime_params.preemptible
    maxRetries: runtime_params.max_retries
    cpu: runtime_params.cpu
  }

  output {
    File filtered_vcf = "~{output_vcf}"
    File filtered_vcf_idx = "~{output_vcf_idx}"
    File filtering_stats = "filtering.stats"
  }
}

task Funcotate {
  input {
    File ref_fasta
    File ref_fai
    File ref_dict
    File input_vcf
    File input_vcf_idx
    String reference_version
    String output_file_base_name
    String output_format
    Boolean compress
    Boolean use_gnomad
    # This should be updated when a new version of the data sources is released
    # TODO: Make this dynamically chosen in the command.
    File? data_sources_tar_gz = "gs://broad-public-datasets/funcotator/funcotator_dataSources.v1.7.20200521s.tar.gz"
    String? control_id
    String? case_id
    String? sequencing_center
    String? sequence_source
    String? transcript_selection_mode
    File? transcript_selection_list
    Array[String]? annotation_defaults
    Array[String]? annotation_overrides
    Array[String]? funcotator_excluded_fields
    Boolean? filter_funcotations
    File? interval_list

    String? extra_args
    String? gcs_project_for_requester_pays

    # ==============
    Runtime runtime_params
    Int? disk_space   #override to request more disk than default small task params

    # You may have to change the following two parameter values depending on the task requirements
    Int default_ram_mb = 3000
    # WARNING: In the workflow, you should calculate the disk space as an input to this task (disk_space_gb).  Please see [TODO: Link from Jose] for examples.
    Int default_disk_space_gb = 100
  }

  # ==============
  # Process input args:
  String output_maf = output_file_base_name + ".maf"
  String output_maf_index = output_maf + ".idx"
  String output_vcf = output_file_base_name + if compress then ".vcf.gz" else ".vcf"
  String output_vcf_idx = output_vcf +  if compress then ".tbi" else ".idx"
  String output_file = if output_format == "MAF" then output_maf else output_vcf
  String output_file_index = if output_format == "MAF" then output_maf_index else output_vcf_idx
  String transcript_selection_arg = if defined(transcript_selection_list) then " --transcript-list " else ""
  String annotation_def_arg = if defined(annotation_defaults) then " --annotation-default " else ""
  String annotation_over_arg = if defined(annotation_overrides) then " --annotation-override " else ""
  String filter_funcotations_args = if defined(filter_funcotations) && (filter_funcotations) then " --remove-filtered-variants " else ""
  String excluded_fields_args = if defined(funcotator_excluded_fields) then " --exclude-field " else ""
  String interval_list_arg = if defined(interval_list) then " -L " else ""
  String extra_args_arg = select_first([extra_args, ""])

  String dollar = "$"

  parameter_meta{
    ref_fasta: {localization_optional: true}
    ref_fai: {localization_optional: true}
    ref_dict: {localization_optional: true}
    input_vcf: {localization_optional: true}
    input_vcf_idx: {localization_optional: true}
  }

  command <<<
      set -e
      export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.gatk_override}

      # Extract our data sources:
      echo "Extracting data sources zip file..."
      mkdir datasources_dir
      tar zxvf ~{data_sources_tar_gz} -C datasources_dir --strip-components 1
      DATA_SOURCES_FOLDER="$PWD/datasources_dir"

      # Handle gnomAD:
      if ~{use_gnomad} ; then
          echo "Enabling gnomAD..."
          for potential_gnomad_gz in gnomAD_exome.tar.gz gnomAD_genome.tar.gz ; do
              if [[ -f ~{dollar}{DATA_SOURCES_FOLDER}/~{dollar}{potential_gnomad_gz} ]] ; then
                  cd ~{dollar}{DATA_SOURCES_FOLDER}
                  tar -zvxf ~{dollar}{potential_gnomad_gz}
                  cd -
              else
                  echo "ERROR: Cannot find gnomAD folder: ~{dollar}{potential_gnomad_gz}" 1>&2
                  false
              fi
          done
      fi

      # Run Funcotator:
      gatk --java-options "-Xmx~{runtime_params.command_mem}m" Funcotator \
          --data-sources-path $DATA_SOURCES_FOLDER \
          --ref-version ~{reference_version} \
          --output-file-format ~{output_format} \
          -R ~{ref_fasta} \
          -V ~{input_vcf} \
          -O ~{output_file} \
          ~{interval_list_arg} ~{default="" interval_list} \
          --annotation-default normal_barcode:~{default="Unknown" control_id} \
          --annotation-default tumor_barcode:~{default="Unknown" case_id} \
          --annotation-default Center:~{default="Unknown" sequencing_center} \
          --annotation-default source:~{default="Unknown" sequence_source} \
          ~{"--transcript-selection-mode " + transcript_selection_mode} \
          ~{transcript_selection_arg}~{default="" sep=" --transcript-list " transcript_selection_list} \
          ~{annotation_def_arg}~{default="" sep=" --annotation-default " annotation_defaults} \
          ~{annotation_over_arg}~{default="" sep=" --annotation-override " annotation_overrides} \
          ~{excluded_fields_args}~{default="" sep=" --exclude-field " funcotator_excluded_fields} \
          ~{filter_funcotations_args} \
          ~{extra_args_arg} \
          ~{"--gcs-project-for-requester-pays " + gcs_project_for_requester_pays}
      # Make sure we have a placeholder index for MAF files so this workflow doesn't fail:
      if [[ "~{output_format}" == "MAF" ]] ; then
        touch ~{output_maf_index}
      fi
  >>>

  runtime {
    docker: runtime_params.gatk_docker
    bootDiskSizeGb: runtime_params.boot_disk_size
    memory: runtime_params.machine_mem + " MB"
    disks: "local-disk " + select_first([disk_space, runtime_params.disk]) + " HDD"
    preemptible: runtime_params.preemptible
    maxRetries: runtime_params.max_retries
    cpu: runtime_params.cpu
  }

  output {
    File funcotated_output_file = "~{output_file}"
    File funcotated_output_file_index = "~{output_file_index}"
  }
}

task VariantFiltration {
  input {
    File input_vcf
    File input_vcf_index
    String base_name

    File ref_dict
    File ref_fasta
    File ref_fasta_index

    String? gatk_override
    String docker
    Int? preemptible_count
    Int? memory
    Int? disk
  }

  command <<<
    export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

    gatk \
      VariantFiltration \
    --R ~{ref_fasta} \
    --V ~{input_vcf} \
    --window 35 \
    --cluster 3 \
    --filter-name "FS" \
    --filter "FS > 30.0" \
    --filter-name "QD" \
    --filter "QD < 2.0" \
    -O ~{base_name}
  >>>

  output {
    File output_vcf = "~{base_name}"
    File output_vcf_index = "~{base_name}.tbi"
  }

  runtime {
    docker: docker
    bootDiskSizeGb: "32"
    memory: select_first([memory,6])+" GB"
    disks: "local-disk "+select_first([disk, 250])+" HDD"
    preemptible: select_first([preemptible_count, 10])
  }
}

task picard_CleanAfterStar {
  input {
    File input_bam
    File input_bai
    String base_name

    File ref_dict
    File ref_fasta
    File ref_fasta_index

    String? SO="coordinate"
    String? RGID=base_name
    String? RGLB="TruSeq_DNA_prep"
    String? RGPL="illumina"
    String? RGPU="hiseq2000"
    String? RGSM=base_name
    String? VALIDATION_STRINGENCY="SILENT"

    String? gatk_override
    String docker
    Int? preemptible_count
    Int? memory
    Int? disk
  }

  command <<<
    export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

    gatk \
      AddOrReplaceReadGroups \
      -I ~{input_bam} \
      -O ~{base_name}_rg.bam \
      -SO ~{SO} \
      --RGID ~{RGID} \
      --RGLB ~{RGLB} \
      --RGPL ~{RGPL} \
      --RGPU ~{RGPU} \
      --RGSM ~{RGSM} \
      --CREATE_INDEX true

    gatk \
      MarkDuplicates \
      -I ~{base_name}_rg.bam \
      -O ~{base_name}_rg_md.bam \
      --VALIDATION_STRINGENCY ~{VALIDATION_STRINGENCY} \
      --CREATE_INDEX true \
      -M ~{base_name}_rg_md.metrics
  >>>

  output {
    File output_bam = "~{base_name}_rg_md.bam"
    File output_bam_index = "~{base_name}_rg_md.bai"
    File output_metrics = "~{base_name}_rg_md.metrics"
  }

  runtime {
    docker: docker
    bootDiskSizeGb: "32"
    memory: select_first([memory,6])+" GB"
    disks: "local-disk "+select_first([disk, 250])+" HDD"
    preemptible: select_first([preemptible_count, 10])
  }
}
