
import "https://raw.githubusercontent.com/BroadInstitute/ccle_processing/CNV_pipeline/wdl_scripts/BamToUnmappedRGBams.wdl" as BamToUnmappedBam
import "https://raw.githubusercontent.com/gatk-workflows/gatk4-somatic-cnvs/1.3.0/cnv_common_tasks.wdl" as CNVTasks
import "https://raw.githubusercontent.com/gatk-workflows/gatk4-somatic-cnvs/1.3.0/cnv_somatic_oncotator_workflow.wdl" as CNVOncotator
import "https://raw.githubusercontent.com/BroadInstitute/ccle_processing/CNV_pipeline/wdl_scripts/CNV_Somatic_Pair_Workflow.wdl" as CNV_GATK
import "https://raw.githubusercontent.com/BroadInstitute/ccle_processing/CNV_pipeline/wdl_scripts/PreProcessingForVariantDiscovery_GATK4.wdl" as PreProcessingForVariantDiscovery_GATK4
import "https://raw.githubusercontent.com/BroadInstitute/ccle_processing/CNV_pipeline/wdl_scripts/Manta_SomaticSV.wdl" as Manta_SomaticSV
import "https://api.firecloud.org/ga4gh/v1/tools/gatk:cnn-variant-common-tasks/versions/1/plain-WDL/descriptor" as CNNTasks

workflow CCLE_CNV_Pipeline {

	File ref_fasta
	File ref_fasta_index
	File ref_dict

	File input_bam

	String picard_path
	String picard_docker

	Int preemptible_tries

	##REMAPPING
	String sample_name
	String ref_name

	String unmapped_bam_suffix

	String bwa_commandline
	Int compression_level

	File dbSNP_vcf
	File dbSNP_vcf_index
	Array[File] known_indels_sites_VCFs
	Array[File] known_indels_sites_indices

	String gotc_docker
	String gatk_docker
	String python_docker

	String gotc_path
	String gatk_launch_path

	Int flowcell_small_disk
	Int flowcell_medium_disk
	Int agg_small_disk
	Int agg_medium_disk
	Int agg_large_disk

	Int agg_preemptible_tries

	String base_file_name = sample_name + "." + ref_name

	##### CNN HC Mutations
	File resource_fofn               # File of VCF file names of resources of known SNPs and INDELs, (e.g. mills, gnomAD)
	File resource_fofn_index         # File of VCF file indices of resources
	File? architecture_json          # Neural Net configuration for CNNScoreVariants
	File? architecture_hd5           # Pre-Trained weights and architecture for CNNScoreVariants
	Int? inference_batch_size        # Batch size for python in CNNScoreVariants
	Int? transfer_batch_size         # Batch size for java in CNNScoreVariants
	Int? intra_op_threads            # Tensorflow threading within nodes
	Int? inter_op_threads            # Tensorflow threading between nodes
	String output_prefix             # Identifying string for this run will be used to name all output files
	String? tensor_type              # What kind of tensors the Neural Net expects (e.g. reference, read_tensor)
	String info_key                  # The score key for the info field of the vcf (e.g. CNN_1D, CNN_2D)
	String snp_tranches              # Filtering threshold(s) for SNPs in terms of sensitivity to overlapping known variants in resources
	String indel_tranches            # Filtering threshold(s) for INDELs in terms of sensitivity to overlapping known variants in resources
	File? gatk_override              # GATK Jar file to over ride the one included in gatk_docker
	File calling_intervals
	Int scatter_count                # Number of shards for parallelization of HaplotypeCaller and CNNScoreVariants
	String extra_args                # Extra arguments for HaplotypeCaller

	# Runtime parameters
	Int? mem_gb
	Int? preemptible_attempts
	Float? disk_space_gb
	Int? cpu

	Int? increase_disk_size

	### GATK CNV

	##################################
	#### required basic arguments ####
	##################################
	File common_sites
	File intervals
	File? blacklist_intervals
	File? normal_bam
	File? normal_bam_idx
	File read_count_pon

	##################################
	#### optional basic arguments ####
	##################################
    # For running oncotator
	Boolean? is_run_oncotator
	File? gatk4_jar_override
	Int? preemptible_attempts
	# Use as a last resort to increase the disk given to every task in case of ill behaving data
	Int? emergency_extra_disk

	####################################################
	#### optional arguments for PreprocessIntervals ####
	####################################################
	Int? padding
	Int? bin_length
	Int? mem_gb_for_preprocess_intervals

	##############################################
	#### optional arguments for CollectCounts ####
	##############################################
	String? collect_counts_format
	Int? mem_gb_for_collect_counts

	#####################################################
	#### optional arguments for CollectAllelicCounts ####
	#####################################################
	String? minimum_base_quality
	Int? mem_gb_for_collect_allelic_counts

	##################################################
	#### optional arguments for DenoiseReadCounts ####
	##################################################
	Int? number_of_eigensamples
	Int? mem_gb_for_denoise_read_counts

	##############################################
	#### optional arguments for ModelSegments ####
	##############################################
	Int? max_num_segments_per_chromosome
	Int? min_total_allele_count
	Int? min_total_allele_count_normal
	Float? genotyping_homozygous_log_ratio_threshold
	Float? genotyping_base_error_rate
	Float? kernel_variance_copy_ratio
	Float? kernel_variance_allele_fraction
	Float? kernel_scaling_allele_fraction
	Int? kernel_approximation_dimension
	Array[Int]+? window_sizes = [8, 16, 32, 64, 128, 256]
	Float? num_changepoints_penalty_factor
	Float? minor_allele_fraction_prior_alpha
	Int? num_samples_copy_ratio
	Int? num_burn_in_copy_ratio
	Int? num_samples_allele_fraction
	Int? num_burn_in_allele_fraction
	Float? smoothing_threshold_copy_ratio
	Float? smoothing_threshold_allele_fraction
	Int? max_num_smoothing_iterations
	Int? num_smoothing_iterations_per_fit
	Int? mem_gb_for_model_segments

	######################################################
	#### optional arguments for CallCopyRatioSegments ####
	######################################################
	Float? neutral_segment_copy_ratio_lower_bound
	Float? neutral_segment_copy_ratio_upper_bound
	Float? outlier_neutral_segment_copy_ratio_z_score_threshold
	Float? calling_copy_ratio_z_score_threshold
	Int? mem_gb_for_call_copy_ratio_segments

	#########################################
	#### optional arguments for plotting ####
	#########################################
	Int? minimum_contig_length
	Int? mem_gb_for_plotting

	##########################################
	#### optional arguments for Oncotator ####
	##########################################
	String? additional_args_for_oncotator
	String? oncotator_docker
	Int? mem_gb_for_oncotator
	Int? boot_disk_space_gb_for_oncotator

	Int ref_size = ceil(size(ref_fasta, "GB") + size(ref_fasta_dict, "GB") + size(ref_fasta_index, "GB"))
	Int read_count_pon_size = ceil(size(read_count_pon, "GB"))
	Int normal_bam_size = if defined(normal_bam) then ceil(size(normal_bam, "GB") + size(normal_bam_idx, "GB")) else 0

	Int gatk4_override_size = if defined(gatk4_jar_override) then ceil(size(gatk4_jar_override, "GB")) else 0
	# This is added to every task as padding, should increase if systematically you need more disk for every call
	Int disk_pad = 20 + ceil(size(intervals, "GB")) + ceil(size(common_sites, "GB")) + gatk4_override_size + select_first([emergency_extra_disk, 0])

	File final_normal_bam = select_first([normal_bam, "null"])
	File final_normal_bam_idx = select_first([normal_bam_idx, "null"])

	Int preprocess_intervals_disk = ref_size + disk_pad

	Int additional_disk = select_first([increase_disk_size, 20])

	#### MANTA
	File? interval_list

	String manta_docker
	String config_manta

	Boolean is_exome = defined(interval_list)
	Boolean is_cram

	if (is_exome) {
    call ConvertToBedTabix {
      input: interval_list = interval_list
    }
	}

	# Revert input to unmapped
	call BamToUnmappedBam.RevertBamToUnmappedRGBams as RevertBamToUnmappedRGBams {
    input:
			input_bam = input_bam,
			picard_path = picard_path,
			docker_image = picard_docker
	}

	scatter (unmapped_bam in RevertBamToUnmappedRGBams.unmapped_bams) {

    # Get the basename, i.e. strip the filepath and the extension
    String bam_basename = basename(unmapped_bam, ".bam")

    # Sort the BAM records
    call BamToUnmappedBam.SortBamByQueryname as SortBamByQueryname {
			input:
        input_bam = unmapped_bam,
        sorted_bam_name = bam_basename + ".unmapped.bam",
        picard_path = picard_path,
        docker_image = picard_docker,
        preemptible_tries = preemptible_tries
    }

    # ValidateSamFile
    call BamToUnmappedBam.ValidateSamFile as ValidateSamFile {
			input:
        input_bam = SortBamByQueryname.sorted_bam,
        report_filename = bam_basename + ".unmapped.validation_report",
        picard_path = picard_path,
        docker_image = picard_docker,
        preemptible_tries = preemptible_tries
    }
	}


	call CreateTxt {
    input:
      array_of_files = SortBamByQueryname.sorted_bam
      list_name = sample_name
	}

	Array[File] flowcell_unmapped_bams = read_lines(CreateTxt.file_list_name)

	# Get the version of BWA to include in the PG record in the header of the BAM produced
	# by MergeBamAlignment.
	call PreProcessingForVariantDiscovery_GATK4.GetBwaVersion as GetBwaVersion {
    input:
			docker_image = gotc_docker,
			bwa_path = gotc_path,
			preemptible_tries = preemptible_tries
	}

	# Align flowcell-level unmapped input bams in parallel
	scatter (unmapped_bam in flowcell_unmapped_bams) {

    # Get the basename, i.e. strip the filepath and the extension
    String bam_basename = basename(unmapped_bam, unmapped_bam_suffix)

    # Map reads to reference
    call PreProcessingForVariantDiscovery_GATK4.SamToFastqAndBwaMem as SamToFastqAndBwaMem {
			input:
        input_bam = unmapped_bam,
        bwa_commandline = bwa_commandline,
        output_bam_basename = bam_basename + ".unmerged",
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        docker_image = gotc_docker,
        bwa_path = gotc_path,
        picard_path = gotc_path,
        disk_size = flowcell_medium_disk,
        preemptible_tries = preemptible_tries,
        compression_level = compression_level
		 }

    # Merge original uBAM and BWA-aligned BAM
    call PreProcessingForVariantDiscovery_GATK4.MergeBamAlignment as MergeBamAlignment {
			input:
        unmapped_bam = unmapped_bam,
        bwa_commandline = bwa_commandline,
        bwa_version = GetBwaVersion.version,
        aligned_bam = SamToFastqAndBwaMem.output_bam,
        output_bam_basename = bam_basename + ".aligned.unsorted",
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        docker_image = picard_docker,
        picard_path = picard_path,
        disk_size = flowcell_medium_disk,
        preemptible_tries = preemptible_tries,
        compression_level = compression_level
    }
	}

	# Aggregate aligned+merged flowcell BAM files and mark duplicates
	# We take advantage of the tool's ability to take multiple BAM inputs and write out a single output
	# to avoid having to spend time just merging BAM files.
	call PreProcessingForVariantDiscovery_GATK4.call as call MarkDuplicates {
    input:
			input_bams = MergeBamAlignment.output_bam,
			output_bam_basename = base_file_name + ".aligned.unsorted.duplicates_marked",
			metrics_filename = base_file_name + ".duplicate_metrics",
			docker_image = picard_docker,
			picard_path = picard_path,
			disk_size = agg_large_disk,
			compression_level = compression_level,
			preemptible_tries = agg_preemptible_tries
	}

	# Sort aggregated+deduped BAM file and fix tags
	call PreProcessingForVariantDiscovery_GATK4.call as call SortAndFixTags {
    input:
			input_bam = MarkDuplicates.output_bam,
			output_bam_basename = base_file_name + ".aligned.duplicate_marked.sorted",
			ref_dict = ref_dict,
			ref_fasta = ref_fasta,
			ref_fasta_index = ref_fasta_index,
			docker_image = picard_docker,
			picard_path = picard_path,
			disk_size = agg_large_disk,
			preemptible_tries = 0,
			compression_level = compression_level
	}

	# Create list of sequences for scatter-gather parallelization
	call PreProcessingForVariantDiscovery_GATK4.CreateSequenceGroupingTSV as CreateSequenceGroupingTSV {
    input:
			ref_dict = ref_dict,
			docker_image = python_docker,
			preemptible_tries = preemptible_tries
	}

	# Perform Base Quality Score Recalibration (BQSR) on the sorted BAM in parallel
	scatter (subgroup in CreateSequenceGroupingTSV.sequence_grouping) {
    # Generate the recalibration model by interval
    call PreProcessingForVariantDiscovery_GATK4.BaseRecalibrator as BaseRecalibrator {
			input:
        input_bam = SortAndFixTags.output_bam,
        input_bam_index = SortAndFixTags.output_bam_index,
        recalibration_report_filename = base_file_name + ".recal_data.csv",
        sequence_group_interval = subgroup,
        dbSNP_vcf = dbSNP_vcf,
        dbSNP_vcf_index = dbSNP_vcf_index,
        known_indels_sites_VCFs = known_indels_sites_VCFs,
        known_indels_sites_indices = known_indels_sites_indices,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        docker_image = gatk_docker,
        gatk_launch_path = gatk_launch_path,
        disk_size = agg_small_disk,
        preemptible_tries = agg_preemptible_tries
    }
	}

	# Merge the recalibration reports resulting from by-interval recalibration
	call PreProcessingForVariantDiscovery_GATK4.GatherBqsrReports as GatherBqsrReports {
    input:
			input_bqsr_reports = BaseRecalibrator.recalibration_report,
			output_report_filename = base_file_name + ".recal_data.csv",
			docker_image = gatk_docker,
			gatk_launch_path = gatk_launch_path,
			disk_size = flowcell_small_disk,
			preemptible_tries = preemptible_tries
	}

	scatter (subgroup in CreateSequenceGroupingTSV.sequence_grouping_with_unmapped) {

    # Apply the recalibration model by interval
    call PreProcessingForVariantDiscovery_GATK4.ApplyBQSR as ApplyBQSR {
			input:
        input_bam = SortAndFixTags.output_bam,
        input_bam_index = SortAndFixTags.output_bam_index,
        output_bam_basename = base_file_name + ".aligned.duplicates_marked.recalibrated",
        recalibration_report = GatherBqsrReports.output_bqsr_report,
        sequence_group_interval = subgroup,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        docker_image = gatk_docker,
        gatk_launch_path = gatk_launch_path,
        disk_size = agg_small_disk,
        preemptible_tries = agg_preemptible_tries
    }
	}

	# Merge the recalibrated BAM files resulting from by-interval recalibration
	call PreProcessingForVariantDiscovery_GATK4.GatherBamFiles as GatherBamFiles {
    input:
			input_bams = ApplyBQSR.recalibrated_bam,
			output_bam_basename = base_file_name,
			docker_image = picard_docker,
			picard_path = picard_path,
			disk_size = agg_large_disk,
			preemptible_tries = agg_preemptible_tries,
			compression_level = compression_level
	}

	call CNVTasks.PreprocessIntervals {
    input:
      intervals = intervals,
      blacklist_intervals = blacklist_intervals,
      ref_fasta = ref_fasta,
      ref_fasta_fai = ref_fasta_index,
      ref_fasta_dict = ref_fasta_dict,
      padding = padding,
      bin_length = bin_length,
      gatk4_jar_override = gatk4_jar_override,
      gatk_docker = gatk_docker,
      mem_gb = mem_gb_for_preprocess_intervals,
      disk_space_gb = preprocess_intervals_disk,
      preemptible_attempts = preemptible_attempts
	}

	Int tumor_bam_size = ceil(size(GatherBamFiles.output_bam, "GB") + size(GatherBamFiles.output_bam_index, "GB"))
	Int collect_counts_tumor_disk = tumor_bam_size + ceil(size(PreprocessIntervals.preprocessed_intervals, "GB")) + disk_pad

	call CNVTasks.CollectCounts as CollectCountsTumor {
    input:
      intervals = PreprocessIntervals.preprocessed_intervals,
      bam = GatherBamFiles.output_bam,
      bam_idx = GatherBamFiles.output_bam_index,
      ref_fasta = ref_fasta,
      ref_fasta_fai = ref_fasta_index,
      ref_fasta_dict = ref_fasta_dict,
      format = collect_counts_format,
      gatk4_jar_override = gatk4_jar_override,
      gatk_docker = gatk_docker,
      mem_gb = mem_gb_for_collect_counts,
      disk_space_gb = collect_counts_tumor_disk,
      preemptible_attempts = preemptible_attempts
	}

	Int collect_allelic_counts_tumor_disk = tumor_bam_size + ref_size + disk_pad
	call CNVTasks.CollectAllelicCounts as CollectAllelicCountsTumor {
    input:
      common_sites = common_sites,
      bam = GatherBamFiles.output_bam,
      bam_idx = GatherBamFiles.output_bam_index,
      ref_fasta = ref_fasta,
      ref_fasta_dict = ref_fasta_dict,
      ref_fasta_fai = ref_fasta_index,
      minimum_base_quality =  minimum_base_quality,
      gatk4_jar_override = gatk4_jar_override,
      gatk_docker = gatk_docker,
      mem_gb = mem_gb_for_collect_allelic_counts,
      disk_space_gb = collect_allelic_counts_tumor_disk,
      preemptible_attempts = preemptible_attempts
	}

	Int denoise_read_counts_tumor_disk = read_count_pon_size + ceil(size(CollectCountsTumor.counts, "GB")) + disk_pad
	call CNV_GATK.DenoiseReadCounts as DenoiseReadCountsTumor {
    input:
      entity_id = CollectCountsTumor.entity_id,
      read_counts = CollectCountsTumor.counts,
      read_count_pon = read_count_pon,
      number_of_eigensamples = number_of_eigensamples,
      gatk4_jar_override = gatk4_jar_override,
      gatk_docker = gatk_docker,
      mem_gb = mem_gb_for_denoise_read_counts,
      disk_space_gb = denoise_read_counts_tumor_disk,
      preemptible_attempts = preemptible_attempts
  }

	Int model_segments_normal_portion = if defined(normal_bam) then ceil(size(CollectAllelicCountsNormal.allelic_counts, "GB")) else 0
	Int model_segments_tumor_disk = ceil(size(DenoiseReadCountsTumor.denoised_copy_ratios, "GB")) + ceil(size(CollectAllelicCountsTumor.allelic_counts, "GB")) + model_segments_normal_portion + disk_pad
	call CNV_GATK.ModelSegments as ModelSegmentsTumor {
    input:
      entity_id = CollectCountsTumor.entity_id,
      denoised_copy_ratios = DenoiseReadCountsTumor.denoised_copy_ratios,
      allelic_counts = CollectAllelicCountsTumor.allelic_counts,
      normal_allelic_counts = CollectAllelicCountsNormal.allelic_counts,
      max_num_segments_per_chromosome = max_num_segments_per_chromosome,
      min_total_allele_count = min_total_allele_count,
      min_total_allele_count_normal = min_total_allele_count_normal,
      genotyping_homozygous_log_ratio_threshold = genotyping_homozygous_log_ratio_threshold,
      genotyping_base_error_rate = genotyping_base_error_rate,
      kernel_variance_copy_ratio = kernel_variance_copy_ratio,
      kernel_variance_allele_fraction = kernel_variance_allele_fraction,
      kernel_scaling_allele_fraction = kernel_scaling_allele_fraction,
      kernel_approximation_dimension = kernel_approximation_dimension,
      window_sizes = window_sizes,
      num_changepoints_penalty_factor = num_changepoints_penalty_factor,
      minor_allele_fraction_prior_alpha = minor_allele_fraction_prior_alpha,
      num_samples_copy_ratio = num_samples_copy_ratio,
      num_burn_in_copy_ratio = num_burn_in_copy_ratio,
      num_samples_allele_fraction = num_samples_allele_fraction,
      num_burn_in_allele_fraction = num_burn_in_allele_fraction,
      smoothing_threshold_copy_ratio = smoothing_threshold_copy_ratio,
      smoothing_threshold_allele_fraction = smoothing_threshold_allele_fraction,
      max_num_smoothing_iterations = max_num_smoothing_iterations,
      num_smoothing_iterations_per_fit = num_smoothing_iterations_per_fit,
      gatk4_jar_override = gatk4_jar_override,
      gatk_docker = gatk_docker,
      mem_gb = mem_gb_for_model_segments,
      disk_space_gb = model_segments_tumor_disk,
      preemptible_attempts = preemptible_attempts
	}

	Int copy_ratio_segments_tumor_disk = ceil(size(DenoiseReadCountsTumor.denoised_copy_ratios, "GB")) + ceil(size(ModelSegmentsTumor.copy_ratio_only_segments, "GB")) + disk_pad
	call CNV_GATK.CallCopyRatioSegments as CallCopyRatioSegmentsTumor {
    input:
      entity_id = CollectCountsTumor.entity_id,
      copy_ratio_segments = ModelSegmentsTumor.copy_ratio_only_segments,
      neutral_segment_copy_ratio_lower_bound = neutral_segment_copy_ratio_lower_bound,
      neutral_segment_copy_ratio_upper_bound = neutral_segment_copy_ratio_upper_bound,
      outlier_neutral_segment_copy_ratio_z_score_threshold = outlier_neutral_segment_copy_ratio_z_score_threshold,
      calling_copy_ratio_z_score_threshold = calling_copy_ratio_z_score_threshold,
      gatk4_jar_override = gatk4_jar_override,
      gatk_docker = gatk_docker,
      mem_gb = mem_gb_for_call_copy_ratio_segments,
      disk_space_gb = copy_ratio_segments_tumor_disk,
      preemptible_attempts = preemptible_attempts
	}

	# The F=files from other tasks are small enough to just combine into one disk variable and pass to the tumor plotting tasks
	Int plot_tumor_disk = ref_size + ceil(size(DenoiseReadCountsTumor.standardized_copy_ratios, "GB")) + ceil(size(DenoiseReadCountsTumor.denoised_copy_ratios, "GB")) + ceil(size(ModelSegmentsTumor.het_allelic_counts, "GB")) + ceil(size(ModelSegmentsTumor.modeled_segments, "GB")) + disk_pad
	call CNV_GATK.PlotDenoisedCopyRatios as PlotDenoisedCopyRatiosTumor {
    input:
      entity_id = CollectCountsTumor.entity_id,
      standardized_copy_ratios = DenoiseReadCountsTumor.standardized_copy_ratios,
      denoised_copy_ratios = DenoiseReadCountsTumor.denoised_copy_ratios,
      ref_fasta_dict = ref_fasta_dict,
      minimum_contig_length = minimum_contig_length,
      gatk4_jar_override = gatk4_jar_override,
      gatk_docker = gatk_docker,
      mem_gb = mem_gb_for_plotting,
      disk_space_gb = plot_tumor_disk,
      preemptible_attempts = preemptible_attempts
	}

	call CNV_GATK.PlotModeledSegments as PlotModeledSegmentsTumor {
    input:
      entity_id = CollectCountsTumor.entity_id,
      denoised_copy_ratios = DenoiseReadCountsTumor.denoised_copy_ratios,
      het_allelic_counts = ModelSegmentsTumor.het_allelic_counts,
      modeled_segments = ModelSegmentsTumor.modeled_segments,
      ref_fasta_dict = ref_fasta_dict,
      minimum_contig_length = minimum_contig_length,
      gatk4_jar_override = gatk4_jar_override,
      gatk_docker = gatk_docker,
      mem_gb = mem_gb_for_plotting,
      disk_space_gb = plot_tumor_disk,
      preemptible_attempts = preemptible_attempts
	}

	Int collect_counts_normal_disk = normal_bam_size + ceil(size(PreprocessIntervals.preprocessed_intervals, "GB")) + disk_pad
	if (defined(normal_bam)) {
    call CNVTasks.CollectCounts as CollectCountsNormal {
      input:
        intervals = PreprocessIntervals.preprocessed_intervals,
        bam = final_normal_bam,
        bam_idx = final_normal_bam_idx,
        ref_fasta = ref_fasta,
        ref_fasta_fai = ref_fasta_index,
        ref_fasta_dict = ref_fasta_dict,
        format = collect_counts_format,
        gatk4_jar_override = gatk4_jar_override,
        gatk_docker = gatk_docker,
        mem_gb = mem_gb_for_collect_counts,
        disk_space_gb = collect_counts_normal_disk,
        preemptible_attempts = preemptible_attempts
			}

			Int collect_allelic_counts_normal_disk = normal_bam_size + ref_size + disk_pad
			call CNVTasks.CollectAllelicCounts as CollectAllelicCountsNormal {
        input:
          common_sites = common_sites,
          bam = final_normal_bam,
          bam_idx = final_normal_bam_idx,
          ref_fasta = ref_fasta,
          ref_fasta_dict = ref_fasta_dict,
          ref_fasta_fai = ref_fasta_index,
          minimum_base_quality =  minimum_base_quality,
          gatk4_jar_override = gatk4_jar_override,
          gatk_docker = gatk_docker,
          mem_gb = mem_gb_for_collect_allelic_counts,
          disk_space_gb = collect_allelic_counts_normal_disk,
          preemptible_attempts = preemptible_attempts
			}

			Int denoise_read_counts_normal_disk = read_count_pon_size + ceil(size(CollectCountsNormal.counts, "GB")) + disk_pad
			call CNV_GATK.DenoiseReadCounts as DenoiseReadCountsNormal {
        input:
          entity_id = CollectCountsNormal.entity_id,
          read_counts = CollectCountsNormal.counts,
          read_count_pon = read_count_pon,
          number_of_eigensamples = number_of_eigensamples,
          gatk4_jar_override = gatk4_jar_override,
          gatk_docker = gatk_docker,
          mem_gb = mem_gb_for_denoise_read_counts,
          disk_space_gb = denoise_read_counts_normal_disk,
          preemptible_attempts = preemptible_attempts
			}

			Int model_segments_normal_disk = ceil(size(DenoiseReadCountsNormal.denoised_copy_ratios, "GB")) + ceil(size(CollectAllelicCountsNormal.allelic_counts, "GB")) + disk_pad
			call CNV_GATK.ModelSegments as ModelSegmentsNormal {
        input:
          entity_id = CollectCountsNormal.entity_id,
          denoised_copy_ratios = DenoiseReadCountsNormal.denoised_copy_ratios,
          allelic_counts = CollectAllelicCountsNormal.allelic_counts,
          max_num_segments_per_chromosome = max_num_segments_per_chromosome,
          min_total_allele_count = min_total_allele_count_normal,
          genotyping_homozygous_log_ratio_threshold = genotyping_homozygous_log_ratio_threshold,
          genotyping_base_error_rate = genotyping_base_error_rate,
          kernel_variance_copy_ratio = kernel_variance_copy_ratio,
          kernel_variance_allele_fraction = kernel_variance_allele_fraction,
          kernel_scaling_allele_fraction = kernel_scaling_allele_fraction,
          kernel_approximation_dimension = kernel_approximation_dimension,
          window_sizes = window_sizes,
          num_changepoints_penalty_factor = num_changepoints_penalty_factor,
          minor_allele_fraction_prior_alpha = minor_allele_fraction_prior_alpha,
          num_samples_copy_ratio = num_samples_copy_ratio,
          num_burn_in_copy_ratio = num_burn_in_copy_ratio,
          num_samples_allele_fraction = num_samples_allele_fraction,
          num_burn_in_allele_fraction = num_burn_in_allele_fraction,
          smoothing_threshold_copy_ratio = smoothing_threshold_copy_ratio,
          smoothing_threshold_allele_fraction = smoothing_threshold_allele_fraction,
          max_num_smoothing_iterations = max_num_smoothing_iterations,
          num_smoothing_iterations_per_fit = num_smoothing_iterations_per_fit,
          gatk4_jar_override = gatk4_jar_override,
          gatk_docker = gatk_docker,
          mem_gb = mem_gb_for_model_segments,
          disk_space_gb = model_segments_normal_disk,
          preemptible_attempts = preemptible_attempts
			}

			Int copy_ratio_segments_normal_disk = ceil(size(DenoiseReadCountsNormal.denoised_copy_ratios, "GB")) + ceil(size(ModelSegmentsNormal.copy_ratio_only_segments, "GB")) + disk_pad
			call CNV_GATK.CallCopyRatioSegments as CallCopyRatioSegmentsNormal {
        input:
          entity_id = CollectCountsNormal.entity_id,
          copy_ratio_segments = ModelSegmentsNormal.copy_ratio_only_segments,
          neutral_segment_copy_ratio_lower_bound = neutral_segment_copy_ratio_lower_bound,
          neutral_segment_copy_ratio_upper_bound = neutral_segment_copy_ratio_upper_bound,
          outlier_neutral_segment_copy_ratio_z_score_threshold = outlier_neutral_segment_copy_ratio_z_score_threshold,
          calling_copy_ratio_z_score_threshold = calling_copy_ratio_z_score_threshold,
          gatk4_jar_override = gatk4_jar_override,
          gatk_docker = gatk_docker,
          mem_gb = mem_gb_for_call_copy_ratio_segments,
          disk_space_gb = copy_ratio_segments_normal_disk,
          preemptible_attempts = preemptible_attempts
			}

			# The files from other tasks are small enough to just combine into one disk variable and pass to the normal plotting tasks
			Int plot_normal_disk = ref_size + ceil(size(DenoiseReadCountsNormal.standardized_copy_ratios, "GB")) + ceil(size(DenoiseReadCountsNormal.denoised_copy_ratios, "GB")) + ceil(size(ModelSegmentsNormal.het_allelic_counts, "GB")) + ceil(size(ModelSegmentsNormal.modeled_segments, "GB")) + disk_pad
			call CNV_GATK.PlotDenoisedCopyRatios as PlotDenoisedCopyRatiosNormal {
        input:
          entity_id = CollectCountsNormal.entity_id,
          standardized_copy_ratios = DenoiseReadCountsNormal.standardized_copy_ratios,
          denoised_copy_ratios = DenoiseReadCountsNormal.denoised_copy_ratios,
          ref_fasta_dict = ref_fasta_dict,
          minimum_contig_length = minimum_contig_length,
          gatk4_jar_override = gatk4_jar_override,
          gatk_docker = gatk_docker,
          mem_gb = mem_gb_for_plotting,
          disk_space_gb = plot_normal_disk,
          preemptible_attempts = preemptible_attempts
			}

			call CNV_GATK.PlotModeledSegments as PlotModeledSegmentsNormal {
        input:
          entity_id = CollectCountsNormal.entity_id,
          denoised_copy_ratios = DenoiseReadCountsNormal.denoised_copy_ratios,
          het_allelic_counts = ModelSegmentsNormal.het_allelic_counts,
          modeled_segments = ModelSegmentsNormal.modeled_segments,
          ref_fasta_dict = ref_fasta_dict,
          minimum_contig_length = minimum_contig_length,
          gatk4_jar_override = gatk4_jar_override,
          gatk_docker = gatk_docker,
          mem_gb = mem_gb_for_plotting,
          disk_space_gb = plot_normal_disk,
          preemptible_attempts = preemptible_attempts
			}
	}

	if (select_first([is_run_oncotator, false])) {
    call CNVOncotator.CNVOncotatorWorkflow as CNVOncotatorWorkflow {
			input:
        called_file = CallCopyRatioSegmentsTumor.called_copy_ratio_segments,
        additional_args = additional_args_for_oncotator,
        oncotator_docker = oncotator_docker,
        mem_gb_for_oncotator = mem_gb_for_oncotator,
        boot_disk_space_gb_for_oncotator = boot_disk_space_gb_for_oncotator,
        preemptible_attempts = preemptible_attempts
			}
	}

	call Manta_SomaticSV.Manta as Manta {
    input:
      sample_name = sample_name,
      tumor_bam = GatherBamFiles.output_bam,
      tumor_bam_index = GatherBamFiles.output_bam_index,
      normal_bam = normal_bam,
      normal_bam_index = normal_bam_idx,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      interval_bed = ConvertToBedTabix.out_interval,
      interval_bed_index = ConvertToBedTabix.out_interval_index,
      manta_docker = manta_docker,
      config_manta = config_manta,
      is_cram = is_cram
	}

	# Clunky check to see if the input is a BAM or a CRAM
	if (basename(GatherBamFiles.output_bam) == basename(GatherBamFiles.output_bam, ".bam")){
    call CNNTasks.CramToBam {
      input:
        reference_fasta = ref_fasta,
        reference_dict = ref_dict,
        reference_fasta_index = ref_fasta_index,
        cram_file = GatherBamFiles.output_bam,
        output_prefix = output_prefix,
        disk_space_gb = round(4*size(GatherBamFiles.output_bam, "GB") + ref_size + additional_disk),
        preemptible_attempts = preemptible_attempts
		}
	}

	call CNNTasks.SplitIntervals {
    input:
      gatk_override = gatk_override,
      scatter_count = scatter_count,
      intervals = calling_intervals,
      ref_fasta = ref_fasta,
      ref_dict = ref_dict,
      ref_fai = ref_fasta_index,
      gatk_docker = gatk_docker,
      disk_space = round(additional_disk + ref_size)
	}

	String mybam = select_first([CramToBam.output_bam, GatherBamFiles.output_bam])
	Float bam_size = size(mybam, "GB")

	scatter (calling_interval in SplitIntervals.interval_files) {
    call CNNTasks.RunHC4 {
      input:
        input_bam = mybam,
        input_bam_index = select_first([CramToBam.output_bam_index, GatherBamFiles.output_bam_index]),
        reference_fasta = ref_fasta,
        reference_dict = ref_dict,
        reference_fasta_index = ref_fasta_index,
        output_prefix = output_prefix,
        interval_list = calling_interval,
        gatk_docker = gatk_docker,
        gatk_override = gatk_override,
        preemptible_attempts = preemptible_attempts,
        extra_args = extra_args,
        disk_space_gb = round(bam_size + ref_size + additional_disk)
    }

    call CNNTasks.CNNScoreVariants {
      input:
        input_vcf = RunHC4.raw_vcf,
        input_vcf_index = RunHC4.raw_vcf_index,
        bam_file = RunHC4.bamout,
        bam_file_index = RunHC4.bamout_index,
        architecture_json = architecture_json,
        architecture_hd5 = architecture_hd5,
        reference_fasta = ref_fasta,
        tensor_type = tensor_type,
        inference_batch_size = inference_batch_size,
        transfer_batch_size = transfer_batch_size,
        intra_op_threads = intra_op_threads,
        inter_op_threads = inter_op_threads,
        reference_dict = ref_dict,
        reference_fasta_index = ref_fasta_index,
        output_prefix = output_prefix,
        interval_list = calling_interval,
        gatk_override = gatk_override,
        gatk_docker = gatk_docker,
        preemptible_attempts = preemptible_attempts,
        mem_gb = mem_gb,
        disk_space_gb = round((bam_size/scatter_count) + ref_size + additional_disk)
    }
	}

	call CNNTasks.MergeVCFs as MergeVCF_HC4 {
    input:
      input_vcfs = CNNScoreVariants.cnn_annotated_vcf,
      output_prefix = output_prefix,
      gatk_override = gatk_override,
      preemptible_attempts = preemptible_attempts,
      gatk_docker = gatk_docker,
      disk_space_gb = additional_disk
	}

	call CNNTasks.FilterVariantTranches {
    input:
      input_vcf = MergeVCF_HC4.merged_vcf,
      input_vcf_index = MergeVCF_HC4.merged_vcf_index,
      resource_fofn = resource_fofn,
      resource_fofn_index = resource_fofn_index,
      output_prefix = output_prefix,
      snp_tranches = snp_tranches,
      indel_tranches = indel_tranches,
      info_key = info_key,
      gatk_override = gatk_override,
      preemptible_attempts = preemptible_attempts,
      gatk_docker = gatk_docker,
      disk_space_gb = additional_disk
	}

	call CNNTasks.SamtoolsMergeBAMs {
    input:
      input_bams = RunHC4.bamout,
      output_prefix = output_prefix,
      disk_space_gb = round(bam_size + ref_size + additional_disk)
	}

	output {

    ## Unmapping
    Array[File] validatesam_out = ValidateSamFile.report

    ### CNN Variant Filter
    File FilterVariantTranches.*

    ### REMAPPING
    File duplication_metrics = MarkDuplicates.duplicate_metrics
    File bqsr_report = GatherBqsrReports.output_bqsr_report
    File analysis_ready_bam = GatherBamFiles.output_bam
    File analysis_ready_bam_index = GatherBamFiles.output_bam_index
    File analysis_ready_bam_md5 = GatherBamFiles.output_bam_md5

    ## GATK SOMOATIC CNV CALL
    File preprocessed_intervals = PreprocessIntervals.preprocessed_intervals

    File read_counts_entity_id_tumor = CollectCountsTumor.entity_id
    File read_counts_tumor = CollectCountsTumor.counts
    File allelic_counts_entity_id_tumor = CollectAllelicCountsTumor.entity_id
    File allelic_counts_tumor = CollectAllelicCountsTumor.allelic_counts
    File denoised_copy_ratios_tumor = DenoiseReadCountsTumor.denoised_copy_ratios
    File standardized_copy_ratios_tumor = DenoiseReadCountsTumor.standardized_copy_ratios
    File het_allelic_counts_tumor = ModelSegmentsTumor.het_allelic_counts
    File normal_het_allelic_counts_tumor = ModelSegmentsTumor.normal_het_allelic_counts
    File copy_ratio_only_segments_tumor = ModelSegmentsTumor.copy_ratio_only_segments
    File copy_ratio_legacy_segments_tumor = ModelSegmentsTumor.copy_ratio_legacy_segments
    File allele_fraction_legacy_segments_tumor = ModelSegmentsTumor.allele_fraction_legacy_segments
    File modeled_segments_begin_tumor = ModelSegmentsTumor.modeled_segments_begin
    File copy_ratio_parameters_begin_tumor = ModelSegmentsTumor.copy_ratio_parameters_begin
    File allele_fraction_parameters_begin_tumor = ModelSegmentsTumor.allele_fraction_parameters_begin
    File modeled_segments_tumor = ModelSegmentsTumor.modeled_segments
    File copy_ratio_parameters_tumor = ModelSegmentsTumor.copy_ratio_parameters
    File allele_fraction_parameters_tumor = ModelSegmentsTumor.allele_fraction_parameters
    File called_copy_ratio_segments_tumor = CallCopyRatioSegmentsTumor.called_copy_ratio_segments
    File called_copy_ratio_legacy_segments_tumor = CallCopyRatioSegmentsTumor.called_copy_ratio_legacy_segments
    File denoised_copy_ratios_plot_tumor = PlotDenoisedCopyRatiosTumor.denoised_copy_ratios_plot
    File denoised_copy_ratios_lim_4_plot_tumor = PlotDenoisedCopyRatiosTumor.denoised_copy_ratios_lim_4_plot
    File standardized_MAD_tumor = PlotDenoisedCopyRatiosTumor.standardized_MAD
    Float standardized_MAD_value_tumor = PlotDenoisedCopyRatiosTumor.standardized_MAD_value
    File denoised_MAD_tumor = PlotDenoisedCopyRatiosTumor.denoised_MAD
    Float denoised_MAD_value_tumor = PlotDenoisedCopyRatiosTumor.denoised_MAD_value
    File delta_MAD_tumor = PlotDenoisedCopyRatiosTumor.delta_MAD
    Float delta_MAD_value_tumor = PlotDenoisedCopyRatiosTumor.delta_MAD_value
    File scaled_delta_MAD_tumor = PlotDenoisedCopyRatiosTumor.scaled_delta_MAD
    Float scaled_delta_MAD_value_tumor = PlotDenoisedCopyRatiosTumor.scaled_delta_MAD_value
    File modeled_segments_plot_tumor = PlotModeledSegmentsTumor.modeled_segments_plot

    File? read_counts_entity_id_normal = CollectCountsNormal.entity_id
    File? read_counts_normal = CollectCountsNormal.counts
    File? allelic_counts_entity_id_normal = CollectAllelicCountsNormal.entity_id
    File? allelic_counts_normal = CollectAllelicCountsNormal.allelic_counts
    File? denoised_copy_ratios_normal = DenoiseReadCountsNormal.denoised_copy_ratios
    File? standardized_copy_ratios_normal = DenoiseReadCountsNormal.standardized_copy_ratios
    File? het_allelic_counts_normal = ModelSegmentsNormal.het_allelic_counts
    File? normal_het_allelic_counts_normal = ModelSegmentsNormal.normal_het_allelic_counts
    File? copy_ratio_only_segments_normal = ModelSegmentsNormal.copy_ratio_only_segments
    File? copy_ratio_legacy_segments_normal = ModelSegmentsNormal.copy_ratio_legacy_segments
    File? allele_fraction_legacy_segments_normal = ModelSegmentsNormal.allele_fraction_legacy_segments
    File? modeled_segments_begin_normal = ModelSegmentsNormal.modeled_segments_begin
    File? copy_ratio_parameters_begin_normal = ModelSegmentsNormal.copy_ratio_parameters_begin
    File? allele_fraction_parameters_begin_normal = ModelSegmentsNormal.allele_fraction_parameters_begin
    File? modeled_segments_normal = ModelSegmentsNormal.modeled_segments
    File? copy_ratio_parameters_normal = ModelSegmentsNormal.copy_ratio_parameters
    File? allele_fraction_parameters_normal = ModelSegmentsNormal.allele_fraction_parameters
    File? called_copy_ratio_segments_normal = CallCopyRatioSegmentsNormal.called_copy_ratio_segments
    File? called_copy_ratio_legacy_segments_normal = CallCopyRatioSegmentsNormal.called_copy_ratio_legacy_segments
    File? denoised_copy_ratios_plot_normal = PlotDenoisedCopyRatiosNormal.denoised_copy_ratios_plot
    File? denoised_copy_ratios_lim_4_plot_normal = PlotDenoisedCopyRatiosNormal.denoised_copy_ratios_lim_4_plot
    File? standardized_MAD_normal = PlotDenoisedCopyRatiosNormal.standardized_MAD
    Float? standardized_MAD_value_normal = PlotDenoisedCopyRatiosNormal.standardized_MAD_value
    File? denoised_MAD_normal = PlotDenoisedCopyRatiosNormal.denoised_MAD
    Float? denoised_MAD_value_normal = PlotDenoisedCopyRatiosNormal.denoised_MAD_value
    File? delta_MAD_normal = PlotDenoisedCopyRatiosNormal.delta_MAD
    Float? delta_MAD_value_normal = PlotDenoisedCopyRatiosNormal.delta_MAD_value
    File? scaled_delta_MAD_normal = PlotDenoisedCopyRatiosNormal.scaled_delta_MAD
    Float? scaled_delta_MAD_value_normal = PlotDenoisedCopyRatiosNormal.scaled_delta_MAD_value
    File? modeled_segments_plot_normal = PlotModeledSegmentsNormal.modeled_segments_plot

    File oncotated_called_file_tumor = select_first([CNVOncotatorWorkflow.oncotated_called_file, "null"])
    File oncotated_called_gene_list_file_tumor = select_first([CNVOncotatorWorkflow.oncotated_called_gene_list_file, "null"])

    #### MANTA
    File Manta.*
	}
}


# Create the txt file from the list of files
task CreateTxt {
	# Command parameters
	Array[String] array_of_files
	String list_name

	command {
    mv ${write_lines(array_of_files)} ${list_name}.txt
	}
	output {
    File file_list_name = "${list_name}.txt"
	}
	runtime {
    docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1500064817"
    preemptible: 3
	}
}
