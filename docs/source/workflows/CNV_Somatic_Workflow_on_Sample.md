##########

## CNVSomaticPairWorkflow

### Inputs

#### Required

  * `common_sites` (File, **required**)
  * `gatk_docker` (String, **required**)
  * `intervals` (File, **required**)
  * `read_count_pon` (File, **required**)
  * `ref_fasta` (File, **required**)
  * `ref_fasta_dict` (File, **required**)
  * `ref_fasta_fai` (File, **required**)
  * `tumor_bam` (File, **required**)
  * `tumor_bam_idx` (File, **required**)

#### Optional

  * `additional_args_for_oncotator` (String?)
  * `bin_length` (Int?)
  * `blacklist_intervals` (File?)
  * `boot_disk_space_gb_for_oncotator` (Int?)
  * `calling_copy_ratio_z_score_threshold` (Float?)
  * `collect_counts_format` (String?)
  * `emergency_extra_disk` (Int?)
  * `gatk4_jar_override` (File?)
  * `genotyping_base_error_rate` (Float?)
  * `genotyping_homozygous_log_ratio_threshold` (Float?)
  * `is_run_oncotator` (Boolean?)
  * `kernel_approximation_dimension` (Int?)
  * `kernel_scaling_allele_fraction` (Float?)
  * `kernel_variance_allele_fraction` (Float?)
  * `kernel_variance_copy_ratio` (Float?)
  * `max_num_segments_per_chromosome` (Int?)
  * `max_num_smoothing_iterations` (Int?)
  * `mem_gb_for_call_copy_ratio_segments` (Int?)
  * `mem_gb_for_collect_allelic_counts` (Int?)
  * `mem_gb_for_collect_counts` (Int?)
  * `mem_gb_for_denoise_read_counts` (Int?)
  * `mem_gb_for_model_segments` (Int?)
  * `mem_gb_for_oncotator` (Int?)
  * `mem_gb_for_plotting` (Int?)
  * `mem_gb_for_preprocess_intervals` (Int?)
  * `min_total_allele_count` (Int?)
  * `min_total_allele_count_normal` (Int?)
  * `minimum_base_quality` (String?)
  * `minimum_contig_length` (Int?)
  * `minor_allele_fraction_prior_alpha` (Float?)
  * `neutral_segment_copy_ratio_lower_bound` (Float?)
  * `neutral_segment_copy_ratio_upper_bound` (Float?)
  * `normal_bam` (File?)
  * `normal_bam_idx` (File?)
  * `num_burn_in_allele_fraction` (Int?)
  * `num_burn_in_copy_ratio` (Int?)
  * `num_changepoints_penalty_factor` (Float?)
  * `num_samples_allele_fraction` (Int?)
  * `num_samples_copy_ratio` (Int?)
  * `num_smoothing_iterations_per_fit` (Int?)
  * `number_of_eigensamples` (Int?)
  * `oncotator_docker` (String?)
  * `outlier_neutral_segment_copy_ratio_z_score_threshold` (Float?)
  * `padding` (Int?)
  * `preemptible_attempts` (Int?)
  * `smoothing_threshold_allele_fraction` (Float?)
  * `smoothing_threshold_copy_ratio` (Float?)
  * `window_sizes` (Array[Int]+?)
  * `CallCopyRatioSegmentsNormal.cpu` (Int?)
  * `CallCopyRatioSegmentsTumor.cpu` (Int?)
  * `CollectAllelicCountsNormal.cpu` (Int?)
  * `CollectAllelicCountsTumor.cpu` (Int?)
  * `CollectCountsNormal.cpu` (Int?)
  * `CollectCountsTumor.cpu` (Int?)
  * `DenoiseReadCountsNormal.cpu` (Int?)
  * `DenoiseReadCountsTumor.cpu` (Int?)
  * `ModelSegmentsNormal.cpu` (Int?)
  * `ModelSegmentsNormal.min_total_allele_count_normal` (Int?)
  * `ModelSegmentsNormal.normal_allelic_counts` (File?)
  * `ModelSegmentsNormal.output_dir` (String?)
  * `ModelSegmentsTumor.cpu` (Int?)
  * `ModelSegmentsTumor.output_dir` (String?)
  * `PlotDenoisedCopyRatiosNormal.cpu` (Int?)
  * `PlotDenoisedCopyRatiosNormal.output_dir` (String?)
  * `PlotDenoisedCopyRatiosTumor.cpu` (Int?)
  * `PlotDenoisedCopyRatiosTumor.output_dir` (String?)
  * `PlotModeledSegmentsNormal.cpu` (Int?)
  * `PlotModeledSegmentsNormal.output_dir` (String?)
  * `PlotModeledSegmentsTumor.cpu` (Int?)
  * `PlotModeledSegmentsTumor.output_dir` (String?)
  * `PreprocessIntervals.cpu` (Int?)
  * `CNVOncotatorWorkflow.OncotateSegments.cpu` (Int?)
  * `CNVOncotatorWorkflow.OncotateSegments.disk_space_gb` (Int?)

#### Defaults

  * `collect_allelic_counts_tumor_disk` (Int, default=tumor_bam_size + ref_size + disk_pad)
  * `collect_counts_normal_disk` (Int, default=normal_bam_size + ceil(size(PreprocessIntervals.preprocessed_intervals,"GB")) + disk_pad)
  * `collect_counts_tumor_disk` (Int, default=tumor_bam_size + ceil(size(PreprocessIntervals.preprocessed_intervals,"GB")) + disk_pad)
  * `copy_ratio_segments_tumor_disk` (Int, default=ceil(size(DenoiseReadCountsTumor.denoised_copy_ratios,"GB")) + ceil(size(ModelSegmentsTumor.copy_ratio_only_segments,"GB")) + disk_pad)
  * `denoise_read_counts_tumor_disk` (Int, default=read_count_pon_size + ceil(size(CollectCountsTumor.counts,"GB")) + disk_pad)
  * `disk_pad` (Int, default=20 + ceil(size(intervals,"GB")) + ceil(size(common_sites,"GB")) + gatk4_override_size + select_first([emergency_extra_disk, 0]))
  * `final_normal_bam` (File, default=select_first([normal_bam, "null"]))
  * `final_normal_bam_idx` (File, default=select_first([normal_bam_idx, "null"]))
  * `gatk4_override_size` (Int, default=if defined(gatk4_jar_override) then ceil(size(gatk4_jar_override,"GB")) else 0)
  * `model_segments_normal_portion` (Int, default=if defined(normal_bam) then ceil(size(CollectAllelicCountsNormal.allelic_counts,"GB")) else 0)
  * `model_segments_tumor_disk` (Int, default=ceil(size(DenoiseReadCountsTumor.denoised_copy_ratios,"GB")) + ceil(size(CollectAllelicCountsTumor.allelic_counts,"GB")) + model_segments_normal_portion + disk_pad)
  * `normal_bam_size` (Int, default=if defined(normal_bam) then ceil((size(normal_bam,"GB") + size(normal_bam_idx,"GB"))) else 0)
  * `plot_tumor_disk` (Int, default=ref_size + ceil(size(DenoiseReadCountsTumor.standardized_copy_ratios,"GB")) + ceil(size(DenoiseReadCountsTumor.denoised_copy_ratios,"GB")) + ceil(size(ModelSegmentsTumor.het_allelic_counts,"GB")) + ceil(size(ModelSegmentsTumor.modeled_segments,"GB")) + disk_pad)
  * `preprocess_intervals_disk` (Int, default=ref_size + disk_pad)
  * `read_count_pon_size` (Int, default=ceil(size(read_count_pon,"GB")))
  * `ref_size` (Int, default=ceil((size(ref_fasta,"GB") + size(ref_fasta_dict,"GB") + size(ref_fasta_fai,"GB"))))
  * `tumor_bam_size` (Int, default=ceil((size(tumor_bam,"GB") + size(tumor_bam_idx,"GB"))))
  * `CallCopyRatioSegmentsNormal.command_mem_mb` (Int, default=machine_mem_mb - 1000)
  * `CallCopyRatioSegmentsNormal.machine_mem_mb` (Int, default=select_first([mem_gb, 7]) * 1000)
  * `CallCopyRatioSegmentsNormal.use_ssd` (Boolean, default=false)
  * `CallCopyRatioSegmentsTumor.command_mem_mb` (Int, default=machine_mem_mb - 1000)
  * `CallCopyRatioSegmentsTumor.machine_mem_mb` (Int, default=select_first([mem_gb, 7]) * 1000)
  * `CallCopyRatioSegmentsTumor.use_ssd` (Boolean, default=false)
  * `CollectAllelicCountsNormal.allelic_counts_filename` (String, default="~{base_filename}.allelicCounts.tsv")
  * `CollectAllelicCountsNormal.base_filename` (String, default=basename(bam,".bam"))
  * `CollectAllelicCountsNormal.command_mem_mb` (Int, default=machine_mem_mb - 1000)
  * `CollectAllelicCountsNormal.machine_mem_mb` (Int, default=select_first([mem_gb, 13]) * 1000)
  * `CollectAllelicCountsNormal.use_ssd` (Boolean, default=false)
  * `CollectAllelicCountsTumor.allelic_counts_filename` (String, default="~{base_filename}.allelicCounts.tsv")
  * `CollectAllelicCountsTumor.base_filename` (String, default=basename(bam,".bam"))
  * `CollectAllelicCountsTumor.command_mem_mb` (Int, default=machine_mem_mb - 1000)
  * `CollectAllelicCountsTumor.machine_mem_mb` (Int, default=select_first([mem_gb, 13]) * 1000)
  * `CollectAllelicCountsTumor.use_ssd` (Boolean, default=false)
  * `CollectCountsNormal.base_filename` (String, default=basename(bam,".bam"))
  * `CollectCountsNormal.command_mem_mb` (Int, default=machine_mem_mb - 1000)
  * `CollectCountsNormal.counts_filename` (String, default=if !defined(format) then "~{base_filename}.counts.hdf5" else "~{base_filename}.counts.tsv")
  * `CollectCountsNormal.machine_mem_mb` (Int, default=select_first([mem_gb, 7]) * 1000)
  * `CollectCountsNormal.use_ssd` (Boolean, default=false)
  * `CollectCountsTumor.base_filename` (String, default=basename(bam,".bam"))
  * `CollectCountsTumor.command_mem_mb` (Int, default=machine_mem_mb - 1000)
  * `CollectCountsTumor.counts_filename` (String, default=if !defined(format) then "~{base_filename}.counts.hdf5" else "~{base_filename}.counts.tsv")
  * `CollectCountsTumor.machine_mem_mb` (Int, default=select_first([mem_gb, 7]) * 1000)
  * `CollectCountsTumor.use_ssd` (Boolean, default=false)
  * `DenoiseReadCountsNormal.command_mem_mb` (Int, default=machine_mem_mb - 1000)
  * `DenoiseReadCountsNormal.machine_mem_mb` (Int, default=select_first([mem_gb, 13]) * 1000)
  * `DenoiseReadCountsNormal.use_ssd` (Boolean, default=false)
  * `DenoiseReadCountsTumor.command_mem_mb` (Int, default=machine_mem_mb - 1000)
  * `DenoiseReadCountsTumor.machine_mem_mb` (Int, default=select_first([mem_gb, 13]) * 1000)
  * `DenoiseReadCountsTumor.use_ssd` (Boolean, default=false)
  * `ModelSegmentsNormal.command_mem_mb` (Int, default=machine_mem_mb - 3000)
  * `ModelSegmentsNormal.default_min_total_allele_count` (Int, default=if defined(normal_allelic_counts) then 0 else 30)
  * `ModelSegmentsNormal.machine_mem_mb` (Int, default=select_first([mem_gb, 13]) * 1000)
  * `ModelSegmentsNormal.min_total_allele_count_` (Int, default=select_first([min_total_allele_count, default_min_total_allele_count]))
  * `ModelSegmentsNormal.output_dir_` (String, default=select_first([output_dir, "out"]))
  * `ModelSegmentsNormal.use_ssd` (Boolean, default=false)
  * `ModelSegmentsTumor.command_mem_mb` (Int, default=machine_mem_mb - 3000)
  * `ModelSegmentsTumor.default_min_total_allele_count` (Int, default=if defined(normal_allelic_counts) then 0 else 30)
  * `ModelSegmentsTumor.machine_mem_mb` (Int, default=select_first([mem_gb, 13]) * 1000)
  * `ModelSegmentsTumor.min_total_allele_count_` (Int, default=select_first([min_total_allele_count, default_min_total_allele_count]))
  * `ModelSegmentsTumor.output_dir_` (String, default=select_first([output_dir, "out"]))
  * `ModelSegmentsTumor.use_ssd` (Boolean, default=false)
  * `PlotDenoisedCopyRatiosNormal.command_mem_mb` (Int, default=machine_mem_mb - 1000)
  * `PlotDenoisedCopyRatiosNormal.machine_mem_mb` (Int, default=select_first([mem_gb, 7]) * 1000)
  * `PlotDenoisedCopyRatiosNormal.output_dir_` (String, default=select_first([output_dir, "out"]))
  * `PlotDenoisedCopyRatiosNormal.use_ssd` (Boolean, default=false)
  * `PlotDenoisedCopyRatiosTumor.command_mem_mb` (Int, default=machine_mem_mb - 1000)
  * `PlotDenoisedCopyRatiosTumor.machine_mem_mb` (Int, default=select_first([mem_gb, 7]) * 1000)
  * `PlotDenoisedCopyRatiosTumor.output_dir_` (String, default=select_first([output_dir, "out"]))
  * `PlotDenoisedCopyRatiosTumor.use_ssd` (Boolean, default=false)
  * `PlotModeledSegmentsNormal.command_mem_mb` (Int, default=machine_mem_mb - 1000)
  * `PlotModeledSegmentsNormal.machine_mem_mb` (Int, default=select_first([mem_gb, 7]) * 1000)
  * `PlotModeledSegmentsNormal.output_dir_` (String, default=select_first([output_dir, "out"]))
  * `PlotModeledSegmentsNormal.use_ssd` (Boolean, default=false)
  * `PlotModeledSegmentsTumor.command_mem_mb` (Int, default=machine_mem_mb - 1000)
  * `PlotModeledSegmentsTumor.machine_mem_mb` (Int, default=select_first([mem_gb, 7]) * 1000)
  * `PlotModeledSegmentsTumor.output_dir_` (String, default=select_first([output_dir, "out"]))
  * `PlotModeledSegmentsTumor.use_ssd` (Boolean, default=false)
  * `PreprocessIntervals.base_filename` (String, default=basename(filename,".interval_list"))
  * `PreprocessIntervals.command_mem_mb` (Int, default=machine_mem_mb - 500)
  * `PreprocessIntervals.filename` (String, default=select_first([intervals, "wgs"]))
  * `PreprocessIntervals.machine_mem_mb` (Int, default=select_first([mem_gb, 2]) * 1000)
  * `PreprocessIntervals.use_ssd` (Boolean, default=false)
  * `CNVOncotatorWorkflow.OncotateSegments.basename_called_file` (String, default=basename(called_file))
  * `CNVOncotatorWorkflow.OncotateSegments.machine_mem_mb` (Int, default=select_first([mem_gb, 3]) * 1000)
  * `CNVOncotatorWorkflow.OncotateSegments.use_ssd` (Boolean, default=false)

### Outputs

  * `preprocessed_intervals` (File)
  * `read_counts_entity_id_tumor` (File)
  * `read_counts_tumor` (File)
  * `allelic_counts_entity_id_tumor` (File)
  * `allelic_counts_tumor` (File)
  * `denoised_copy_ratios_tumor` (File)
  * `standardized_copy_ratios_tumor` (File)
  * `het_allelic_counts_tumor` (File)
  * `normal_het_allelic_counts_tumor` (File)
  * `copy_ratio_only_segments_tumor` (File)
  * `copy_ratio_legacy_segments_tumor` (File)
  * `allele_fraction_legacy_segments_tumor` (File)
  * `modeled_segments_begin_tumor` (File)
  * `copy_ratio_parameters_begin_tumor` (File)
  * `allele_fraction_parameters_begin_tumor` (File)
  * `modeled_segments_tumor` (File)
  * `copy_ratio_parameters_tumor` (File)
  * `allele_fraction_parameters_tumor` (File)
  * `called_copy_ratio_segments_tumor` (File)
  * `called_copy_ratio_legacy_segments_tumor` (File)
  * `denoised_copy_ratios_plot_tumor` (File)
  * `denoised_copy_ratios_lim_4_plot_tumor` (File)
  * `standardized_MAD_tumor` (File)
  * `standardized_MAD_value_tumor` (Float)
  * `denoised_MAD_tumor` (File)
  * `denoised_MAD_value_tumor` (Float)
  * `delta_MAD_tumor` (File)
  * `delta_MAD_value_tumor` (Float)
  * `scaled_delta_MAD_tumor` (File)
  * `scaled_delta_MAD_value_tumor` (Float)
  * `modeled_segments_plot_tumor` (File)
  * `read_counts_entity_id_normal` (File?)
  * `read_counts_normal` (File?)
  * `allelic_counts_entity_id_normal` (File?)
  * `allelic_counts_normal` (File?)
  * `denoised_copy_ratios_normal` (File?)
  * `standardized_copy_ratios_normal` (File?)
  * `het_allelic_counts_normal` (File?)
  * `normal_het_allelic_counts_normal` (File?)
  * `copy_ratio_only_segments_normal` (File?)
  * `copy_ratio_legacy_segments_normal` (File?)
  * `allele_fraction_legacy_segments_normal` (File?)
  * `modeled_segments_begin_normal` (File?)
  * `copy_ratio_parameters_begin_normal` (File?)
  * `allele_fraction_parameters_begin_normal` (File?)
  * `modeled_segments_normal` (File?)
  * `copy_ratio_parameters_normal` (File?)
  * `allele_fraction_parameters_normal` (File?)
  * `called_copy_ratio_segments_normal` (File?)
  * `called_copy_ratio_legacy_segments_normal` (File?)
  * `denoised_copy_ratios_plot_normal` (File?)
  * `denoised_copy_ratios_lim_4_plot_normal` (File?)
  * `standardized_MAD_normal` (File?)
  * `standardized_MAD_value_normal` (Float?)
  * `denoised_MAD_normal` (File?)
  * `denoised_MAD_value_normal` (Float?)
  * `delta_MAD_normal` (File?)
  * `delta_MAD_value_normal` (Float?)
  * `scaled_delta_MAD_normal` (File?)
  * `scaled_delta_MAD_value_normal` (Float?)
  * `modeled_segments_plot_normal` (File?)
  * `oncotated_called_file_tumor` (File)
  * `oncotated_called_gene_list_file_tumor` (File)

## DenoiseReadCounts

### Inputs

#### Required

  * `entity_id` (String, **required**)
  * `gatk_docker` (String, **required**)
  * `read_count_pon` (File, **required**)
  * `read_counts` (File, **required**)

#### Optional

  * `cpu` (Int?)
  * `disk_space_gb` (Int?)
  * `gatk4_jar_override` (File?)
  * `mem_gb` (Int?)
  * `number_of_eigensamples` (Int?)
  * `preemptible_attempts` (Int?)

#### Defaults

  * `command_mem_mb` (Int, default=machine_mem_mb - 1000)
  * `machine_mem_mb` (Int, default=select_first([mem_gb, 13]) * 1000)
  * `use_ssd` (Boolean, default=false)

### Outputs

  * `standardized_copy_ratios` (File)
  * `denoised_copy_ratios` (File)

## ModelSegments

### Inputs

#### Required

  * `allelic_counts` (File, **required**)
  * `denoised_copy_ratios` (File, **required**)
  * `entity_id` (String, **required**)
  * `gatk_docker` (String, **required**)

#### Optional

  * `cpu` (Int?)
  * `disk_space_gb` (Int?)
  * `gatk4_jar_override` (File?)
  * `genotyping_base_error_rate` (Float?)
  * `genotyping_homozygous_log_ratio_threshold` (Float?)
  * `kernel_approximation_dimension` (Int?)
  * `kernel_scaling_allele_fraction` (Float?)
  * `kernel_variance_allele_fraction` (Float?)
  * `kernel_variance_copy_ratio` (Float?)
  * `max_num_segments_per_chromosome` (Int?)
  * `max_num_smoothing_iterations` (Int?)
  * `mem_gb` (Int?)
  * `min_total_allele_count` (Int?)
  * `min_total_allele_count_normal` (Int?)
  * `minor_allele_fraction_prior_alpha` (Float?)
  * `normal_allelic_counts` (File?)
  * `num_burn_in_allele_fraction` (Int?)
  * `num_burn_in_copy_ratio` (Int?)
  * `num_changepoints_penalty_factor` (Float?)
  * `num_samples_allele_fraction` (Int?)
  * `num_samples_copy_ratio` (Int?)
  * `num_smoothing_iterations_per_fit` (Int?)
  * `output_dir` (String?)
  * `preemptible_attempts` (Int?)
  * `smoothing_threshold_allele_fraction` (Float?)
  * `smoothing_threshold_copy_ratio` (Float?)
  * `window_sizes` (Array[Int]+?)

#### Defaults

  * `command_mem_mb` (Int, default=machine_mem_mb - 3000)
  * `default_min_total_allele_count` (Int, default=if defined(normal_allelic_counts) then 0 else 30)
  * `machine_mem_mb` (Int, default=select_first([mem_gb, 13]) * 1000)
  * `min_total_allele_count_` (Int, default=select_first([min_total_allele_count, default_min_total_allele_count]))
  * `output_dir_` (String, default=select_first([output_dir, "out"]))
  * `use_ssd` (Boolean, default=false)

### Outputs

  * `het_allelic_counts` (File)
  * `normal_het_allelic_counts` (File)
  * `copy_ratio_only_segments` (File)
  * `copy_ratio_legacy_segments` (File)
  * `allele_fraction_legacy_segments` (File)
  * `modeled_segments_begin` (File)
  * `copy_ratio_parameters_begin` (File)
  * `allele_fraction_parameters_begin` (File)
  * `modeled_segments` (File)
  * `copy_ratio_parameters` (File)
  * `allele_fraction_parameters` (File)

## CallCopyRatioSegments

### Inputs

#### Required

  * `copy_ratio_segments` (File, **required**)
  * `entity_id` (String, **required**)
  * `gatk_docker` (String, **required**)

#### Optional

  * `calling_copy_ratio_z_score_threshold` (Float?)
  * `cpu` (Int?)
  * `disk_space_gb` (Int?)
  * `gatk4_jar_override` (File?)
  * `mem_gb` (Int?)
  * `neutral_segment_copy_ratio_lower_bound` (Float?)
  * `neutral_segment_copy_ratio_upper_bound` (Float?)
  * `outlier_neutral_segment_copy_ratio_z_score_threshold` (Float?)
  * `preemptible_attempts` (Int?)

#### Defaults

  * `command_mem_mb` (Int, default=machine_mem_mb - 1000)
  * `machine_mem_mb` (Int, default=select_first([mem_gb, 7]) * 1000)
  * `use_ssd` (Boolean, default=false)

### Outputs

  * `called_copy_ratio_segments` (File)
  * `called_copy_ratio_legacy_segments` (File)

## PlotDenoisedCopyRatios

### Inputs

#### Required

  * `denoised_copy_ratios` (File, **required**)
  * `entity_id` (String, **required**)
  * `gatk_docker` (String, **required**)
  * `ref_fasta_dict` (File, **required**)
  * `standardized_copy_ratios` (File, **required**)

#### Optional

  * `cpu` (Int?)
  * `disk_space_gb` (Int?)
  * `gatk4_jar_override` (File?)
  * `mem_gb` (Int?)
  * `minimum_contig_length` (Int?)
  * `output_dir` (String?)
  * `preemptible_attempts` (Int?)

#### Defaults

  * `command_mem_mb` (Int, default=machine_mem_mb - 1000)
  * `machine_mem_mb` (Int, default=select_first([mem_gb, 7]) * 1000)
  * `output_dir_` (String, default=select_first([output_dir, "out"]))
  * `use_ssd` (Boolean, default=false)

### Outputs

  * `denoised_copy_ratios_plot` (File)
  * `denoised_copy_ratios_lim_4_plot` (File)
  * `standardized_MAD` (File)
  * `standardized_MAD_value` (Float)
  * `denoised_MAD` (File)
  * `denoised_MAD_value` (Float)
  * `delta_MAD` (File)
  * `delta_MAD_value` (Float)
  * `scaled_delta_MAD` (File)
  * `scaled_delta_MAD_value` (Float)

## PlotModeledSegments

### Inputs

#### Required

  * `denoised_copy_ratios` (File, **required**)
  * `entity_id` (String, **required**)
  * `gatk_docker` (String, **required**)
  * `het_allelic_counts` (File, **required**)
  * `modeled_segments` (File, **required**)
  * `ref_fasta_dict` (File, **required**)

#### Optional

  * `cpu` (Int?)
  * `disk_space_gb` (Int?)
  * `gatk4_jar_override` (File?)
  * `mem_gb` (Int?)
  * `minimum_contig_length` (Int?)
  * `output_dir` (String?)
  * `preemptible_attempts` (Int?)

#### Defaults

  * `command_mem_mb` (Int, default=machine_mem_mb - 1000)
  * `machine_mem_mb` (Int, default=select_first([mem_gb, 7]) * 1000)
  * `output_dir_` (String, default=select_first([output_dir, "out"]))
  * `use_ssd` (Boolean, default=false)

### Outputs

  * `modeled_segments_plot` (File)
