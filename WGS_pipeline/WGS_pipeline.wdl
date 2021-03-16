import "https://raw.githubusercontent.com/broadinstitute/ccle_processing/reorganization/WGS_pipeline/BamToUnmappedRGBams.wdl" as BamToUnmappedRGBams 
import "https://raw.githubusercontent.com/broadinstitute/ccle_processing/reorganization/WGS_pipeline/ArrayOfFilesToTxt.wdl" as ArrayOfFilesToTxt 
import "https://raw.githubusercontent.com/broadinstitute/ccle_processing/reorganization/WGS_pipeline/PreProcessingForVariantDiscovery_GATK4.wdl" as PreProcessingForVariantDiscovery_GATK4 
import "https://raw.githubusercontent.com/broadinstitute/ccle_processing/reorganization/WGS_pipeline/CNV_Somatic_Workflow_on_Sample.wdl" as CNV_Somatic_Workflow_on_Sample 

workflow WGS_pipeline {

	#BamToUnmappedRGBams
	File ref_fasta
	File ref_fasta_index

	File input_bam

	String picard_path
	String picard_docker

	Int preemptible_tries

	#ArrayOfFilesToTxt
	String sample_name #this.name

	#PreProcessingForVariantDiscovery_GATK4
	String ref_name
	
	String unmapped_bam_suffix

	File ref_dict

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

	#CNVSomaticPairWorkflow
	File common_sites
	File intervals
	File read_count_pon

	call BamToUnmappedRGBams.BamToUnmappedRGBamsWf as BamToUnmappedRGBamsWf {
		input:
			ref_fasta=ref_fasta,
			ref_fasta_index=ref_fasta_index,
			input_bam=input_bam,
			picard_path=picard_path,
			picard_docker=picard_docker,
			preemptible_tries=preemptible_tries 
	}

	call ArrayOfFilesToTxt.CreateTxt as CreateTxt {
		input:
			array_of_files=BamToUnmappedRGBamsWf.sortsam_out,
			list_name=sample_name
	}

	call PreProcessingForVariantDiscovery_GATK4.PreProcessingForVariantDiscovery_GATK4 as PreProcessingForVariantDiscovery_GATK4 {
		input:
			sample_name=sample_name,
			ref_name=ref_name,
			flowcell_unmapped_bams_list=CreateTxt.file_list_name,
			unmapped_bam_suffix=unmapped_bam_suffix,
			ref_dict=ref_dict,
			ref_fasta=ref_fasta,
			ref_fasta_index=ref_fasta_index,
			bwa_commandline=bwa_commandline,
			compression_level=compression_level,
			dbSNP_vcf=dbSNP_vcf,
			dbSNP_vcf_index=dbSNP_vcf_index,
			known_indels_sites_VCFs=known_indels_sites_VCFs,
			known_indels_sites_indices=known_indels_sites_indices,
			gotc_docker=gotc_docker,
			picard_docker=picard_docker,
			gatk_docker=gatk_docker,
			python_docker=python_docker,
			gotc_path=gotc_path,
			picard_path=picard_path,
			gatk_launch_path=gatk_launch_path,
			flowcell_small_disk=flowcell_small_disk,
			flowcell_medium_disk=flowcell_medium_disk,
			agg_small_disk=agg_small_disk,
			agg_medium_disk=agg_medium_disk,
			agg_large_disk=agg_large_disk,
			preemptible_tries=preemptible_tries,
			agg_preemptible_tries=agg_preemptible_tries,
	
			# SamToFastqAndBwaMem.java_opt
			# SamToFastqAndBwaMem.mem_size
			# SamToFastqAndBwaMem.num_cpu
			# SamToFastqAndBwaMem.ref_amb
			# SamToFastqAndBwaMem.ref_ann
			# SamToFastqAndBwaMem.ref_bwt
			# SamToFastqAndBwaMem.ref_pac
			# SamToFastqAndBwaMem.ref_sa
			# SamToFastqAndBwaMem.ref_alt

			# SortAndFixTags.java_opt_fix
			# SortAndFixTags.java_opt_sort
			# SortAndFixTags.mem_size

			# ApplyBQSR.java_opt
			# ApplyBQSR.mem_size

			# BaseRecalibrator.java_opt
			# BaseRecalibrator.mem_size

			# CreateSequenceGroupingTSV.mem_size

			# GatherBamFiles.java_opt
			# GatherBamFiles.mem_size

			# GatherBqsrReports.java_opt
			# GatherBqsrReports.mem_size

			# GetBwaVersion.mem_size

			# MarkDuplicates.java_opt
			# MarkDuplicates.mem_size

			# MergeBamAlignment.java_opt
			# MergeBamAlignment.mem_size
	}

	call CNV_Somatic_Workflow_on_Sample.CNVSomaticPairWorkflow as CNVSomaticPairWorkflow {
		input:
			common_sites=common_sites,
			gatk_docker=gatk_docker,
			intervals=intervals,
			read_count_pon=read_count_pon,
			ref_fasta=ref_fasta,
			ref_fasta_dict=ref_dict,
			ref_fasta_fai=ref_fasta_index,
			tumor_bam=PreProcessingForVariantDiscovery_GATK4.analysis_ready_bam,
			tumor_bam_idx=PreProcessingForVariantDiscovery_GATK4.analysis_ready_bam_index
	}

	output {
		#BamToUnmappedRGBams
		#createTxT
		#PreProcessingForVariantDiscovery_GATK4
		File duplication_metrics = PreProcessingForVariantDiscovery_GATK4.duplication_metrics
		File bqsr_report = PreProcessingForVariantDiscovery_GATK4.bqsr_report
		File analysis_ready_bam = PreProcessingForVariantDiscovery_GATK4.analysis_ready_bam
		File analysis_ready_bam_index = PreProcessingForVariantDiscovery_GATK4.analysis_ready_bam_index
		File analysis_ready_bam_md5 = PreProcessingForVariantDiscovery_GATK4.analysis_ready_bam_md5
		#CNVSomaticPairWorkflow
		File preprocessed_intervals = CNVSomaticPairWorkflow.preprocessed_intervals
		File read_counts_entity_id_tumor = CNVSomaticPairWorkflow.read_counts_entity_id_tumor
		File read_counts_tumor = CNVSomaticPairWorkflow.read_counts_tumor
		File allelic_counts_entity_id_tumor = CNVSomaticPairWorkflow.allelic_counts_entity_id_tumor
		File allelic_counts_tumor = CNVSomaticPairWorkflow.allelic_counts_tumor
		File denoised_copy_ratios_tumor = CNVSomaticPairWorkflow.denoised_copy_ratios_tumor
		File standardized_copy_ratios_tumor = CNVSomaticPairWorkflow.standardized_copy_ratios_tumor
		File het_allelic_counts_tumor = CNVSomaticPairWorkflow.het_allelic_counts_tumor
		File normal_het_allelic_counts_tumor = CNVSomaticPairWorkflow.normal_het_allelic_counts_tumor
		File copy_ratio_only_segments_tumor = CNVSomaticPairWorkflow.copy_ratio_only_segments_tumor
		File copy_ratio_legacy_segments_tumor = CNVSomaticPairWorkflow.copy_ratio_legacy_segments_tumor
		File allele_fraction_legacy_segments_tumor = CNVSomaticPairWorkflow.allele_fraction_legacy_segments_tumor
		File modeled_segments_begin_tumor = CNVSomaticPairWorkflow.modeled_segments_begin_tumor
		File copy_ratio_parameters_begin_tumor = CNVSomaticPairWorkflow.copy_ratio_parameters_begin_tumor
		File allele_fraction_parameters_begin_tumor = CNVSomaticPairWorkflow.allele_fraction_parameters_begin_tumor
		File modeled_segments_tumor = CNVSomaticPairWorkflow.modeled_segments_tumor
		File copy_ratio_parameters_tumor = CNVSomaticPairWorkflow.copy_ratio_parameters_tumor
		File allele_fraction_parameters_tumor = CNVSomaticPairWorkflow.allele_fraction_parameters_tumor
		File called_copy_ratio_segments_tumor = CNVSomaticPairWorkflow.called_copy_ratio_segments_tumor
		File called_copy_ratio_legacy_segments_tumor = CNVSomaticPairWorkflow.called_copy_ratio_legacy_segments_tumor
		File denoised_copy_ratios_plot_tumor = CNVSomaticPairWorkflow.denoised_copy_ratios_plot_tumor
		File denoised_copy_ratios_lim_4_plot_tumor = CNVSomaticPairWorkflow.denoised_copy_ratios_lim_4_plot_tumor
		File standardized_MAD_tumor = CNVSomaticPairWorkflow.standardized_MAD_tumor
		Float standardized_MAD_value_tumor = CNVSomaticPairWorkflow.standardized_MAD_value_tumor
		File denoised_MAD_tumor = CNVSomaticPairWorkflow.denoised_MAD_tumor
		Float denoised_MAD_value_tumor = CNVSomaticPairWorkflow.denoised_MAD_value_tumor
		File delta_MAD_tumor = CNVSomaticPairWorkflow.delta_MAD_tumor
		Float delta_MAD_value_tumor = CNVSomaticPairWorkflow.delta_MAD_value_tumor
		File scaled_delta_MAD_tumor = CNVSomaticPairWorkflow.scaled_delta_MAD_tumor
		Float scaled_delta_MAD_value_tumor = CNVSomaticPairWorkflow.scaled_delta_MAD_value_tumor
		File modeled_segments_plot_tumor = CNVSomaticPairWorkflow.modeled_segments_plot_tumor
		File? read_counts_entity_id_normal = CNVSomaticPairWorkflow.read_counts_entity_id_normal
		File? read_counts_normal = CNVSomaticPairWorkflow.read_counts_normal
		File? allelic_counts_entity_id_normal = CNVSomaticPairWorkflow.allelic_counts_entity_id_normal
		File? allelic_counts_normal = CNVSomaticPairWorkflow.allelic_counts_normal
		File? denoised_copy_ratios_normal = CNVSomaticPairWorkflow.denoised_copy_ratios_normal
		File? standardized_copy_ratios_normal = CNVSomaticPairWorkflow.standardized_copy_ratios_normal
		File? het_allelic_counts_normal = CNVSomaticPairWorkflow.het_allelic_counts_normal
		File? normal_het_allelic_counts_normal = CNVSomaticPairWorkflow.normal_het_allelic_counts_normal
		File? copy_ratio_only_segments_normal = CNVSomaticPairWorkflow.copy_ratio_only_segments_normal
		File? copy_ratio_legacy_segments_normal = CNVSomaticPairWorkflow.copy_ratio_legacy_segments_normal
		File? allele_fraction_legacy_segments_normal = CNVSomaticPairWorkflow.allele_fraction_legacy_segments_normal
		File? modeled_segments_begin_normal = CNVSomaticPairWorkflow.modeled_segments_begin_normal
		File? copy_ratio_parameters_begin_normal = CNVSomaticPairWorkflow.copy_ratio_parameters_begin_normal
		File? allele_fraction_parameters_begin_normal = CNVSomaticPairWorkflow.allele_fraction_parameters_begin_normal
		File? modeled_segments_normal = CNVSomaticPairWorkflow.modeled_segments_normal
		File? copy_ratio_parameters_normal = CNVSomaticPairWorkflow.copy_ratio_parameters_normal
		File? allele_fraction_parameters_normal = CNVSomaticPairWorkflow.allele_fraction_parameters_normal
		File? called_copy_ratio_segments_normal = CNVSomaticPairWorkflow.called_copy_ratio_segments_normal
		File? called_copy_ratio_legacy_segments_normal = CNVSomaticPairWorkflow.called_copy_ratio_legacy_segments_normal
		File? denoised_copy_ratios_plot_normal = CNVSomaticPairWorkflow.denoised_copy_ratios_plot_normal
		File? denoised_copy_ratios_lim_4_plot_normal = CNVSomaticPairWorkflow.denoised_copy_ratios_lim_4_plot_normal
		File? standardized_MAD_normal = CNVSomaticPairWorkflow.standardized_MAD_normal
		Float? standardized_MAD_value_normal = CNVSomaticPairWorkflow.standardized_MAD_value_normal
		File? denoised_MAD_normal = CNVSomaticPairWorkflow.denoised_MAD_normal
		Float? denoised_MAD_value_normal = CNVSomaticPairWorkflow.denoised_MAD_value_normal
		File? delta_MAD_normal = CNVSomaticPairWorkflow.delta_MAD_normal
		Float? delta_MAD_value_normal = CNVSomaticPairWorkflow.delta_MAD_value_normal
		File? scaled_delta_MAD_normal = CNVSomaticPairWorkflow.scaled_delta_MAD_normal
		Float? scaled_delta_MAD_value_normal = CNVSomaticPairWorkflow.scaled_delta_MAD_value_normal
		File? modeled_segments_plot_normal = CNVSomaticPairWorkflow.modeled_segments_plot_normal
		File oncotated_called_file_tumor = CNVSomaticPairWorkflow.oncotated_called_file_tumor
		File oncotated_called_gene_list_file_tumor = CNVSomaticPairWorkflow.oncotated_called_gene_list_file_tumor
	}
}