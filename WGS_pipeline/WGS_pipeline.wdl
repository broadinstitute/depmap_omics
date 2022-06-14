version 1.0

import "BamToUnmappedRGBams.wdl" as BamToUnmappedRGBams
import "PreProcessingForVariantDiscovery_GATK4.wdl" as PreProcessingForVariantDiscovery_GATK4
import "CNV_Somatic_Workflow_on_Sample.wdl" as CNV_Somatic_Workflow_on_Sample
import "Manta_SomaticSV.wdl" as Manta_SomaticSV
import "ArrayOfFilesToTxt.wdl" as ArrayOfFilesToTxt
import "gatk_mutect2_v21.wdl" as mutect2
import "bcftools.wdl" as setGT
import "fix_mutect2.wdl" as fixmutect2
import "remove_filtered.wdl" as removeFiltered
import "vcf_to_depmap.wdl" as vcf_to_depmap
import "PureCN.wdl" as PureCN
# import "opencravat.wdl" as openCravat

workflow WGS_pipeline {

    input {
        #BamToUnmappedRGBams
        #"broadinstitute/genomes-in-the-cloud:2.3.1-1504795437"
        File ref_fasta
        File ref_fasta_index
        File ref_dict

        File input_bam
        File input_bam_index

        String picard_path
        String picard_docker

        Int preemptible_tries

        String sample_name #this.name

        #PreProcessingForVariantDiscovery_GATK4
        #"broadinstitute/gatk:4.beta.3"
        #"broadinstitute/genomes-in-the-cloud:2.3.0-1501082129"
        String ref_name


        File dbSNP_vcf
        File dbSNP_vcf_index

        Array[File] known_indels_sites_VCFs
        Array[File] known_indels_sites_indices

        String gotc_docker
        String gatk_docker="broadinstitute/gatk:4.2.6.0"
        String python_docker

        String gcs_project_for_requester_pays

        String gotc_path

        #CNV_Somatic_Workflow_on_Sample
        #us.gcr.io/broad-gatk/gatk:4.1.5.0
        File common_sites
        File intervals

        # mutect2
        Int M2scatter=10

        File gnomad="gs://gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz"
        File gnomad_idx="gs://gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi"
        String m2_extra_args="--genotype-germline-sites true --genotype-pon-sites true"
        String? m2_filter_args
        File pon="gs://gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz"
        File pon_idx="gs://gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz.tbi"

        # PureCN
        File purecn_intervals
        File call_wgd_and_cin_script

    }

    call BamToUnmappedRGBams.BamToUnmappedRGBamsWf as BamToUnmappedRGBamsWf {
        input:
            ref_fasta=ref_fasta,
            ref_fasta_index=ref_fasta_index,
            input_bam=input_bam,
            picard_path=picard_path,
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
            ref_dict=ref_dict,
            ref_fasta=ref_fasta,
            ref_fasta_index=ref_fasta_index,
            dbSNP_vcf=dbSNP_vcf,
            dbSNP_vcf_index=dbSNP_vcf_index,
            known_indels_sites_VCFs=known_indels_sites_VCFs,
            known_indels_sites_indices=known_indels_sites_indices,
            gotc_docker=gotc_docker,
            python_docker=python_docker,
            gotc_path=gotc_path,
            picard_path=picard_path,
            preemptible_tries=preemptible_tries
    }

    call CNV_Somatic_Workflow_on_Sample.CNVSomaticPairWorkflow as CNVSomaticPairWorkflow {
        input:
            common_sites=common_sites,
            intervals=intervals,
            ref_fasta=ref_fasta,
            ref_fasta_dict=ref_dict,
            ref_fasta_fai=ref_fasta_index,
            tumor_bam=PreProcessingForVariantDiscovery_GATK4.analysis_ready_bam,
            tumor_bam_idx=PreProcessingForVariantDiscovery_GATK4.analysis_ready_bam_index
    }

    call Manta_SomaticSV.MantaSomaticSV as MantaSomaticSV {
        input:
            sample_name=sample_name,
            tumor_bam=PreProcessingForVariantDiscovery_GATK4.analysis_ready_bam,
            tumor_bam_index=PreProcessingForVariantDiscovery_GATK4.analysis_ready_bam_index,
            ref_fasta=ref_fasta,
            ref_fasta_index=ref_fasta_index,
            is_cram=false
    }

    call mutect2.Mutect2 as mutect2 {
        input:
            gatk_docker=gatk_docker,
            ref_dict=ref_dict,
            ref_fai=ref_fasta_index,
            ref_fasta=ref_fasta,
            scatter_count=M2scatter,
            tumor_reads=PreProcessingForVariantDiscovery_GATK4.analysis_ready_bam,
            tumor_reads_index=PreProcessingForVariantDiscovery_GATK4.analysis_ready_bam_index,
            intervals=intervals,
            gcs_project_for_requester_pays=gcs_project_for_requester_pays,
            compress_vcfs=true,
            filter_funcotations=false,
            funco_compress=true,
            funco_filter_funcotations=false,
            funco_output_format="VCF",
            funco_reference_version="hg38",
            gnomad=gnomad,
            gnomad_idx=gnomad_idx,
            m2_extra_args=m2_extra_args,
            m2_extra_filtering_args=m2_filter_args,
            make_bamout=false,
            pon=pon,
            pon_idx=pon_idx,
            run_funcotator=true,
            run_orientation_bias_mixture_model_filter=true
    }

    call PureCN.PureCN as PureCN {
        input:
            sample_id=sample_name,
            segFile=CNVSomaticPairWorkflow.modeled_segments_tumor,
            vcf=mutect2.funcotated_output_file,
            intervals=purecn_intervals,
            call_wgd_and_cin_script=call_wgd_and_cin_script,
    } 

    call setGT.bcftools_fix_ploidy as set_GT {
        input:
            sample_id=sample_name,
            vcf=select_first([mutect2.funcotated_file, mutect2.funcotated_output_file]),
    }

    call fixmutect2.fix_mutect2 as fix_mutect2 {
        input:
            sample_id=sample_name,
            vcf_file=set_GT.vcf_fixedploid
    }

    call removeFiltered.RemoveFiltered as RemoveFiltered {
        input:
            sample_id=sample_name,
            input_vcf=fix_mutect2.vcf_fixed
    }

    call vcf_to_depmap.vcf_to_depmap as my_vcf_to_depmap {
        input:
            input_vcf=RemoveFiltered.output_vcf,
            sample_id=sample_name,
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
        File oncotated_called_file_tumor = CNVSomaticPairWorkflow.oncotated_called_file_tumor
        File oncotated_called_gene_list_file_tumor = CNVSomaticPairWorkflow.oncotated_called_gene_list_file_tumor
        #MantaSomaticSV
        File candidate_indel_vcf= MantaSomaticSV.candidate_indel_vcf
        File candidate_indel_vcf_index= MantaSomaticSV.candidate_indel_vcf_index
        File candidate_sv_vcf= MantaSomaticSV.candidate_sv_vcf
        File candidate_sv_vcf_index= MantaSomaticSV.candidate_sv_vcf_index
        File germline_sv_vcf= MantaSomaticSV.germline_sv_vcf
        File germline_sv_vcf_index= MantaSomaticSV.germline_sv_vcf_index
        File somatic_sv_vcf= MantaSomaticSV.somatic_sv_vcf
        File somatic_sv_vcf_index= MantaSomaticSV.somatic_sv_vcf_index
        # omics_mutect2
        File omics_mutect2_out_vcf=fix_mutect2.vcf_fixed
        # PureCN
        File PureCN_solutions_pdf = PureCN.solutions_pdf
        File chromosomes_pdf = PureCN.chromosomes_pdf
        File PureCN_rds = PureCN.rds
        File PureCN_dnacopy = PureCN.dnacopy
        File PureCN_variants = PureCN.variants
        File PureCN_loh = PureCN.loh
        File PureCN_genes = PureCN.genes
        File PureCN_segmentation = PureCN.segmentation
        File PureCN_log = PureCN.log
        File PureCN_selected_solution = PureCN.selected_solution
        File PureCN_local_optima_pdf = PureCN.local_optima_pdf
        String PureCN_purity = PureCN.purity
        String PureCN_ploidy = PureCN.ploidy
        String PureCN_contamination = PureCN.contamination
        String PureCN_flagged = PureCN.flagged
        String PureCN_curated = PureCN.curated
        String PureCN_comment = PureCN.comment
        String PureCN_wgd = PureCN.wgd
        String PureCN_loh_fraction = PureCN.loh_fraction
        String PureCN_cin = PureCN.cin
        String PureCN_cin_allele_specific = PureCN.cin_allele_specific
        String PureCN_cin_ploidy_robust = PureCN.cin_ploidy_robust
        String PureCN_cin_allele_specific_ploidy_robust = PureCN.cin_allele_specific_ploidy_robust
        # vcf_to_depmap
        Array[File] main_output=my_vcf_to_depmap.full_file
        File somatic_maf=my_vcf_to_depmap.depmap_maf
    }
}
