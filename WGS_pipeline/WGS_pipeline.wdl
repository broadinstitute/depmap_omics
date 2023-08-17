version 1.0

import "cnv_somatic_pair_workflow.wdl" as CNV_Somatic_Workflow_on_Sample
import "Manta_SomaticSV_v1_0.wdl" as Manta_SomaticSV
import "manta_annot.wdl" as manta_annot
import "mutect2_v4.2.6.1.wdl" as mutect2
import "bcftools.wdl" as setGT
import "fix_mutect2.wdl" as fixmutect2
import "remove_filtered.wdl" as removeFiltered
import "../data/23Q2/WGSconfig/WGS_pipeline/vcf_to_depmap.wdl" as vcf_to_depmap
import "PureCN_pipeline/PureCN.wdl" as PureCN
import "msisensor2.wdl" as msisensor2
import "opencravat_dm.wdl" as openCravat

workflow WGS_pipeline {

    input {
        #BamToUnmappedRGBams
        #"broadinstitute/genomes-in-the-cloud:2.3.1-1504795437"
        File ref_fasta = "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta"
        File ref_fasta_index = "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai"
        File ref_dict = "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dict"

        File input_bam
        File input_bam_index

        String picard_path = "/usr/gitc/"
        String picard_docker = "broadinstitute/genomes-in-the-cloud:2.3.1-1504795437"

        Int preemptible_tries = 0

        String sample_name #this.name

        #PreProcessingForVariantDiscovery_GATK4
        #"broadinstitute/gatk:4.beta.3"
        #"broadinstitute/genomes-in-the-cloud:2.3.0-1501082129"
        String ref_name = "hg38"

        String gatk_docker="broadinstitute/gatk:4.2.6.1"
        # To be updated after 22Q4:
        String gatk_docker_cnv="us.gcr.io/broad-gatk/gatk:4.1.5.0"

        String gcs_project_for_requester_pays

        #CNV_Somatic_Workflow_on_Sample
        #us.gcr.io/broad-gatk/gatk:4.1.5.0
        File common_sites = "gs://ccleparams/references/intervals/common_sites_hg38_lifted.list"
        File intervals
        File read_count_pon
        Boolean is_run_funcotator_for_cnv
        Int cnv_preemptible_attempts = 1
        Float num_changepoints_penalty_factor = 5
        String funcotator_ref_version = "hg38"
        File funcotator_data_sources_tar_gz = "gs://broad-public-datasets/funcotator/funcotator_dataSources.v1.7.20200521s.tar.gz" # same is used in mutect2
        File? blacklist_intervals

        String manta_docker = "wwliao/manta:latest"
        String config_manta="/opt/conda/pkgs/manta-1.2.1-py27_0/bin/configManta.py"

        # mutect2
        Int M2scatter=30

        File? mutect2_intervals
        File? mutect2_intervals_for_contamination
        File gnomad="gs://gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz"
        File gnomad_idx="gs://gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi"
        String m2_extra_args="--genotype-germline-sites true --genotype-pon-sites true"
        String? m2_filter_args
        File pon="gs://gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz"
        File pon_idx="gs://gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz.tbi"
        String bcftools_exclude_string='FILTER~"weak_evidence" || FILTER~"map_qual" || FILTER~"strand_bias" || FILTER~"slippage" || FILTER~"clustered_events" || FILTER~"base_qual"'

        # PureCN
        File purecn_intervals = "gs://ccleparams/references/PureCN_intervals/wgs_hg38_2_percent_intervals.txt"

        #opencravat
        Array[String] annotators = ["cscape", "civic", "brca1_func_assay", "provean", "dann", "revel", "spliceai", "gtex", "funseq2", "pharmgkb", "dida", "gwas_catalog", "ccre_screen", "alfa"]
        File oc_modules = "gs://ccleparams/oc_modules.tar"
    }

    call CNV_Somatic_Workflow_on_Sample.CNVSomaticPairWorkflow as CNVSomaticPairWorkflow {
        input:
            common_sites=common_sites,
            intervals=intervals,
            ref_fasta=ref_fasta,
            ref_fasta_dict=ref_dict,
            ref_fasta_fai=ref_fasta_index,
            tumor_bam=input_bam,
            tumor_bam_idx=input_bam_index,
            read_count_pon=read_count_pon,
            gatk_docker=gatk_docker_cnv,
            is_run_funcotator=is_run_funcotator_for_cnv,
            funcotator_ref_version=funcotator_ref_version,
            gcs_project_for_requester_pays=gcs_project_for_requester_pays,
            funcotator_data_sources_tar_gz=funcotator_data_sources_tar_gz,
            blacklist_intervals=blacklist_intervals,
            preemptible_attempts = cnv_preemptible_attempts,
            num_changepoints_penalty_factor = num_changepoints_penalty_factor
    }

    call Manta_SomaticSV.MantaSomaticSV as MantaSomaticSV {
        input:
            sample_name=sample_name,
            tumor_bam=input_bam,
            tumor_bam_index=input_bam_index,
            ref_fasta=ref_fasta,
            ref_fasta_index=ref_fasta_index,
            is_cram=false,
            config_manta=config_manta,
            manta_docker=manta_docker
    }

    call manta_annot.run_manta_annotator as manta_annotator{
        input:
            sv = MantaSomaticSV.somatic_sv_vcf
    }

    call mutect2.Mutect2 as mutect2 {
        input:
            gatk_docker=gatk_docker,
            ref_dict=ref_dict,
            ref_fai=ref_fasta_index,
            ref_fasta=ref_fasta,
            scatter_count=M2scatter,
            intervals=mutect2_intervals,
            intervals_for_contamination=mutect2_intervals_for_contamination,
            tumor_name=sample_name,
            tumor_reads=input_bam,
            tumor_reads_index=input_bam_index,
            gcs_project_for_requester_pays=gcs_project_for_requester_pays,
            compress_vcfs=true,
            filter_funcotations=false,
            funco_compress=true,
            funco_filter_funcotations=false,
            funco_output_format="VCF",
            funco_reference_version="hg38",
            funco_data_sources_tar_gz=funcotator_data_sources_tar_gz,
            gnomad=gnomad,
            gnomad_idx=gnomad_idx,
            m2_extra_args=m2_extra_args,
            m2_extra_filtering_args=m2_filter_args,
            make_bamout=false,
            pon=pon,
            pon_idx=pon_idx,
            run_funcotator=true,
            run_orientation_bias_mixture_model_filter=true,
    }

    # call PureCN.PureCN as PureCN {
    #     input:
    #         sampleID=sample_name,
    #         segFile=CNVSomaticPairWorkflow.modeled_segments_tumor,
    #         vcf=mutect2.base_vcf,
    #         intervals=purecn_intervals,
    # }

    call msisensor2.msisensor2_workflow as msisensor2_workflow{
        input:
            sample_id=sample_name,
            bam=input_bam,
            bai=input_bam_index
    }

    call setGT.bcftools_fix_ploidy as set_GT {
        input:
            sample_id=sample_name,
            vcf=select_first([mutect2.funcotated_file, mutect2.base_vcf]),
    }

    call fixmutect2.fix_mutect2 as fix_mutect2 {
        input:
            sample_id=sample_name,
            vcf_file=set_GT.vcf_fixedploid
    }

    call removeFiltered.RemoveFiltered as RemoveFiltered {
        input:
            sample_id=sample_name,
            input_vcf=fix_mutect2.vcf_fixed,
            bcftools_exclude_string=bcftools_exclude_string
    }

    call openCravat.opencravat as open_cravat {
        input:
            vcf=RemoveFiltered.output_vcf,
            annotators_to_use=annotators,
            oc_modules=oc_modules,
    }

    call vcf_to_depmap.vcf_to_depmap as my_vcf_to_depmap {
        input:
            input_vcf=open_cravat.oc_main_file,
            sample_id=sample_name,
    }

    output {
        # #CNVSomaticPairWorkflow
        File read_counts_entity_id_tumor = CNVSomaticPairWorkflow.read_counts_entity_id_tumor
        File read_counts_tumor = CNVSomaticPairWorkflow.read_counts_tumor
        File allelic_counts_entity_id_tumor = CNVSomaticPairWorkflow.allelic_counts_entity_id_tumor
        File allelic_counts_tumor = CNVSomaticPairWorkflow.allelic_counts_tumor
        File denoised_copy_ratios_tumor = CNVSomaticPairWorkflow.denoised_copy_ratios_tumor
        File standardized_copy_ratios_tumor = CNVSomaticPairWorkflow.standardized_copy_ratios_tumor
        File het_allelic_counts_tumor = CNVSomaticPairWorkflow.het_allelic_counts_tumor
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
        File denoised_copy_ratios_plot_tumor = CNVSomaticPairWorkflow.denoised_copy_ratios_plot_tumor
        # File denoised_copy_ratios_lim_4_plot_tumor = CNVSomaticPairWorkflow.denoised_copy_ratios_lim_4_plot_tumor
        File standardized_MAD_tumor = CNVSomaticPairWorkflow.standardized_MAD_tumor
        Float standardized_MAD_value_tumor = CNVSomaticPairWorkflow.standardized_MAD_value_tumor
        File denoised_MAD_tumor = CNVSomaticPairWorkflow.denoised_MAD_tumor
        Float denoised_MAD_value_tumor = CNVSomaticPairWorkflow.denoised_MAD_value_tumor
        File delta_MAD_tumor = CNVSomaticPairWorkflow.delta_MAD_tumor
        Float delta_MAD_value_tumor = CNVSomaticPairWorkflow.delta_MAD_value_tumor
        File scaled_delta_MAD_tumor = CNVSomaticPairWorkflow.scaled_delta_MAD_tumor
        Float scaled_delta_MAD_value_tumor = CNVSomaticPairWorkflow.scaled_delta_MAD_value_tumor
        File modeled_segments_plot_tumor = CNVSomaticPairWorkflow.modeled_segments_plot_tumor
        #MantaSomaticSV
        File candidate_indel_vcf= MantaSomaticSV.candidate_indel_vcf
        File candidate_indel_vcf_index= MantaSomaticSV.candidate_indel_vcf_index
        File candidate_sv_vcf= MantaSomaticSV.candidate_sv_vcf
        File candidate_sv_vcf_index= MantaSomaticSV.candidate_sv_vcf_index
        File germline_sv_vcf= MantaSomaticSV.germline_sv_vcf
        File germline_sv_vcf_index= MantaSomaticSV.germline_sv_vcf_index
        File somatic_sv_vcf= MantaSomaticSV.somatic_sv_vcf
        File somatic_sv_vcf_index= MantaSomaticSV.somatic_sv_vcf_index
        #manta_annot
        File somatic_annotated_sv = manta_annotator.somatic_annotated_sv 
        File filtered_annotated_sv = manta_annotator.filtered_annotated_sv 
        File dropped_sv = manta_annotator.dropped 
        # omics_mutect2
        File omics_mutect2_out_vcf=fix_mutect2.vcf_fixed
        File mutect2_base_vcf = mutect2.base_vcf
        File full_vcf_idx=select_first([mutect2.funcotated_file_index, mutect2.base_vcf_idx])
        # # PureCN
        # File PureCN_solutions_pdf = PureCN.solutions_pdf
        # File chromosomes_pdf = PureCN.chromosomes_pdf
        # File PureCN_rds = PureCN.rds
        # File PureCN_dnacopy = PureCN.dnacopy
        # File PureCN_variants = PureCN.variants
        # File PureCN_loh = PureCN.loh
        # File PureCN_genes = PureCN.genes
        # File PureCN_segmentation = PureCN.segmentation
        # File PureCN_log = PureCN.log
        # File PureCN_selected_solution = PureCN.selected_solution
        # File PureCN_local_optima_pdf = PureCN.local_optima_pdf
        # String PureCN_purity = PureCN.purity
        # String PureCN_ploidy = PureCN.ploidy
        # String PureCN_contamination = PureCN.contamination
        # String PureCN_flagged = PureCN.flagged
        # String PureCN_curated = PureCN.curated
        # String PureCN_comment = PureCN.comment
        # String PureCN_wgd = PureCN.wgd
        # String PureCN_loh_fraction = PureCN.loh_fraction
        # String PureCN_cin = PureCN.cin
        # String PureCN_cin_allele_specific = PureCN.cin_allele_specific
        # String PureCN_cin_ploidy_robust = PureCN.cin_ploidy_robust
        # String PureCN_cin_allele_specific_ploidy_robust = PureCN.cin_allele_specific_ploidy_robust
        # msisensor2
        Float msisensor2_score=msisensor2_workflow.msisensor2_score
        File msisensor2_output=msisensor2_workflow.msisensor2_output
        File msisensor2_output_dis=msisensor2_workflow.msisensor2_output_dis
        File msisensor2_output_somatic=msisensor2_workflow.msisensor2_output_somatic
        # opencravat
        File oc_error_file=open_cravat.oc_error_file
        File oc_log_file=open_cravat.oc_log_file
        File oc_main_file=open_cravat.oc_main_file
        # File oc_sql_files=open_cravat.oc_sql_file
        # vcf_to_depmap
        Array[File] dna_pipeline_main_parquet=my_vcf_to_depmap.full_file
        File somatic_maf=my_vcf_to_depmap.depmap_maf
    }
}
