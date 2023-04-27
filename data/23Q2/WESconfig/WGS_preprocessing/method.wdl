version 1.0
# converts hg19 bams to hg38, and preprocesses it for the following steps (WGS_pipeline)
import "bam-to-unmapped-bams.wdl" as https://raw.githubusercontent.com/gatk-workflows/seq-format-conversion/3.0.0/bam-to-unmapped-bams.wdl
import "processing-for-variant-discovery-gatk4.wdl" as https://raw.githubusercontent.com/gatk-workflows/gatk4-data-processing/2.1.1/processing-for-variant-discovery-gatk4.wdl
import "ArrayOfFilesToTxt.wdl" as https://raw.githubusercontent.com/broadinstitute/depmap_omics/dev/WGS_pipeline/ArrayOfFilesToTxt.wdl

workflow WGS_preprocessing {

    input {
        #BamToUnmappedRGBams

        File input_bam
        File input_bam_index

        String sample_name #this.name

        #PreProcessingForVariantDiscovery_GATK4
        #https://github.com/gatk-workflows/gatk4-data-processing/blob/c44603c464fe3cb7d9b82da2a95f844fdeb20e3c/processing-for-variant-discovery-gatk4.wdl
        String ref_name
        File ref_fasta
        File ref_fasta_index
        File ref_dict
        File? ref_alt
        File ref_sa
        File ref_ann
        File ref_bwt
        File ref_pac
        File ref_amb

        Int preemptible_tries
        
        File dbSNP_vcf
        File dbSNP_vcf_index

        Array[File] known_indels_sites_VCFs
        Array[File] known_indels_sites_indices

        String unmapped_bam_suffix = ".bam"
    }

    call BamToUnmappedRGBams.BamToUnmappedBams as BamToUnmappedBams {
        input:
            input_bam=input_bam,
    }

    call ArrayOfFilesToTxt.CreateTxt as CreateTxt {
        input:
            array_of_files=BamToUnmappedBams.output_bams,
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
            dbSNP_vcf=dbSNP_vcf,
            dbSNP_vcf_index=dbSNP_vcf_index,
            known_indels_sites_VCFs=known_indels_sites_VCFs,
            known_indels_sites_indices=known_indels_sites_indices,
            preemptible_tries=preemptible_tries,
            ref_sa = ref_sa,
            ref_ann = ref_ann,
            ref_bwt = ref_bwt,
            ref_pac = ref_pac,
            ref_amb = ref_amb,
            ref_alt = ref_alt
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

    }







}
