version 1.0

import "https://raw.githubusercontent.com/NCIP/ctat-mutations/Terra-3.3.0/WDL/ctat_mutations.wdl" as CTAT_Mutations_wf




struct Ctat_mutations_config {


  File gtf
  File ref_bed
  File ref_fasta
  File ref_fasta_index

  File ref_dict
  
  String genome_version
    
  File cravat_lib_tar_gz

  File db_snp_vcf
  File db_snp_vcf_index
  
  File cosmic_vcf
  File cosmic_vcf_index
  
  File gnomad_vcf
  File gnomad_vcf_index

  File ref_splice_adj_regions_bed

  File repeat_mask_bed

  File rna_editing_vcf
  File rna_editing_vcf_index
  
  File star_reference
}

  

workflow ctat_mutations_Terra {


  input {
    String docker
    String sample_id
    File left
    File? right
    File? intervals
    Boolean annotate_variants = true
    String boosting_method = "none"
    Ctat_mutations_config pipe_inputs_config
  }
  
  call CTAT_Mutations_wf.ctat_mutations as CM_wf {

    input:
      docker = docker,
      sample_id = sample_id,
      left = left,
      right = right,

      intervals = intervals,
      annotate_variants = annotate_variants,
      boosting_method = boosting_method,
      
      gtf = pipe_inputs_config.gtf,
      ref_bed = pipe_inputs_config.ref_bed,
      ref_fasta = pipe_inputs_config.ref_fasta,
      ref_fasta_index = pipe_inputs_config.ref_fasta_index,
      ref_dict = pipe_inputs_config.ref_dict,
      genome_version = pipe_inputs_config.genome_version,
      cravat_lib_tar_gz = pipe_inputs_config.cravat_lib_tar_gz,
      db_snp_vcf = pipe_inputs_config.db_snp_vcf,
      db_snp_vcf_index = pipe_inputs_config.db_snp_vcf_index,
      cosmic_vcf = pipe_inputs_config.cosmic_vcf,
      cosmic_vcf_index = pipe_inputs_config.cosmic_vcf_index,
      gnomad_vcf = pipe_inputs_config.gnomad_vcf,
      gnomad_vcf_index = pipe_inputs_config.gnomad_vcf_index,
      ref_splice_adj_regions_bed = pipe_inputs_config.ref_splice_adj_regions_bed,
      repeat_mask_bed = pipe_inputs_config.repeat_mask_bed,
      rna_editing_vcf = pipe_inputs_config.rna_editing_vcf,
      rna_editing_vcf_index = pipe_inputs_config.rna_editing_vcf_index,
      star_reference = pipe_inputs_config.star_reference

  }


    output {
        File? haplotype_caller_vcf = CM_wf.haplotype_caller_vcf
        File? annotated_vcf = CM_wf.annotated_vcf
        File? filtered_vcf = CM_wf.filtered_vcf
        File? aligned_bam = CM_wf.aligned_bam
        File? output_log_final =  CM_wf.output_log_final
        File? output_SJ =  CM_wf.output_SJ
        File? recalibrated_bam = CM_wf.recalibrated_bam
        File? recalibrated_bam_index = CM_wf.recalibrated_bam_index
        File? cancer_igv_report = CM_wf.cancer_igv_report
        File? cancer_variants_tsv = CM_wf.cancer_variants_tsv
        File? cancer_vcf = CM_wf.cancer_vcf
    }
}


