version 1.0


import "https://raw.githubusercontent.com/NCIP/ctat-mutations/Terra-3.3.0/WDL/Terra/ctat_mutations.Terra.wdl" as CTAT_Mutations_Terra


workflow ctat_mutations_Terra_hg38 {

  input {

    String docker = "trinityctat/ctat_mutations"
    String sample_id
    File left
    File? right
    File? intervals
    Boolean annotate_variants = true
    String boosting_method = "none"
  

    String gs_base_url = "gs://ctat_genome_libs/__genome_libs_StarFv1.10/GRCh38_gencode_v22_CTAT_lib_Mar012021.plug-n-play"


    Ctat_mutations_config pipe_inputs_config = {
      "genome_version" : "hg38",
      "gtf" : gs_base_url + "/ref_annot.gtf",
      "ref_bed" : gs_base_url + "/ctat_mutation_lib/refGene.sort.bed.gz",
      "ref_fasta" : gs_base_url + "/ref_genome.fa",
      "ref_fasta_index" : gs_base_url + "/ref_genome.fa.fai",
      "ref_dict" : gs_base_url + "/ref_genome.dict",
      "cravat_lib_tar_gz" : gs_base_url + "/ctat_mutation_lib/cravat.tar.bz2",
      "db_snp_vcf" : gs_base_url + "/ctat_mutation_lib/dbsnp.vcf.gz",
      "db_snp_vcf_index" : gs_base_url + "/ctat_mutation_lib/dbsnp.vcf.gz.tbi",
      "cosmic_vcf" : gs_base_url + "/ctat_mutation_lib/cosmic.vcf.gz",
      "cosmic_vcf_index" : gs_base_url + "/ctat_mutation_lib/cosmic.vcf.gz.csi",
      "gnomad_vcf" : gs_base_url + "/ctat_mutation_lib/gnomad-lite.vcf.gz",
      "gnomad_vcf_index" : gs_base_url + "/ctat_mutation_lib/gnomad-lite.vcf.gz.csi",
      "ref_splice_adj_regions_bed" : gs_base_url + "/ctat_mutation_lib/ref_annot.splice_adj.bed.gz",
      "repeat_mask_bed" : gs_base_url + "/ctat_mutation_lib/repeats_ucsc_gb.bed.gz",
      "rna_editing_vcf" : gs_base_url + "/ctat_mutation_lib/RNAediting.library.vcf.gz",
      "rna_editing_vcf_index" : gs_base_url + "/ctat_mutation_lib/RNAediting.library.vcf.gz.csi",
      "star_reference" : gs_base_url + "/ref_genome.fa.star.idx.tar.bz2"

    }

  }

  call CTAT_Mutations_Terra.ctat_mutations_Terra as CM_Terra_wf {

    input:
      docker = docker,
      sample_id = sample_id,
      left = left,
      right = right,
      intervals = intervals,
      annotate_variants = annotate_variants,
      boosting_method = boosting_method,
      pipe_inputs_config = pipe_inputs_config

  }

 output {
        File? haplotype_caller_vcf = CM_Terra_wf.haplotype_caller_vcf
        File? annotated_vcf = CM_Terra_wf.annotated_vcf
        File? filtered_vcf = CM_Terra_wf.filtered_vcf
        File? aligned_bam = CM_Terra_wf.aligned_bam
        File? output_log_final =  CM_Terra_wf.output_log_final
        File? output_SJ =  CM_Terra_wf.output_SJ
        File? recalibrated_bam = CM_Terra_wf.recalibrated_bam
        File? recalibrated_bam_index = CM_Terra_wf.recalibrated_bam_index
        File? cancer_igv_report = CM_Terra_wf.cancer_igv_report
        File? cancer_variants_tsv = CM_Terra_wf.cancer_variants_tsv
        File? cancer_vcf = CM_Terra_wf.cancer_vcf

 }
}
