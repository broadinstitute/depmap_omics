version 1.0

import "samtofastq_wdl1-0.wdl" as samtofastq_v1
import "star_wdl1-0.wdl" as star_v1
import "rnaseqc2_wdl1-0.wdl" as rnaseqc2_v1
import "rsem_depmap.wdl" as rsem_v1
import "hg38_STAR_fusion_wdl1-0.wdl" as hg38_STAR_fusion
#import "rnaseq_mutect2_tumor_only.wdl" as rna_mutect2


workflow RNA_pipeline {

  input {
  #samtofastq_v1
  File input_bam_cram
  File reference_fasta

  #star_v1
  File star_index

  #star fusion
  Array[File] ctat_genome_lib_build_dir_files
  Array[File] ref_genome_fa_star_idx_files

  #rsem
  File rsem_reference

  #rnaseqc2_v1
  File genes_gtf
  String sample_id

  #rna_mutect2
#   Boolean run_funcotator=false
#   String gatk_docker="broadinstitute/gatk:4.2.2.0"
#   Array[File] knownVcfs
#   Array[File] knownVcfsIndices
#   File dbSnpVcf
#   File dbSnpVcfIndex
#   File ref_fasta
#   File ref_fai
#   File ref_dict
#   String source
#   Boolean compress_vcfs=true
#   Boolean filter_funcotations=false
#   Boolean funco_compress=true
#   String funco_output_format="VCF"
#   String funco_reference_version="hg38"
#   Boolean funco_use_gnomad_AF=true
#   Int scatter_count=20
#   String gcs_project_for_requester_pays="broad-firecloud-ccle"
#   File intervals="gs://gcp-public-data--broad-references/hg38/v0/wgs_coverage_regions.hg38.interval_list"
#   File gnomad = "gs://gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz"
#   File gnomad_idx = "gs://gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi"
#   String m2_extra_args="--genotype-germline-sites true --genotype-pon-sites true"
#   Boolean make_bamout=false
#   File pon="gs://gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz"
#   File pon_idx="gs://gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz.tbi"
  }

  call samtofastq_v1.samtofastq as samtofastq {
    input:
      input_bam_cram=input_bam_cram,
      prefix=sample_id,
      reference_fasta=reference_fasta
  }

  call star_v1.star as star {
    input:
      prefix=sample_id,
      fastq1=samtofastq.fastq1,
      fastq2=samtofastq.fastq2,
      star_index=star_index
  }

  call rnaseqc2_v1.rnaseqc2 as rnaseqc2 {
    input:
      bam_file=star.bam_file,
      genes_gtf=genes_gtf,
      sample_id=sample_id
  }

  call rsem_v1.rsem as rsem {
    input:
      transcriptome_bam=star.transcriptome_bam,
      prefix=sample_id,
      rsem_reference=rsem_reference,
      is_stranded="false",
      paired_end="true"
  }

  call hg38_STAR_fusion.StarFusion as StarFusion {
    input:
      left_fastq=samtofastq.fastq1,
      right_fastq=samtofastq.fastq2,
      prefix=sample_id,
      ctat_genome_lib_build_dir_files=ctat_genome_lib_build_dir_files,
      ref_genome_fa_star_idx_files=ref_genome_fa_star_idx_files
  }

#   call rna_mutect2.RNAseq_mutect2 as RNAseq_mutect2{
#       input:
#         run_funcotator=run_funcotator,
#         gatk_docker=gatk_docker,
#         knownVcfs=knownVcfs,
#         knownVcfsIndices=knownVcfsIndices,
#         dbSnpVcf=dbSnpVcf,
#         dbSnpVcfIndex=dbSnpVcfIndex,
#         intervals=intervals,
#         ref_fasta=ref_fasta,
#         ref_fai=ref_fasta,
#         ref_dict=ref_dict,
#         tumor_bam=star.bam_file,
#         tumor_bai=star.bam_index,
#         tumor_name=sample_id,
#         sequencing_center=source,
#         sequence_source=source,
#         compress_vcfs=compress_vcfs,
#         filter_funcotations=filter_funcotations,
#         funco_compress=funco_compress,
#         funco_output_format=funco_output_format,
#         funco_reference_version=funco_reference_version,
#         funco_use_gnomad_AF=funco_use_gnomad_AF,
#         scatter_count=scatter_count,
#         gcs_project_for_requester_pays=gcs_project_for_requester_pays,
#         gnomad = gnomad,
#         gnomad_idx = gnomad_idx,
#         m2_extra_args=m2_extra_args,
#         make_bamout=make_bamout,
#         pon=pon,
#         pon_idx=pon_idx,
#   }

  output {
    #samtofastq
    #star
    File bam_file=star.bam_file
    File bam_index=star.bam_index
    File transcriptome_bam=star.transcriptome_bam
    File chimeric_junctions=star.chimeric_junctions
    File chimeric_bam_file=star.chimeric_bam_file
    File read_counts=star.read_counts
    File junctions=star.junctions
    File junctions_pass1=star.junctions_pass1
    Array[File] logs=star.logs
    #rnaseqc
    File gene_tpm=rnaseqc2.gene_tpm
    File gene_counts=rnaseqc2.gene_counts
    File exon_counts=rnaseqc2.exon_counts
    File metrics=rnaseqc2.metrics
    File insertsize_distr=rnaseqc2.insertsize_distr
    #rsem
    File genes=rsem.genes
    File isoforms=rsem.isoforms
    #StarFusion
    File fusion_predictions=StarFusion.fusion_predictions
    File fusion_predictions_abridged=StarFusion.fusion_predictions_abridged
    #rna_mutect2
    # File merged_vcf = RNAseq_mutect2.merged_vcf
    # File merged_vcf_index = RNAseq_mutect2.merged_vcf_index
    # File variant_filtered_vcf = RNAseq_mutect2.variant_filtered_vcf
    # File variant_filtered_vcf_index = RNAseq_mutect2.variant_filtered_vcf_index
    # File filtering_stats = RNAseq_mutect2.filtering_stats
    # File mutect_stats = RNAseq_mutect2.mutect_stats

    # File? funcotated_file = RNAseq_mutect2.funcotated_file
    # File? funcotated_file_index = RNAseq_mutect2.funcotated_file_index
  }
}
