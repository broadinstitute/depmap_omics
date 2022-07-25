version 1.0

import "samtofastq_wdl1-0.wdl" as samtofastq_v1
import "star_wdl1-0.wdl" as star_v1
import "rnaseqc2_wdl1-0.wdl" as rnaseqc2_v1
import "rsem_depmap.wdl" as rsem_v1
import "hg38_STAR_fusion.wdl" as hg38_STAR_fusion


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

  #rnaseqc2_v1
  File genes_gtf
  String sample_id
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
      is_stranded="false",
      paired_end="true"
  }

  call hg38_STAR_fusion.StarFusion as StarFusion {
    input:
      left_fastq=samtofastq.fastq1,
      right_fastq=samtofastq.fastq2,
      prefix=sample_id,
      ctat_genome_lib_build_dir_files=ctat_genome_lib_build_dir_files],
      ref_genome_fa_star_idx_files=ref_genome_fa_star_idx_files
  }

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
  }
}

