import "samtofastq_v1-0_BETA_cfg.wdl" as samtofastq_v1
import "star_v1-0_BETA_cfg.wdl" as star_v1
import "rnaseqc2_v1-0_BETA_cfg.wdl" as rnaseqc2_v1
import "rsem_depmap.wdl" as rsem_v1
import "hg38_STAR_fusion.wdl" as hg38_STAR_fusion


workflow RNA_pipeline {

  #samtofastq_v1
  File input_bam_cram
  File reference_fasta

  #rnaseqc2_v1
  File genes_gtf
  String sample_id

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
      fastq2=samtofastq.fastq2
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
      prefix=sample_id
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

