version 1.0

import "rnaseqc2_wdl1-0.wdl" as rnaseqc2_v1
import "rsem_depmap.wdl" as rsem_v1


workflow RNA_pipeline {

  input {
  #samtofastq_v1
  File bam_file
  File star_transcriptome_bam

  String is_stranded = "true"
  String is_stranded_rnaseqc = "RF"

  #rsem
  File rsem_reference

  #rnaseqc2_v1
  File genes_gtf
  String sample_id
  }

 call rnaseqc2_v1.rnaseqc2 as rnaseqc2 {
   input:
     bam_file=bam_file,
     genes_gtf=genes_gtf,
     sample_id=sample_id,
     strandedness=is_stranded_rnaseqc
 }

  call rsem_v1.rsem as rsem {
    input:
      transcriptome_bam=star_transcriptome_bam,
      prefix=sample_id,
      rsem_reference=rsem_reference,
      is_stranded=is_stranded,
      paired_end="true"
  }

  output {
    #rnaseqc
    File gene_tpm=rnaseqc2.gene_tpm
    File gene_counts=rnaseqc2.gene_counts
    File exon_counts=rnaseqc2.exon_counts
    File metrics=rnaseqc2.metrics
    File insertsize_distr=rnaseqc2.insertsize_distr
    #rsem
    File genes=rsem.genes
    File isoforms=rsem.isoforms
  }
}
