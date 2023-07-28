
## rnaseqc2_workflow

### Inputs

#### Required

  * `rnaseqc2.bam_file` (File, **required**)
  * `rnaseqc2.disk_space` (Int, **required**)
  * `rnaseqc2.genes_gtf` (File, **required**)
  * `rnaseqc2.memory` (Int, **required**)
  * `rnaseqc2.num_preempt` (Int, **required**)
  * `rnaseqc2.num_threads` (Int, **required**)
  * `rnaseqc2.sample_id` (String, **required**)

#### Optional

  * `rnaseqc2.flags` (String?)
  * `rnaseqc2.intervals_bed` (File?)
  * `rnaseqc2.strandedness` (String?)

### Outputs

  * `rnaseqc2.gene_tpm` (File)
  * `rnaseqc2.gene_counts` (File)
  * `rnaseqc2.exon_counts` (File)
  * `rnaseqc2.metrics` (File)
  * `rnaseqc2.insertsize_distr` (File)

## rnaseqc2

author
: Francois Aguet

### Inputs

#### Required

  * `bam_file` (File, **required**)
  * `disk_space` (Int, **required**)
  * `genes_gtf` (File, **required**)
  * `memory` (Int, **required**)
  * `num_preempt` (Int, **required**)
  * `num_threads` (Int, **required**)
  * `sample_id` (String, **required**)

#### Optional

  * `flags` (String?)
  * `intervals_bed` (File?)
  * `strandedness` (String?)

### Outputs

  * `gene_tpm` (File)
  * `gene_counts` (File)
  * `exon_counts` (File)
  * `metrics` (File)
  * `insertsize_distr` (File)
