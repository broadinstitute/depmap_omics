
## rnaseqc2_workflow

### Inputs

#### Required

  * `bam_file` (File, **required**)
  * `genes_gtf` (File, **required**)
  * `sample_id` (String, **required**)

#### Optional

  * `flags` (String?)
  * `intervals_bed` (File?)
  * `strandedness` (String?)

#### Defaults

  * `rnaseqc2.disk_space` (Int, default=75)
  * `rnaseqc2.memory` (Int, default=8)
  * `rnaseqc2.num_preempt` (Int, default=1)
  * `rnaseqc2.num_threads` (Int, default=4)

### Outputs

  * `gene_tpm` (File)
  * `gene_counts` (File)
  * `exon_counts` (File)
  * `metrics` (File)
  * `insertsize_distr` (File)

## rnaseqc2

author
: Francois Aguet

### Inputs

#### Required

  * `bam_file` (File, **required**)
  * `genes_gtf` (File, **required**)
  * `sample_id` (String, **required**)

#### Optional

  * `flags` (String?)
  * `intervals_bed` (File?)
  * `strandedness` (String?)

#### Defaults

  * `disk_space` (Int, default=75)
  * `memory` (Int, default=8)
  * `num_preempt` (Int, default=1)
  * `num_threads` (Int, default=4)

### Outputs

  * `gene_tpm` (File)
  * `gene_counts` (File)
  * `exon_counts` (File)
  * `metrics` (File)
  * `insertsize_distr` (File)
