
## samtofastq_workflow

### Inputs

#### Required

  * `disk_space` (Int, **required**)
  * `input_bam_cram` (File, **required**)
  * `memory` (Float, **required**)
  * `num_preempt` (Int, **required**)
  * `num_threads` (Int, **required**)
  * `prefix` (String, **required**)

#### Optional

  * `reference_fasta` (File?)

#### Defaults

  * `java_memory` (Int, default=floor((memory - 0.5)))

### Outputs

  * `fastq1` (File)
  * `fastq2` (File)

## samtofastq

author
: Francois Aguet

### Inputs

#### Required

  * `input_bam_cram` (File, **required**)
  * `prefix` (String, **required**)

#### Optional

  * `reference_fasta` (File?)

#### Defaults

  * `disk_space` (Int, default=400)
  * `java_memory` (Int, default=floor((memory - 0.5)))
  * `memory` (Float, default=64)
  * `num_preempt` (Int, default=5)
  * `num_threads` (Int, default=8)

### Outputs

  * `fastq1` (File)
  * `fastq2` (File)
