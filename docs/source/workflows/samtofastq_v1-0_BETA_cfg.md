
## samtofastq_workflow

### Inputs

#### Required

  * `samtofastq.disk_space` (Int, **required**)
  * `samtofastq.input_bam_cram` (File, **required**)
  * `samtofastq.memory` (Float, **required**)
  * `samtofastq.num_preempt` (Int, **required**)
  * `samtofastq.num_threads` (Int, **required**)
  * `samtofastq.prefix` (String, **required**)

#### Optional

  * `samtofastq.reference_fasta` (File?)

#### Defaults

  * `samtofastq.java_memory` (Int, default=floor((memory - 0.5)))

### Outputs

  * `samtofastq.fastq1` (File)
  * `samtofastq.fastq2` (File)

## samtofastq

author
: Francois Aguet

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
