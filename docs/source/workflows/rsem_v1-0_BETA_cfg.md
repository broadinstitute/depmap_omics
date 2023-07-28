
## rsem_workflow

### Inputs

#### Required

  * `rsem.disk_space` (Int, **required**)
  * `rsem.memory` (Int, **required**)
  * `rsem.num_preempt` (Int, **required**)
  * `rsem.num_threads` (Int, **required**)
  * `rsem.prefix` (String, **required**)
  * `rsem.rsem_reference` (File, **required**)
  * `rsem.transcriptome_bam` (File, **required**)

#### Optional

  * `rsem.estimate_rspd` (String?)
  * `rsem.is_stranded` (String?)
  * `rsem.max_frag_len` (Int?)
  * `rsem.paired_end` (String?)

### Outputs

  * `rsem.genes` (File)
  * `rsem.isoforms` (File)

## rsem

author
: Francois Aguet

### Inputs

#### Required

  * `disk_space` (Int, **required**)
  * `memory` (Int, **required**)
  * `num_preempt` (Int, **required**)
  * `num_threads` (Int, **required**)
  * `prefix` (String, **required**)
  * `rsem_reference` (File, **required**)
  * `transcriptome_bam` (File, **required**)

#### Optional

  * `estimate_rspd` (String?)
  * `is_stranded` (String?)
  * `max_frag_len` (Int?)
  * `paired_end` (String?)

### Outputs

  * `genes` (File)
  * `isoforms` (File)
