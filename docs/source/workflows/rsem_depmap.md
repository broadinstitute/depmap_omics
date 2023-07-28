
## rsem_workflow

### Inputs

#### Required

  * `prefix` (String, **required**)
  * `rsem_reference` (File, **required**)
  * `transcriptome_bam` (File, **required**)

#### Optional

  * `rsem.calc_ci` (String?)
  * `rsem.ci_memory` (Int?)
  * `rsem.is_stranded` (String?)
  * `rsem.paired_end` (String?)

#### Defaults

  * `rsem.disk_space` (Int, default=500)
  * `rsem.estimate_rspd` (String, default="true")
  * `rsem.max_frag_len` (Int, default=1000)
  * `rsem.memory` (Int, default=128)
  * `rsem.num_preempt` (Int, default=1)
  * `rsem.num_threads` (Int, default=32)

### Outputs

  * `genes` (File)
  * `isoforms` (File)

## rsem

author
: David Wu

### Inputs

#### Required

  * `prefix` (String, **required**)
  * `rsem_reference` (File, **required**)
  * `transcriptome_bam` (File, **required**)

#### Optional

  * `calc_ci` (String?)
  * `ci_memory` (Int?)
  * `is_stranded` (String?)
  * `paired_end` (String?)

#### Defaults

  * `disk_space` (Int, default=500)
  * `estimate_rspd` (String, default="true")
  * `max_frag_len` (Int, default=1000)
  * `memory` (Int, default=128)
  * `num_preempt` (Int, default=1)
  * `num_threads` (Int, default=32)

### Outputs

  * `genes` (File)
  * `isoforms` (File)
