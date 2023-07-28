
## msisensor2_workflow

### Inputs

#### Required

  * `bai` (File, **required**)
  * `bam` (File, **required**)
  * `sample_id` (String, **required**)

#### Defaults

  * `disk_space` (Int, default=250)
  * `memory` (Int, default=6)
  * `num_preempt` (Int, default=1)
  * `num_threads` (Int, default=2)

### Outputs

  * `msisensor2_score` (Float)
  * `msisensor2_output` (File)
  * `msisensor2_output_dis` (File)
  * `msisensor2_output_somatic` (File)

## msisensor2

author
: David Wu

### Inputs

#### Required

  * `bai` (File, **required**)
  * `bam` (File, **required**)
  * `sample_id` (String, **required**)

#### Defaults

  * `disk_space` (Int, default=250)
  * `memory` (Int, default=6)
  * `num_preempt` (Int, default=1)
  * `num_threads` (Int, default=2)

### Outputs

  * `msisensor2_score` (Float)
  * `msisensor2_output` (File)
  * `msisensor2_output_dis` (File)
  * `msisensor2_output_somatic` (File)
