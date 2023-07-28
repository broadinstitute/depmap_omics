
## aggregate_vcfs

### Inputs

#### Required

  * `aggregate.disk_space` (Int, **required**)
  * `aggregate.memory` (Int, **required**)
  * `aggregate.num_preempt` (Int, **required**)
  * `aggregate.num_threads` (Int, **required**)
  * `aggregate.option` (String, **required**)
  * `aggregate.sample_set_id` (String, **required**)
  * `aggregate.vcf_files` (Array[File], **required**)

#### Optional

  * `aggregate.vcf_indexes` (Array[File]?)

### Outputs

  * `aggregate.merged_vcf` (File)

## aggregate

author
: Jeremie Kalfon

### Inputs

#### Required

  * `disk_space` (Int, **required**)
  * `memory` (Int, **required**)
  * `num_preempt` (Int, **required**)
  * `num_threads` (Int, **required**)
  * `option` (String, **required**)
  * `sample_set_id` (String, **required**)
  * `vcf_files` (Array[File], **required**)

#### Optional

  * `vcf_indexes` (Array[File]?)

### Outputs

  * `merged_vcf` (File)
