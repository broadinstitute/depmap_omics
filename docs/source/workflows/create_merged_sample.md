
## run_create_merged_sample

### Inputs

#### Required

  * `merge_mode` (String, **required**)
  * `sample_id` (String, **required**)
  * `vcfs` (Array[File], **required**)

#### Defaults

  * `create_merged_sample.boot_disk_size` (Int, default=10)
  * `create_merged_sample.disk_space` (Int, default=50)
  * `create_merged_sample.docker` (String, default="python")
  * `create_merged_sample.memory` (Int, default=4)
  * `create_merged_sample.num_preempt` (Int, default=5)
  * `create_merged_sample.num_threads` (Int, default=1)
  * `create_merged_sample.rna_sample_name` (String, default="none")
  * `create_merged_sample.size_tocheck` (Int, default=100)
  * `merge_vcfs.boot_disk_size` (Int, default=10)
  * `merge_vcfs.disk_space` (Int, default=40)
  * `merge_vcfs.docker` (String, default="biocontainers/bcftools:1.3")
  * `merge_vcfs.memory` (Int, default=4)
  * `merge_vcfs.num_preempt` (Int, default=5)
  * `merge_vcfs.num_threads` (Int, default=1)
  * `merge_vcfs.output_type` (String, default="z")

### Outputs

  * `vcf_updated` (File)

## create_merged_sample

author
: Jeremie Kalfon

### Inputs

#### Required

  * `sample_id` (String, **required**)
  * `vcf_file` (File, **required**)

#### Defaults

  * `boot_disk_size` (Int, default=10)
  * `disk_space` (Int, default=50)
  * `docker` (String, default="python")
  * `memory` (Int, default=4)
  * `num_preempt` (Int, default=5)
  * `num_threads` (Int, default=1)
  * `rna_sample_name` (String, default="none")
  * `size_tocheck` (Int, default=100)

### Outputs

  * `vcf_updated` (File)

## merge_vcfs

author
: Jeremie Kalfon

### Inputs

#### Required

  * `merge_mode` (String, **required**)
  * `sample_id` (String, **required**)
  * `vcfs` (Array[File], **required**)

#### Defaults

  * `boot_disk_size` (Int, default=10)
  * `disk_space` (Int, default=40)
  * `docker` (String, default="biocontainers/bcftools:1.3")
  * `memory` (Int, default=4)
  * `num_preempt` (Int, default=5)
  * `num_threads` (Int, default=1)
  * `output_type` (String, default="z")

### Outputs

  * `vcf_merged` (File)
  * `vcf_merged_index` (File)
