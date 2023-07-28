
## run_fix_mutect2

### Inputs

#### Required

  * `sample_id` (String, **required**)
  * `vcf` (File, **required**)

#### Defaults

  * `fix_mutect2.boot_disk_size` (Int, default=10)
  * `fix_mutect2.disk_space` (Int, default=50)
  * `fix_mutect2.docker` (String, default="python:3.10.7-bullseye")
  * `fix_mutect2.memory` (Int, default=4)
  * `fix_mutect2.num_preempt` (Int, default=5)
  * `fix_mutect2.num_threads` (Int, default=1)
  * `fix_mutect2.size_tocheck` (Int, default=100)

### Outputs

  * `vcf_fixed` (File)

## fix_mutect2

author
: Jeremie Kalfon

### Inputs

#### Required

  * `sample_id` (String, **required**)
  * `vcf_file` (File, **required**)

#### Defaults

  * `boot_disk_size` (Int, default=10)
  * `disk_space` (Int, default=50)
  * `docker` (String, default="python:3.10.7-bullseye")
  * `memory` (Int, default=4)
  * `num_preempt` (Int, default=5)
  * `num_threads` (Int, default=1)
  * `size_tocheck` (Int, default=100)

### Outputs

  * `vcf_fixed` (File)
