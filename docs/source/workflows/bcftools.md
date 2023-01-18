
## run_fix_ploidy

### Inputs

#### Required

  * `sample_id` (String, **required**)
  * `vcf` (File, **required**)

#### Defaults

  * `bcftools_fix_ploidy.boot_disk_size` (Int, default=10)
  * `bcftools_fix_ploidy.disk_space` (Int, default=40)
  * `bcftools_fix_ploidy.docker` (String, default="dceoy/bcftools")
  * `bcftools_fix_ploidy.memory` (Int, default=4)
  * `bcftools_fix_ploidy.num_preempt` (Int, default=5)
  * `bcftools_fix_ploidy.num_threads` (Int, default=1)

### Outputs

  * `bcftools_fix_ploidy.vcf_fixedploid` (File)

## bcftools_fix_ploidy

author
: Jeremie Kalfon

### Inputs

#### Required

  * `sample_id` (String, **required**)
  * `vcf` (File, **required**)

#### Defaults

  * `boot_disk_size` (Int, default=10)
  * `disk_space` (Int, default=40)
  * `docker` (String, default="dceoy/bcftools")
  * `memory` (Int, default=4)
  * `num_preempt` (Int, default=5)
  * `num_threads` (Int, default=1)

### Outputs

  * `vcf_fixedploid` (File)
