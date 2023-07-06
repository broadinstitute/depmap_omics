
## run_RemoveFiltered

### Inputs

#### Required

  * `input_vcf` (File, **required**)
  * `sample_id` (String, **required**)

#### Defaults

  * `RemoveFiltered.bcftools_exclude_string` (String, default='FILTER~"weak_evidence"||FILTER~"map_qual"||FILTER~"strand_bias"||FILTER~"slippage"||FILTER~"clustered_events"||FILTER~"base_qual"')
  * `RemoveFiltered.boot_disk_size` (Int, default=10)
  * `RemoveFiltered.cpu` (Int, default=2)
  * `RemoveFiltered.docker_image` (String, default="dceoy/bcftools")
  * `RemoveFiltered.mem` (Int, default=2)
  * `RemoveFiltered.preemptible` (Int, default=3)

### Outputs

  * `output_vcf` (File)
  * `output_vcf_idx` (File)

## RemoveFiltered

### Inputs

#### Required

  * `input_vcf` (File, **required**)
  * `sample_id` (String, **required**)

#### Defaults

  * `bcftools_exclude_string` (String, default='FILTER~"weak_evidence"||FILTER~"map_qual"||FILTER~"strand_bias"||FILTER~"slippage"||FILTER~"clustered_events"||FILTER~"base_qual"')
  * `boot_disk_size` (Int, default=10)
  * `cpu` (Int, default=2)
  * `docker_image` (String, default="dceoy/bcftools")
  * `mem` (Int, default=2)
  * `preemptible` (Int, default=3)

### Outputs

  * `output_vcf` (File)
  * `output_vcf_idx` (File)
