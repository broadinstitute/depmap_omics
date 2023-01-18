
## run_opencravat

### Inputs

#### Required

  * `vcf` (File, **required**)

#### Optional

  * `opencravat.oncokb_api_key` (File?)

#### Defaults

  * `opencravat.annotators_to_use` (Array[String], default=[])
  * `opencravat.boot_disk_size` (Int, default=100)
  * `opencravat.disk_space` (Int, default=20)
  * `opencravat.docker` (String, default="karchinlab/opencravat")
  * `opencravat.format` (String, default="vcf")
  * `opencravat.genome` (String, default="hg38")
  * `opencravat.memory` (Int, default=16)
  * `opencravat.modules_options` (String, default="vcfreporter.type=separate")
  * `opencravat.num_preempt` (Int, default=2)
  * `opencravat.num_threads` (Int, default=4)
  * `opencravat.retries` (Int, default=1)

### Outputs

  * `oc_error_file` (File)
  * `oc_log_file` (File)
  * `oc_main_file` (File)

## opencravat

author
: Jeremie Kalfon

### Inputs

#### Required

  * `vcf` (File, **required**)

#### Optional

  * `oncokb_api_key` (File?)

#### Defaults

  * `annotators_to_use` (Array[String], default=[])
  * `boot_disk_size` (Int, default=100)
  * `disk_space` (Int, default=20)
  * `docker` (String, default="karchinlab/opencravat")
  * `format` (String, default="vcf")
  * `genome` (String, default="hg38")
  * `memory` (Int, default=16)
  * `modules_options` (String, default="vcfreporter.type=separate")
  * `num_preempt` (Int, default=2)
  * `num_threads` (Int, default=4)
  * `retries` (Int, default=1)

### Outputs

  * `oc_error_file` (File)
  * `oc_log_file` (File)
  * `oc_main_file` (File)
