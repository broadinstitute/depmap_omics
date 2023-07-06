
## update_solution

### Inputs

#### Required

  * `purecn_rds` (File, **required**)
  * `sampleID` (String, **required**)
  * `solution_num` (Int, **required**)

#### Optional

  * `run_update_solution.hardware_disk_size_GB` (Int?)
  * `run_update_solution.hardware_memory_GB` (Int?)
  * `run_update_solution.hardware_preemptible_tries` (Int?)

#### Defaults

  * `run_update_solution.purecn_docker` (String, default="markusriester/purecn:2.2.0")

### Outputs

  * `variants` (File)
  * `loh` (File)
  * `genes` (File)
  * `dnacopy` (File)
  * `selected_solution` (File)
  * `purity` (String)
  * `ploidy` (String)
  * `contamination` (String)
  * `flagged` (String)
  * `comment` (String)
  * `wgd` (String)
  * `loh_fraction` (String)
  * `cin` (String)
  * `cin_allele_specific` (String)
  * `cin_ploidy_robust` (String)
  * `cin_allele_specific_ploidy_robust` (String)

## run_update_solution

### Inputs

#### Required

  * `purecn_rds` (File, **required**)
  * `sampleID` (String, **required**)
  * `solution_num` (Int, **required**)

#### Optional

  * `hardware_disk_size_GB` (Int?)
  * `hardware_memory_GB` (Int?)
  * `hardware_preemptible_tries` (Int?)

#### Defaults

  * `purecn_docker` (String, default="markusriester/purecn:2.2.0")

### Outputs

  * `dnacopy` (File)
  * `variants` (File)
  * `loh` (File)
  * `genes` (File)
  * `selected_solution` (File)
  * `table` (Array[String])
  * `purity` (String)
  * `ploidy` (String)
  * `contamination` (String)
  * `flagged` (String)
  * `comment` (String)
  * `wgd_table` (Array[String])
  * `wgd` (String)
  * `loh_fraction` (String)
  * `cin` (String)
  * `cin_allele_specific` (String)
  * `cin_ploidy_robust` (String)
  * `cin_allele_specific_ploidy_robust` (String)
