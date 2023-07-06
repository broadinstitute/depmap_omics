
## run_PureCN

### Inputs

#### Required

  * `genome` (String, **required**)
  * `intervals` (File, **required**)
  * `sampleID` (String, **required**)
  * `segFile` (File, **required**)
  * `vcf` (File, **required**)

#### Defaults

  * `PureCN.funSegmentation` (String, default="Hclust")
  * `PureCN.hardware_disk_size_GB` (Int, default=250)
  * `PureCN.hardware_memory_GB` (Int, default=32)
  * `PureCN.hardware_preemptible_tries` (Int, default=2)
  * `PureCN.maxCopyNumber` (Int, default=8)
  * `PureCN.maxPurity` (Float, default=0.99)
  * `PureCN.maxSegments` (Int, default=1000)
  * `PureCN.max_retries` (Int, default=0)
  * `PureCN.minPurity` (Float, default=0.9)
  * `PureCN.num_threads` (Int, default=1)
  * `PureCN.otherArguments` (String, default="--post-optimize --model-homozygous --min-total-counts 20")
  * `PureCN.purecn_docker` (String, default="markusriester/purecn:2.2.0")

### Outputs

  * `solutions_pdf` (File)
  * `chromosomes_pdf` (File)
  * `rds` (File)
  * `dnacopy` (File)
  * `variants` (File)
  * `loh` (File)
  * `genes` (File)
  * `segmentation` (File)
  * `log` (File)
  * `selected_solution` (File)
  * `local_optima_pdf` (File)
  * `purity` (String)
  * `ploidy` (String)
  * `contamination` (String)
  * `flagged` (String)
  * `curated` (String)
  * `comment` (String)
  * `wgd` (String)
  * `loh_fraction` (String)
  * `cin` (String)
  * `cin_allele_specific` (String)
  * `cin_ploidy_robust` (String)
  * `cin_allele_specific_ploidy_robust` (String)

## PureCN

### Inputs

#### Required

  * `intervals` (File, **required**)
  * `sampleID` (String, **required**)
  * `segFile` (File, **required**)
  * `vcf` (File, **required**)

#### Defaults

  * `funSegmentation` (String, default="Hclust")
  * `genome` (String, default="hg38")
  * `hardware_disk_size_GB` (Int, default=250)
  * `hardware_memory_GB` (Int, default=32)
  * `hardware_preemptible_tries` (Int, default=2)
  * `maxCopyNumber` (Int, default=8)
  * `maxPurity` (Float, default=0.99)
  * `maxSegments` (Int, default=1000)
  * `max_retries` (Int, default=0)
  * `minPurity` (Float, default=0.9)
  * `num_threads` (Int, default=1)
  * `otherArguments` (String, default="--post-optimize --model-homozygous --min-total-counts 20")
  * `purecn_docker` (String, default="markusriester/purecn:2.2.0")

### Outputs

  * `solutions_pdf` (File)
  * `chromosomes_pdf` (File)
  * `rds` (File)
  * `dnacopy` (File)
  * `variants` (File)
  * `loh` (File)
  * `genes` (File)
  * `segmentation` (File)
  * `log` (File)
  * `selected_solution` (File)
  * `local_optima_pdf` (File)
  * `table` (Array[String])
  * `purity` (String)
  * `ploidy` (String)
  * `contamination` (String)
  * `flagged` (String)
  * `curated` (String)
  * `comment` (String)
  * `wgd_table` (Array[String])
  * `wgd` (String)
  * `loh_fraction` (String)
  * `cin` (String)
  * `cin_allele_specific` (String)
  * `cin_ploidy_robust` (String)
  * `cin_allele_specific_ploidy_robust` (String)
