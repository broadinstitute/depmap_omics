
## trinity_cleaned

### Inputs

#### Required

  * `ctat_genome_lib_build_dir_files` (Array[File], **required**)
  * `fastq1` (File, **required**)
  * `fastq2` (File, **required**)
  * `prefix` (String, **required**)
  * `ref_genome_fa_star_idx_files` (Array[File], **required**)

#### Defaults

  * `StarFusion.disk_space` (Int, default=500)
  * `StarFusion.docker` (String, default="trinityctat/starfusion:1.7.0")
  * `StarFusion.memory` (Int, default=64)
  * `StarFusion.num_preempt` (Int, default=2)
  * `StarFusion.num_threads` (Int, default=8)

### Outputs

  * `fusion_predictions` (File)
  * `fusion_predictions_abridged` (File)

## StarFusion

### Inputs

#### Required

  * `ctat_genome_lib_build_dir_files` (Array[File], **required**)
  * `left_fastq` (File, **required**)
  * `prefix` (String, **required**)
  * `ref_genome_fa_star_idx_files` (Array[File], **required**)
  * `right_fastq` (File, **required**)

#### Defaults

  * `disk_space` (Int, default=500)
  * `docker` (String, default="trinityctat/starfusion:1.7.0")
  * `memory` (Int, default=64)
  * `num_preempt` (Int, default=2)
  * `num_threads` (Int, default=8)

### Outputs

  * `fusion_predictions` (File)
  * `fusion_predictions_abridged` (File)
