
## trinity_cleaned

### Inputs

#### Required

  * `prefix` (String, **required**)
  * `StarFusion.ctat_genome_lib_build_dir_files` (Array[File], **required**)
  * `StarFusion.disk_space` (Int, **required**)
  * `StarFusion.docker` (String, **required**)
  * `StarFusion.left_fastq` (File, **required**)
  * `StarFusion.memory` (Int, **required**)
  * `StarFusion.num_preempt` (Int, **required**)
  * `StarFusion.num_threads` (Int, **required**)
  * `StarFusion.ref_genome_fa_star_idx_files` (Array[File], **required**)
  * `StarFusion.right_fastq` (File, **required**)

### Outputs

  * `StarFusion.fusion_predictions` (File)
  * `StarFusion.fusion_predictions_abridged` (File)

## StarFusion

### Inputs

#### Required

  * `ctat_genome_lib_build_dir_files` (Array[File], **required**)
  * `disk_space` (Int, **required**)
  * `docker` (String, **required**)
  * `left_fastq` (File, **required**)
  * `memory` (Int, **required**)
  * `num_preempt` (Int, **required**)
  * `num_threads` (Int, **required**)
  * `prefix` (String, **required**)
  * `ref_genome_fa_star_idx_files` (Array[File], **required**)
  * `right_fastq` (File, **required**)

### Outputs

  * `fusion_predictions` (File)
  * `fusion_predictions_abridged` (File)
