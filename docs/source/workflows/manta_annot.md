
## run_manta_annotator

### Inputs

#### Required

  * `sv` (File, **required**)

#### Defaults

  * `exon_annot` (File, default="gs://ccleparams/hg38_ensembl_exonlocations_formatted.txt")
  * `manta_annotator.cores` (Int, default=8)
  * `manta_annotator.disk_size` (Int, default=30)
  * `manta_annotator.docker_image` (String, default="jkobject/manta_annot")
  * `manta_annotator.mem_size` (Int, default=4)
  * `manta_annotator.preemptible_tries` (Int, default=3)

### Outputs

  * `somatic_annotated_sv` (File)
  * `filtered_annotated_sv` (File)
  * `dropped` (File)

## manta_annotator

### Inputs

#### Required

  * `exon_annot` (File, **required**)
  * `sv` (File, **required**)

#### Defaults

  * `cores` (Int, default=8)
  * `disk_size` (Int, default=30)
  * `docker_image` (String, default="jkobject/manta_annot")
  * `mem_size` (Int, default=4)
  * `preemptible_tries` (Int, default=3)

### Outputs

  * `somatic_annotated_sv` (File)
  * `filtered_annotated_sv` (File)
  * `dropped` (File)
