
## aggregate_CN_segments_wrkflw

### Inputs

#### Required

  * `aggregate_CN_segments.aggregate_seg_files_script` (File, **required**)
  * `aggregate_CN_segments.disk_space` (Int, **required**)
  * `aggregate_CN_segments.memory` (Int, **required**)
  * `aggregate_CN_segments.num_preempt` (Int, **required**)
  * `aggregate_CN_segments.sample_seg_files` (Array[File], **required**)
  * `aggregate_CN_segments.sample_set_id` (String, **required**)

### Outputs

  * `aggregate_CN_segments.combined_cn_file` (File)

## aggregate_CN_segments

author
: Guillaume Kugener

### Inputs

#### Required

  * `aggregate_seg_files_script` (File, **required**)
  * `disk_space` (Int, **required**)
  * `memory` (Int, **required**)
  * `num_preempt` (Int, **required**)
  * `sample_seg_files` (Array[File], **required**)
  * `sample_set_id` (String, **required**)

### Outputs

  * `combined_cn_file` (File)
