import "https://raw.githubusercontent.com/broadinstitute/ccle_processing/reorganization/WGS_pipeline/Aggregate_CN_seg_files.wdl" as Aggregate_CN_seg_files 

workflow WGS_aggregate {
  call Aggregate_CN_seg_files.aggregate_CN_segments_wrkflw as aggregate_CN_segments_wrkflw
}