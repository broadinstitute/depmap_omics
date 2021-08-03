import "https://raw.githubusercontent.com/broadinstitute/ccle_processing/simone-sandbox/RNA_pipeline/rsem_aggregate_results.wdl" as rsem_aggregate_results
import "https://raw.githubusercontent.com/broadinstitute/ccle_processing/simone-sandbox/RNA_pipeline/Aggregate_Fusion_Calls.wdl" as Aggregate_Fusion_Calls


workflow RNA_aggregate {

  #rsem_aggregate_results

  call rsem_aggregate_results.rsem_aggregate_results

  call Aggregate_Fusion_Calls.aggregate_set_files

}