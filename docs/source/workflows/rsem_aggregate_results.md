
## rsem_aggregate_results_workflow

### Inputs

#### Required

  * `rsem_aggregate_results.disk_space` (Int, **required**)
  * `rsem_aggregate_results.memory` (Int, **required**)
  * `rsem_aggregate_results.num_preempt` (Int, **required**)
  * `rsem_aggregate_results.num_threads` (Int, **required**)
  * `rsem_aggregate_results.prefix` (String, **required**)
  * `rsem_aggregate_results.rsem_genes` (Array[File], **required**)
  * `rsem_aggregate_results.rsem_isoforms` (Array[File], **required**)

### Outputs

  * `rsem_aggregate_results.transcripts_tpm` (File)
  * `rsem_aggregate_results.transcripts_isopct` (File)
  * `rsem_aggregate_results.transcripts_expected_count` (File)
  * `rsem_aggregate_results.genes_tpm` (File)
  * `rsem_aggregate_results.genes_expected_count` (File)

## rsem_aggregate_results

author
: Francois Aguet

### Inputs

#### Required

  * `disk_space` (Int, **required**)
  * `memory` (Int, **required**)
  * `num_preempt` (Int, **required**)
  * `num_threads` (Int, **required**)
  * `prefix` (String, **required**)
  * `rsem_genes` (Array[File], **required**)
  * `rsem_isoforms` (Array[File], **required**)

### Outputs

  * `transcripts_tpm` (File)
  * `transcripts_isopct` (File)
  * `transcripts_expected_count` (File)
  * `genes_tpm` (File)
  * `genes_expected_count` (File)
