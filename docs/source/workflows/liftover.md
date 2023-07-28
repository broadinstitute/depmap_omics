
## liftover_workflows

### Inputs

#### Required

  * `liftover.chain_file` (File, **required**)
  * `liftover.input_vcf` (File, **required**)
  * `liftover.ref_dict` (File, **required**)
  * `liftover.ref_fasta` (File, **required**)
  * `liftover.ref_fasta_idx` (File, **required**)

#### Optional

  * `liftover.disk_space` (String?)
  * `liftover.docker` (String?)
  * `liftover.gatk_loc` (String?)
  * `liftover.max_records_in_ram` (Int?)
  * `liftover.memory` (String?)
  * `liftover.num_preempt` (Int?)
  * `liftover.num_threads` (Int?)
  * `liftover.sample_name` (String?)

### Outputs

  * `liftover.lifted_vcf` (File)
  * `liftover.rejected_mutations_fromliftover` (File)

## liftover

author
: Jeremie Kalfon

### Inputs

#### Required

  * `chain_file` (File, **required**)
  * `input_vcf` (File, **required**)
  * `ref_dict` (File, **required**)
  * `ref_fasta` (File, **required**)
  * `ref_fasta_idx` (File, **required**)

#### Optional

  * `disk_space` (String?)
  * `docker` (String?)
  * `gatk_loc` (String?)
  * `max_records_in_ram` (Int?)
  * `memory` (String?)
  * `num_preempt` (Int?)
  * `num_threads` (Int?)
  * `sample_name` (String?)

### Outputs

  * `lifted_vcf` (File)
  * `rejected_mutations_fromliftover` (File)
