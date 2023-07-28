
## MantaSomaticSV

### Inputs

#### Required

  * `config_manta` (String, **required**)
  * `is_cram` (Boolean, **required**)
  * `manta_docker` (String, **required**)
  * `ref_fasta` (File, **required**)
  * `ref_fasta_index` (File, **required**)
  * `sample_name` (String, **required**)
  * `tumor_bam` (File, **required**)
  * `tumor_bam_index` (File, **required**)

#### Optional

  * `interval_list` (File?)
  * `normal_bam` (File?)
  * `normal_bam_index` (File?)
  * `Manta.cpu_num` (Int?)
  * `Manta.disk_size` (Int?)
  * `Manta.mem_size` (Int?)
  * `Manta.preemptible_attempts` (Int?)

#### Defaults

  * `is_exome` (Boolean, default=defined(interval_list))
  * `ConvertToBedTabix.output_bed` (String, default="interval.bed")

### Outputs

  * `germline_sv_vcf` (File)
  * `germline_sv_vcf_index` (File)
  * `somatic_sv_vcf` (File)
  * `somatic_sv_vcf_index` (File)
  * `candidate_sv_vcf` (File)
  * `candidate_sv_vcf_index` (File)
  * `candidate_indel_vcf` (File)
  * `candidate_indel_vcf_index` (File)

## Manta

### Inputs

#### Required

  * `config_manta` (String, **required**)
  * `is_cram` (Boolean, **required**)
  * `manta_docker` (String, **required**)
  * `ref_fasta` (File, **required**)
  * `ref_fasta_index` (File, **required**)
  * `sample_name` (String, **required**)
  * `tumor_bam` (File, **required**)
  * `tumor_bam_index` (File, **required**)

#### Optional

  * `cpu_num` (Int?)
  * `disk_size` (Int?)
  * `interval_bed` (File?)
  * `interval_bed_index` (File?)
  * `mem_size` (Int?)
  * `normal_bam` (File?)
  * `normal_bam_index` (File?)
  * `preemptible_attempts` (Int?)

### Outputs

  * `germline_sv_vcf` (File)
  * `germline_sv_vcf_index` (File)
  * `somatic_sv_vcf` (File)
  * `somatic_sv_vcf_index` (File)
  * `candidate_sv_vcf` (File)
  * `candidate_sv_vcf_index` (File)
  * `candidate_indel_vcf` (File)
  * `candidate_indel_vcf_index` (File)

## ConvertToBedTabix

### Inputs

#### Optional

  * `interval_list` (File?)

#### Defaults

  * `output_bed` (String, default="interval.bed")

### Outputs

  * `out_interval` (File)
  * `out_interval_index` (File)
