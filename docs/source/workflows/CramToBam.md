
## CramToBam

author
: Ted Brookings

email
: tbrookin@broadinstitute.org

### Inputs

#### Required

  * `cram_file` (File, **required**): .bam or .cram file to search for SVs. bams are preferable, crams will be converted to bams.
  * `reference_fasta` (File, **required**): .fasta file with reference used to align bam or cram file
  * `samtools_docker` (String, **required**)

#### Optional

  * `reference_index` (File?): [optional] reference index file. If omitted, the WDL will look for an index by appending .fai to the .fasta file
  * `runtime_attr_override` (RuntimeAttr?)

#### Defaults

  * `requester_pays` (Boolean, default=false)

### Outputs

  * `bam_file` (File)
  * `bam_index` (File)

## RunCramToBam

### Inputs

#### Required

  * `cram_file` (File, **required**); **localization_optional**: false
  * `reference_fasta` (File, **required**)
  * `samtools_docker` (String, **required**)

#### Optional

  * `reference_index` (File?)
  * `runtime_attr_override` (RuntimeAttr?)

### Outputs

  * `bam_file` (File)
  * `bam_index` (File)

## RunCramToBamRequesterPays

### Inputs

#### Required

  * `cram_file` (File, **required**)
  * `reference_fasta` (File, **required**)
  * `samtools_docker` (String, **required**)

#### Optional

  * `reference_index` (File?)
  * `runtime_attr_override` (RuntimeAttr?)

### Outputs

  * `bam_file` (File)
  * `bam_index` (File)
