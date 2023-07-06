
## star_workflow

### Inputs

#### Required

  * `star.disk_space` (Int, **required**)
  * `star.fastq1` (File, **required**)
  * `star.memory` (Int, **required**)
  * `star.num_preempt` (Int, **required**)
  * `star.num_threads` (Int, **required**)
  * `star.prefix` (String, **required**)
  * `star.star_index` (File, **required**)

#### Optional

  * `star.alignIntronMax` (Int?)
  * `star.alignIntronMin` (Int?)
  * `star.alignMatesGapMax` (Int?)
  * `star.alignSJDBoverhangMin` (Int?)
  * `star.alignSJoverhangMin` (Int?)
  * `star.alignSoftClipAtReferenceEnds` (String?)
  * `star.chimJunctionOverhangMin` (Int?)
  * `star.chimMainSegmentMultNmax` (Int?)
  * `star.chimOutJunctionFormat` (Int?)
  * `star.chimOutType` (String?)
  * `star.chimSegmentMin` (Int?)
  * `star.fastq2` (File?)
  * `star.limitSjdbInsertNsj` (Int?)
  * `star.outFilterIntronMotifs` (String?)
  * `star.outFilterMatchNminOverLread` (Float?)
  * `star.outFilterMismatchNmax` (Int?)
  * `star.outFilterMismatchNoverLmax` (Float?)
  * `star.outFilterMultimapNmax` (Int?)
  * `star.outFilterScoreMinOverLread` (Float?)
  * `star.outFilterType` (String?)
  * `star.outSAMattrRGline` (String?)
  * `star.outSAMattributes` (String?)
  * `star.outSAMstrandField` (String?)
  * `star.quantMode` (String?)
  * `star.sjdbFileChrStartEnd` (File?)
  * `star.varVCFfile` (File?)
  * `star.waspOutputMode` (String?)

### Outputs

  * `star.bam_file` (File)
  * `star.bam_index` (File)
  * `star.transcriptome_bam` (File)
  * `star.chimeric_junctions` (File)
  * `star.chimeric_bam_file` (File)
  * `star.chimeric_bam_index` (File)
  * `star.read_counts` (File)
  * `star.junctions` (File)
  * `star.junctions_pass1` (File)
  * `star.logs` (Array[File])

## star

author
: Francois Aguet

### Inputs

#### Required

  * `disk_space` (Int, **required**)
  * `fastq1` (File, **required**)
  * `memory` (Int, **required**)
  * `num_preempt` (Int, **required**)
  * `num_threads` (Int, **required**)
  * `prefix` (String, **required**)
  * `star_index` (File, **required**)

#### Optional

  * `alignIntronMax` (Int?)
  * `alignIntronMin` (Int?)
  * `alignMatesGapMax` (Int?)
  * `alignSJDBoverhangMin` (Int?)
  * `alignSJoverhangMin` (Int?)
  * `alignSoftClipAtReferenceEnds` (String?)
  * `chimJunctionOverhangMin` (Int?)
  * `chimMainSegmentMultNmax` (Int?)
  * `chimOutJunctionFormat` (Int?)
  * `chimOutType` (String?)
  * `chimSegmentMin` (Int?)
  * `fastq2` (File?)
  * `limitSjdbInsertNsj` (Int?)
  * `outFilterIntronMotifs` (String?)
  * `outFilterMatchNminOverLread` (Float?)
  * `outFilterMismatchNmax` (Int?)
  * `outFilterMismatchNoverLmax` (Float?)
  * `outFilterMultimapNmax` (Int?)
  * `outFilterScoreMinOverLread` (Float?)
  * `outFilterType` (String?)
  * `outSAMattrRGline` (String?)
  * `outSAMattributes` (String?)
  * `outSAMstrandField` (String?)
  * `quantMode` (String?)
  * `sjdbFileChrStartEnd` (File?)
  * `varVCFfile` (File?)
  * `waspOutputMode` (String?)

### Outputs

  * `bam_file` (File)
  * `bam_index` (File)
  * `transcriptome_bam` (File)
  * `chimeric_junctions` (File)
  * `chimeric_bam_file` (File)
  * `chimeric_bam_index` (File)
  * `read_counts` (File)
  * `junctions` (File)
  * `junctions_pass1` (File)
  * `logs` (Array[File])
