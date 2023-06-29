
## star_workflow

### Inputs

#### Required

  * `fastq1` (File, **required**)
  * `prefix` (String, **required**)
  * `star_index` (File, **required**)

#### Optional

  * `fastq2` (File?)
  * `star.sjdbFileChrStartEnd` (File?)
  * `star.varVCFfile` (File?)
  * `star.waspOutputMode` (String?)

#### Defaults

  * `star.alignIntronMax` (Int, default=1000000)
  * `star.alignIntronMin` (Int, default=20)
  * `star.alignMatesGapMax` (Int, default=1000000)
  * `star.alignSJDBoverhangMin` (Int, default=1)
  * `star.alignSJoverhangMin` (Int, default=8)
  * `star.alignSoftClipAtReferenceEnds` (String, default="Yes")
  * `star.chimJunctionOverhangMin` (Int, default=15)
  * `star.chimMainSegmentMultNmax` (Int, default=1)
  * `star.chimOutJunctionFormat` (Int, default=1)
  * `star.chimOutType` (String, default="Junctions WithinBAM SoftClip")
  * `star.chimSegmentMin` (Int, default=15)
  * `star.disk_space` (Int, default=500)
  * `star.limitSjdbInsertNsj` (Int, default=1200000)
  * `star.memory` (Int, default=64)
  * `star.num_preempt` (Int, default=1)
  * `star.num_threads` (Int, default=16)
  * `star.outFilterIntronMotifs` (String, default="None")
  * `star.outFilterMatchNminOverLread` (Float, default=0.33)
  * `star.outFilterMismatchNmax` (Int, default=999)
  * `star.outFilterMismatchNoverLmax` (Float, default=0.1)
  * `star.outFilterMultimapNmax` (Int, default=20)
  * `star.outFilterScoreMinOverLread` (Float, default=0.33)
  * `star.outFilterType` (String, default="BySJout")
  * `star.outSAMattrRGline` (String, default="ID:rg1 SM:sm1")
  * `star.outSAMattributes` (String, default="NH HI AS nM NM ch")
  * `star.outSAMstrandField` (String, default="intronMotif")
  * `star.quantMode` (String, default="TranscriptomeSAM GeneCounts")

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

## star

author
: Francois Aguet

### Inputs

#### Required

  * `fastq1` (File, **required**)
  * `prefix` (String, **required**)
  * `star_index` (File, **required**)

#### Optional

  * `fastq2` (File?)
  * `sjdbFileChrStartEnd` (File?)
  * `varVCFfile` (File?)
  * `waspOutputMode` (String?)

#### Defaults

  * `alignIntronMax` (Int, default=1000000)
  * `alignIntronMin` (Int, default=20)
  * `alignMatesGapMax` (Int, default=1000000)
  * `alignSJDBoverhangMin` (Int, default=1)
  * `alignSJoverhangMin` (Int, default=8)
  * `alignSoftClipAtReferenceEnds` (String, default="Yes")
  * `chimJunctionOverhangMin` (Int, default=15)
  * `chimMainSegmentMultNmax` (Int, default=1)
  * `chimOutJunctionFormat` (Int, default=1)
  * `chimOutType` (String, default="Junctions WithinBAM SoftClip")
  * `chimSegmentMin` (Int, default=15)
  * `disk_space` (Int, default=500)
  * `limitSjdbInsertNsj` (Int, default=1200000)
  * `memory` (Int, default=64)
  * `num_preempt` (Int, default=1)
  * `num_threads` (Int, default=16)
  * `outFilterIntronMotifs` (String, default="None")
  * `outFilterMatchNminOverLread` (Float, default=0.33)
  * `outFilterMismatchNmax` (Int, default=999)
  * `outFilterMismatchNoverLmax` (Float, default=0.1)
  * `outFilterMultimapNmax` (Int, default=20)
  * `outFilterScoreMinOverLread` (Float, default=0.33)
  * `outFilterType` (String, default="BySJout")
  * `outSAMattrRGline` (String, default="ID:rg1 SM:sm1")
  * `outSAMattributes` (String, default="NH HI AS nM NM ch")
  * `outSAMstrandField` (String, default="intronMotif")
  * `quantMode` (String, default="TranscriptomeSAM GeneCounts")

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
