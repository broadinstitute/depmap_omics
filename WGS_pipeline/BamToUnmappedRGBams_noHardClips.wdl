version 1.0
## Copyright Broad Institute, 2018
## https://raw.githubusercontent.com/gatk-workflows/seq-format-conversion/3.0.0/bam-to-unmapped-bams.wdl
## 
## This WDL converts BAM  to unmapped BAMs
##
## Requirements/expectations :
## - BAM file
##
## Outputs :
## - Sorted Unmapped BAMs
##
## Cromwell version support
## - Successfully tested on v47
## - Does not work on versions < v23 due to output syntax
##
## Runtime parameters are optimized for Broad's Google Cloud Platform implementation. 
## For program versions, see docker containers. 
##
## LICENSING : 
## This script is released under the WDL source code license (BSD-3) (see LICENSE in 
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may 
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker 
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.

# WORKFLOW DEFINITION
workflow BamToUnmappedBams {
  input {
    File input_bam

    Int additional_disk_size = 20
    String gatk_docker = "broadinstitute/gatk:4.2.6.1"
    String gatk_path = "/gatk/gatk"
  }
    Float input_size = size(input_bam, "GB")
    
  call RevertSam {
    input:
      input_bam = input_bam,
      disk_size = ceil(input_size * 3) + additional_disk_size,
      docker = gatk_docker,
      gatk_path = gatk_path
  }

  scatter (unmapped_bam in RevertSam.unmapped_bams) {
    String output_basename = basename(unmapped_bam, ".coord.sorted.unmapped.bam")
    Float unmapped_bam_size = size(unmapped_bam, "GB")

    call SortSam {
      input:
        input_bam = unmapped_bam,
        sorted_bam_name = output_basename + ".unmapped.bam",
        disk_size = ceil(unmapped_bam_size * 6) + additional_disk_size,
        docker = gatk_docker,
        gatk_path = gatk_path
    }
  }

  output {
    Array[File] output_bams = SortSam.sorted_bam
  }
}

task RevertSam {
  input {
    #Command parameters
    File input_bam
    String gatk_path

    #Runtime parameters
    Int disk_size
    String docker
    Int machine_mem_gb = 2
    Int preemptible_attempts = 3
  }
    Int command_mem_gb = machine_mem_gb - 1    ####Needs to occur after machine_mem_gb is set 

  command { 
 
    ~{gatk_path} --java-options "-Xmx~{command_mem_gb}g" \
    RevertSam \
    --INPUT ~{input_bam} \
    --OUTPUT ./ \
    --OUTPUT_BY_READGROUP true \
    --VALIDATION_STRINGENCY LENIENT \
    --ATTRIBUTE_TO_CLEAR FT \
    --ATTRIBUTE_TO_CLEAR CO \
    --SORT_ORDER coordinate \
    --RESTORE_HARDCLIPS false
  }
  runtime {
    docker: docker
    disks: "local-disk " + disk_size + " SSD"
    memory: machine_mem_gb + " GB"
    preemptible: preemptible_attempts
  }
  output {
    Array[File] unmapped_bams = glob("*.bam")
  }
}

task SortSam {
  input {
    #Command parameters
    File input_bam
    String sorted_bam_name
    #Runtime parameters
    String gatk_path
    Int disk_size
    String docker
    Int machine_mem_gb = 4
    Int preemptible_attempts = 3
  }
    Int command_mem_gb = machine_mem_gb - 1    ####Needs to occur after machine_mem_gb is set 

  command {
    ~{gatk_path} --java-options "-Xmx~{command_mem_gb}g" \
    SortSam \
    --INPUT ~{input_bam} \
    --OUTPUT ~{sorted_bam_name} \
    --SORT_ORDER queryname \
    --MAX_RECORDS_IN_RAM 1000000
  }
  runtime {
    docker: docker
    disks: "local-disk " + disk_size + " SSD"
    memory: machine_mem_gb + " GB"
    preemptible: preemptible_attempts
  }
  output {
    File sorted_bam = "~{sorted_bam_name}"
  }
}
