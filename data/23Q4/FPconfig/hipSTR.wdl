task run_hipstr {
  
  # File-related inputs
  String sample_name
  File bam_file
  File bam_index
  File ref
  File ref_index
  File bed
  
  # Hardware-related inputs
  Int? hardware_disk_size_GB = 100
  Int? hardware_memory_GB = 16
  Int? hardware_preemptible_tries = 2
  
  command {
  	HipSTR --bams "${bam_file}" \
    		--fasta "${ref}" \
            --regions "${bed}" \
            --str-vcf "${sample_name}.vcf.gz" \
            --min-reads 20 --def-stutter-model --max-str-len 150
  }
  
  output {
    File vcf = "${sample_name}.vcf.gz"
  }
  
  runtime {
    docker: "philpalmer/hipstr:latest"
    bootDiskSizeGb: 32
    disks: "local-disk ${hardware_disk_size_GB} HDD"
    memory: "${hardware_memory_GB}GB"
    cpu: 1
    preemptible: hardware_preemptible_tries
    maxRetries: 0
  }
  
}

workflow hipstr {
  
  call run_hipstr
  
  output {
    File vcf = run_hipstr.vcf
    }
}