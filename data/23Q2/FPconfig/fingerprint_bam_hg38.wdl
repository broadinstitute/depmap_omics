task run_fingerprint {
  
  # File-related inputs
  File bam
  File bam_index
  String output_name
  File haplotype_map
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  
  # Method configuration inputs
  Int? contamination
  String? sample_alias
  Int? locus_max_reads
  
  # Hardware-related inputs
  Int? hardware_disk_size_GB = 300
  Int? hardware_memory_GB = 16
  Int? hardware_preemptible_tries = 2
  
  command {
  	java "-Xmx${hardware_memory_GB}g" -jar /usr/picard/picard.jar ExtractFingerprint \
      --INPUT "${bam}" \
      --OUTPUT "${output_name}.vcf" \
      --HAPLOTYPE_MAP "${haplotype_map}" \
      --REFERENCE_SEQUENCE "${ref_fasta}" \
      ${"--CONTAMINATION " + contamination} \
      ${"--SAMPLE_ALIAS " + sample_alias} \
      ${"--LOCUS_MAX_READS " + locus_max_reads}
  }
  
  output {
    File out = "${output_name}.vcf"
  }
  
  runtime {
    docker: "broadinstitute/picard:2.25.0"
    bootDiskSizeGb: 32
    disks: "local-disk ${hardware_disk_size_GB} HDD"
    memory: "${hardware_memory_GB}GB"
    cpu: 1
    preemptible: hardware_preemptible_tries
    maxRetries: 0
  }
  
}

workflow fingerprint {
  call run_fingerprint
}