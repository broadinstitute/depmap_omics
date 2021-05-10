task run_fingerprint {
  
  # File-related inputs
  File bam
  File bam_index
  String output_name
  File haplotype_map
  File bam_ref_fasta
  File bam_ref_fasta_index
  File bam_ref_dict
  File liftover_chain
  File liftover_ref_fasta
  File liftover_ref_fasta_index
  File liftover_ref_dict
  
  
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
      --OUTPUT "${output_name}_fingerprint.vcf" \
      --HAPLOTYPE_MAP "${haplotype_map}" \
      --REFERENCE_SEQUENCE "${bam_ref_fasta}" \
      ${"--CONTAMINATION " + contamination} \
      ${"--SAMPLE_ALIAS " + sample_alias} \
      ${"--LOCUS_MAX_READS " + locus_max_reads}
      
    java "-Xmx${hardware_memory_GB}g" -jar /usr/picard/picard.jar LiftoverVcf \
      --INPUT "${output_name}_fingerprint.vcf" \
      --OUTPUT "${output_name}.vcf" \
      --REFERENCE_SEQUENCE "${liftover_ref_fasta}" \
      --CHAIN "${liftover_chain}" \
      --REJECT reject.vcf \
      --RECOVER_SWAPPED_REF_ALT
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