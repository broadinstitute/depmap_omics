task run_crosscheck {
  
  # File-related inputs
  Array[File] vcf_input_array
  Array[File] vcf_second_input_array
  String output_name
  File haplotype_map
  
  # Method configuration inputs
  Int? genotyping_error_rate
  String? tumor_aware_results = "False"
  Int? loss_of_het_rate
  
  # Hardware-related inputs
  Int? hardware_disk_size_GB = 128
  Int? hardware_memory_GB = 16
  Int? hardware_preemptible_tries = 2
  
  command {
  	java "-Xmx${hardware_memory_GB}g" -jar /usr/picard/picard.jar CrosscheckFingerprints \
		-I ${sep=" -I " vcf_input_array} \
        -SI ${sep=" -SI " vcf_second_input_array} \
        -H "${haplotype_map}" \
		--CALCULATE_TUMOR_AWARE_RESULTS "${tumor_aware_results}" \
		--EXIT_CODE_WHEN_MISMATCH 0 \
        --CROSSCHECK_BY FILE \
        --CROSSCHECK_MODE CHECK_ALL_OTHERS \
		-O "${output_name}_crosscheck" \
      	${"--GENOTYPING_ERROR_RATE " + genotyping_error_rate} \
      	${"--LOSS_OF_HET_RATE " + loss_of_het_rate}
  }
  
  output {
    File out = "${output_name}_crosscheck"
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

workflow crosscheck {
  call run_crosscheck
}