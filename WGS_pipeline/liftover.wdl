
task liftover{

	String? gatk_loc = "gatk"
    String? docker = "broadinstitute/gatk:latest"
    
    File input_vcf 
    String? sample_name = "sample"
    File ref_fasta
    File ref_dict
    File ref_fasta_idx
    File chain_file
    
    String? memory="16"
    String? disk_space="200"
    Int? num_threads=1
    Int? num_preempt=2
    Int? max_records_in_ram=100000

	command {
    	java -Xmx12g -jar /root/gatk.jar LiftoverVcf \
        	-I ${input_vcf} \
            --MAX_RECORDS_IN_RAM ${max_records_in_ram}\
            -O ${sample_name}_liftover.vcf \
            -C ${chain_file} \
            --REJECT ${sample_name}_rejected_variants.vcf \
            -R ${ref_fasta}       
    }
    
    output {
        File lifted_vcf = "${sample_name}_liftover.vcf"
        File rejected_mutations_fromliftover = "${sample_name}_rejected_variants.vcf"
    }

    runtime {
        docker: "${docker}"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} SSD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    meta {
        author: "Jeremie Kalfon"
    }  
}

workflow liftover_workflows {
    call liftover
}
