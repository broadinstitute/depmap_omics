version 1.0

# Given a set of samples, combine segment files into a single file
# more information available at https://open-cravat.readthedocs.io/en/latest/2.-Command-line-usage.html
import "filter_to_maf.vcf" as filter_to_maf

workflow run_merge_mutect2_replicates {
    input {
        File vcf
        String sample_id 
    }
    
    call merge_replicates {
        input:
            vcf_file=vcf,
            sample_id=sample_id
    }
    call filter_to_maf {
        input:

    }
    output {
        File merged_vcf=merge_replicates.merged_vcf
    }
}

task merge_replicates {
    input {
        File vcf_file
        String sample_id

        Int memory = 2
        Int boot_disk_size = 10
        Int num_threads = 1
        Int num_preempt = 5
        Int disk_space = 10
        String docker = "python"
    }

    command<<<
    mv ${vcf_file} torun.vcf.gz

    >>>

    output {
        File merged_vcf="${sample_id}_fixedcolumn.vcf.gz"
    }

    runtime {
        docker: docker
        bootDiskSizeGb: "${boot_disk_size}"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    meta {
        author: "Jeremie Kalfon"
    }
}