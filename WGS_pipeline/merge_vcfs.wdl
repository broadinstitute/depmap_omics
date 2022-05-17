version 1.0

# Given a set of samples, combine segment files into a single file
# more information available at https://open-cravat.readthedocs.io/en/latest/2.-Command-line-usage.html
workflow run_merge_vcfs {
    input {
        Array[File] vcfs
        String sample_id
        String merge_mode
    }

    call merge_vcfs {
        input:
            vcfs=vcfs,
            sample_id=sample_id,
            merge_mode=merge_mode,
    }

    output {
        File merged_vcf = merge_vcfs.vcf_merged
        File merged_vcf_index = merge_vcfs.vcf_merged_index
    }
}

task merge_vcfs {
    input {
        Array[File] vcfs
        String sample_id
        String merge_mode
        String output_type="z"
    
        Int memory = 4
        Int boot_disk_size = 10
        Int num_threads = 1
        Int num_preempt = 5
        Int disk_space = 40
        String docker = "biocontainers/bcftools:1.3"
    }

    command {
        set -euo pipefail

        for vcf in "${sep=' ' vcfs}"; do
            bcftools index "$vcf"
        done
 
        bcftools merge \
            --missing-to-ref \
            --no-index \
            --threads ${num_threads} \
            --output-type ${output_type} \
            --output ${sample_id}.vcf.gz \
            --merge ${merge_mode}
            ${sep=" " vcfs}
        
        bcftools index ${sample_id}.vcf.gz
    }

    output {
        File vcf_merged="${sample_id}.vcf.gz"
        File vcf_merged_index="${sample_id}.vcf.gz.csi"
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