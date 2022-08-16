version 1.0

# Given a set of samples, combine segment files into a single file
# more information available at https://open-cravat.readthedocs.io/en/latest/2.-Command-line-usage.html
workflow run_merge_vcfs {
    input {
        Array[File] vcfs
        String sample_id
    }

    call merge_vcfs {
        input:
            vcfs=vcfs,
            sample_id=sample_id
    }
}

task merge_vcfs {
    input {
        Array[File] vcfs
        String sample_id
    
        Int memory = 4
        Int boot_disk_size = 10
        Int num_threads = 1
        Int num_preempt = 5
        Int disk_space = 40
        String docker = "biocontainers/bcftools:1.3"
    }

    command {
        set -euo pipefail

        for vcf in "${vcfs}"; do
            bcftools index "$vcf"
        done
 
        bcftools merge_vcfs 
    }

    output {
        File vcf_fixedploid="${sample_id}.vcf.gz"
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