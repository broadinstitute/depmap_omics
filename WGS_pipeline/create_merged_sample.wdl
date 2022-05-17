version 1.0

# Given a set of samples, combine segment files into a single file
# more information available at https://open-cravat.readthedocs.io/en/latest/2.-Command-line-usage.html
workflow run_create_merged_sample {
    input {
        File vcf # a multi sample vcf file from "bcftools merge"
        String sample_id # the name of the merged sample
    }
    
    call create_merged_sample {
        input:
        vcf_file=vcf,
        sample_id=sample_id
    }
    output {
        File vcf_updated=create_merged_sample.vcf_updated
    }
}

task create_merged_sample {
    input {
        File vcf_file
        String sample_id
        String rna_sample_name="none"

        Int memory = 4
        Int boot_disk_size = 10
        Int num_threads = 1
        Int num_preempt = 5
        Int disk_space = 50
        Int size_tocheck = 100
        String docker = "python"
    }

    command {
        git clone https://github.com/broadinstitute/depmap_omics.git
        pip install bgzip

        python depmap_omics/WGS_pipeline/create_merged_sample.py ${vcf_file} ${sample_id} ${rna_sample_name}
    }

    output {
        File vcf_updated="${sample_id}_multi.vcf.gz"
    }

    runtime {
        docker: docker
        bootDiskSizeGb: "${boot_disk_size}"
        memory: "${memory} GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    meta {
        author: "Jeremie Kalfon"
    }
}