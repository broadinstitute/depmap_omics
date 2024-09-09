version 1.0

#It is getting a merged vcf file as input (CHROM POS ... SAMPLE1 SAMPLE2).
#It tries to add an additional sample representing the merged SAMPLE fields of the vcf.
workflow run_create_merged_sample {
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
    
    call create_merged_sample {
        input:
        vcf_file=merge_vcfs.vcf_merged,
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
        disks: "local-disk ${disk_space} SSD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    meta {
        author: "Jeremie Kalfon"
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
        disks: "local-disk ${disk_space} SSD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    meta {
        author: "Jeremie Kalfon"
    }
}
