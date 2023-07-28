version 1.0

workflow run_RemoveFiltered {
    input {
        String sample_id
        File input_vcf
    }

    call RemoveFiltered {
        input:
            input_vcf=input_vcf,
            sample_id=sample_id
    }

    output {
        File output_vcf=RemoveFiltered.output_vcf
        File output_vcf_idx=RemoveFiltered.output_vcf_idx
    }
}

task RemoveFiltered {
    input {
        File input_vcf
        String sample_id
        String bcftools_exclude_string = 'FILTER~"weak_evidence"||FILTER~"map_qual"||FILTER~"strand_bias"||FILTER~"slippage"||FILTER~"clustered_events"||FILTER~"base_qual"'

        String docker_image="dceoy/bcftools"
        Int preemptible=3
        Int boot_disk_size=10
        Int cpu = 2
        Int mem = 2
    }

    command {
        bcftools view --exclude '~{bcftools_exclude_string}' --output-type b -o ~{sample_id}.filtered.vcf.gz --threads ~{cpu} ~{input_vcf} && \
        bcftools index ~{sample_id}.filtered.vcf.gz --threads ~{cpu}
    }

    runtime {
        disks: "local-disk 20 HDD"
        memory: "~{mem} GB"
        cpu: cpu
        preemptible: preemptible
        bootDiskSizeGb: boot_disk_size
        docker: docker_image
    }

    output {
        File output_vcf = "~{sample_id}.filtered.vcf.gz"
        File output_vcf_idx = "~{sample_id}.filtered.vcf.gz.csi"
        
    }
}