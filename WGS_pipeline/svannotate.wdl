version 1.0
# partially stolen from https://github.com/broadinstitute/long-read-pipelines/blob/ym_testing/wdl/SVAnnotate.wdl

workflow SVAnnotate_workflow {
    input {
        String sample_id
        File vcf
        File coding_gtf
        File noncoding_bed
    }

    call SVAnnotate {
        input:
            sample_id=sample_id,
            vcf=vcf, 
            coding_gtf=coding_gtf, 
            noncoding_bed=noncoding_bed
    }

    output {
        File vcf_anno = Standardize.vcf_anno
        File vcf_anno_tbi = Standardize.vcf_anno_tbi
    }
}

task SVAnnotate {
    input {
        File vcf
        File coding_gtf
        File noncoding_bed
    }

    Int disk_size = 5 + 5*ceil(size(vcf, "GB"))

    command <<<
        set -euxo pipefail

        tabix ~{vcf}

        gatk SVAnnotate -V ~{vcf} --non-coding-bed ~{noncoding_bed} --protein-coding-gtf ~{coding_gtf} -O ~{sample_id}.anno.vcf
        bcftools view -Oz ~{sample_id}.anno.vcf > ~{sample_id}.anno.vcf.gz
        tabix ~{sample_id}.anno.vcf.gz
    >>>

    output {
        File vcf_anno = "~{sample_id}.anno.vcf.gz"
        File vcf_anno_tbi = "~{sample_id}.anno.vcf.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             2,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "quay.io/ymostovoy/lr-svannotate:latest"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}
