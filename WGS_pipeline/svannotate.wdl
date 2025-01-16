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
        String sample_id
        File vcf
        File coding_gtf
        File noncoding_bed

        Int cpu = 2
        Int mem_gb = 8
        Int preemptible = 3
        Int max_retries = 0
        Int additional_disk_gb = 0
        Int boot_disk_size = 60
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

    runtime {
        docker: "quay.io/ymostovoy/lr-svannotate:latest"
        bootDiskSizeGb: boot_disk_size
        memory: "~{mem_gb} GB"
        disks: "local-disk ~{disk_size} SSD"
        preemptible: preemptible
        maxRetries: max_retries
        cpu: cpu
    }
}
