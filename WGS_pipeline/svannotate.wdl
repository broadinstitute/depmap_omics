version 1.0
# partially stolen from https://github.com/broadinstitute/long-read-pipelines/blob/ym_testing/wdl/SVAnnotate.wdl

workflow SVAnnotate_workflow {
    input {
        String sample_id
        File vcf
        File coding_gtf
        File noncoding_bed

        # taken from https://app.terra.bio/#workspaces/broad-firecloud-dsde-methods/GATK-Structural-Variants-Joint-Calling/data
        String standardizeVCF_docker="us.gcr.io/broad-dsde-methods/gatk-sv/sv-pipeline:2024-11-15-v1.0-488d7cb0"
        String svannotate_docker="quay.io/ymostovoy/lr-svannotate:latest"
    }

    call StandardizeVCF {
        input: 
            raw_vcf=vcf,
            sample_id=sample_id,
            caller="manta",
            contigs="gs://gcp-public-data--broad-references/hg38/v0/sv-resources/resources/v1/contig.fai",
            sv_pipeline_docker=standardizeVCF_docker,
    }

    call SVAnnotate {
        input:
            sample_id=sample_id,
            vcf=StandardizeVCF.out, 
            coding_gtf=coding_gtf, 
            noncoding_bed=noncoding_bed,
            docker=docker
    }

    output {
        File vcf_anno = SVAnnotate.vcf_anno
        File vcf_anno_tbi = SVAnnotate.vcf_anno_tbi
    }
}

task StandardizeVCF {
    input {
        File raw_vcf
        String sample_id
        String caller
        File contigs
        String sv_pipeline_docker
        Int min_svsize=50

        Int cpu = 2
        Int mem_gb = 8
        Int preemptible = 3
        Int max_retries = 0
        Int additional_disk_gb = 0
        Int boot_disk_size = 60
    }

    Int disk_size = 5 + 5*ceil(size(raw_vcf, "GB"))

    command <<<
        set -euo pipefail
        mkdir out
        
        svtk standardize --sample-names ${sample_id} --prefix ~{caller}_${sample_id} --contigs ~{contigs} --min-size ~{min_svsize} ~{raw_vcf} tmp.vcf ~{caller}
        bcftools sort tmp.vcf -Oz -o out/std.~{caller}.${sample_id}.vcf.gz
        tar czf ~{sample_id}.tar.gz -C out/ .
    >>>

    output {
        File out = "~{sample_id}.tar.gz"
    }

    runtime {
        docker: sv_pipeline_docker
        bootDiskSizeGb: boot_disk_size
        memory: "~{mem_gb} GB"
        disks: "local-disk ~{disk_size} SSD"
        preemptible: preemptible
        maxRetries: max_retries
        cpu: cpu
    }
}

task SVAnnotate {
    input {
        String sample_id
        File vcf
        String docker
        File coding_gtf="gs://gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/MANE.GRCh38.v1.2.ensembl_genomic.gtf"
        File noncoding_bed="gs://gcp-public-data--broad-references/hg38/v0/sv-resources/resources/v1/noncoding.sort.hg38.bed"

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
        docker: docker
        bootDiskSizeGb: boot_disk_size
        memory: "~{mem_gb} GB"
        disks: "local-disk ~{disk_size} SSD"
        preemptible: preemptible
        maxRetries: max_retries
        cpu: cpu
    }
}
