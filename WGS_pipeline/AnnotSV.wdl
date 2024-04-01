version 1.0
# modified from https://github.com/jmonlong/variant_annotation_wf/blob/main/wdl/workflow.wdl

workflow AnnotSV_workflow {
    input {
        File input_vcf
        String sample_id
        File annotsv_db_tar_gz = "gs://ccleparams/AnnotSV_Annotations_Human_3.4.tar.gz"
    }

    call annotate_sv_annotsv as annotate_sv_annotsv{
        input:
            input_vcf = input_vcf,
            annotsv_db_tar_gz = annotsv_db_tar_gz,
            sample_id = sample_id,
    }

    output {
        File annotated_svs = annotate_sv_annotsv.tsv
        File annotsv_log = annotate_sv_annotsv.log
    }
}


task annotate_sv_annotsv {
    input {
        File input_vcf
        File annotsv_db_tar_gz = "gs://ccleparams/AnnotSV_Annotations_Human_3.4.tar.gz"
        String sample_id
        Int memSizeGB = 8
        Int threadCount = 2
        String docker = "quay.io/biocontainers/annotsv"
        Int diskSizeGB = 5*round(size(input_vcf, "GB") + size(annotsv_db_tar_gz, "GB")) + 30
    }
    
    command <<<
        set -eux -o pipefail

        # annotate SVs
        tar -xzf ~{annotsv_db_tar_gz}
        AnnotSV -annotationsDir . -SvinputFile ~{input_vcf} -outputDir out_AnnotSV | tee ~{sample_id}.AnnotSV.log
    >>>

    output {
        File tsv = "~{sample_id}.somaticSV.annotated.tsv"
        File log = "~{sample_id}.AnnotSV.log"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: docker
        preemptible: 1
    }
}