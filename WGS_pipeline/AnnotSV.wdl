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
        File annotated_svs_full = annotate_sv_annotsv.full_tsv
        File annotsv_log_full = annotate_sv_annotsv.full_log
        File annotated_svs_split = annotate_sv_annotsv.split_tsv
        File annotsv_log_split = annotate_sv_annotsv.split_log
        File annotated_svs_full_vcf = annotate_sv_annotsv.full_vcf
        File unannotated_svs_full = annotate_sv_annotsv.unannotated
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
        Int diskSizeGB = 5*round(size(input_vcf, "GB") + size(annotsv_db_tar_gz, "GB")) + 50
        File cosmic_cna = "gs://cds-cosmic/CosmicCNA_v97/CosmicCompleteCNA.tsv.gz"
    }
    
    command <<<
        set -eux -o pipefail

        # annotate SVs
        tar -xzf ~{annotsv_db_tar_gz}
        mv ~{cosmic_cna} Annotations_Human/FtIncludedInSV/COSMIC/GRCh38/CosmicCompleteCNA.tsv.gz

        AnnotSV -annotationsDir . -annotationMode full -includeCI 0 -SVminSize 1 -SvinputFile ~{input_vcf} -outputFile ~{sample_id}.AnnotSV.full.tsv -outputDir . -vcf 1 | tee ~{sample_id}.AnnotSV.full.log

        AnnotSV -annotationsDir . -annotationMode split -includeCI 0 -SVminSize 1 -SvinputFile ~{input_vcf} -outputFile ~{sample_id}.AnnotSV.split.tsv -outputDir . -vcf 1 | tee ~{sample_id}.AnnotSV.split.log
    >>>

    output {
        File full_tsv = "~{sample_id}.AnnotSV.full.tsv"
        File full_log = "~{sample_id}.AnnotSV.full.log"
        File split_tsv = "~{sample_id}.AnnotSV.split.tsv"
        File split_log = "~{sample_id}.AnnotSV.split.log"
        File full_vcf = "~{sample_id}.AnnotSV.full.vcf"
        File unannotated = "~{sample_id}.AnnotSV.full.unannotated.tsv"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: docker
        preemptible: 1
    }
}