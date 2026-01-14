version 1.0

workflow run_mask_variants {
    input {
        File vcf
        String sample_id
    }

    call bcftools_annotate {
        input:
            vcf=vcf,
            sample_id=sample_id
    }

    output {
        File mask_annotated = bcftools_annotate.mask_annotated
    }
}

task bcftools_annotate {
    input {
        File vcf
        String sample_id
        File segdup_bed="gs://ccleparams/segDup_majorAllele_withAltContigs_98pcFracMatch_merged_forBcftools.bed.gz"
        File segdup_bed_index="gs://ccleparams/segDup_majorAllele_withAltContigs_98pcFracMatch_merged_forBcftools.bed.gz.tbi"
        File repeatmasker_bed="gs://ccleparams/repeatMasker_max10_noAlt_merged_forBcftools.bed.gz"
        File repeatmasker_bed_index="gs://ccleparams/repeatMasker_max10_noAlt_merged_forBcftools.bed.gz.tbi"
    
        Int memory = 4
        Int boot_disk_size = 10
        Int num_threads = 1
        Int num_preempt = 5
        Int disk_space = 40
        String docker = "dceoy/bcftools"
    }

    command {

        bcftools annotate \
            -a ~{segdup_bed} \
            -c CHROM,FROM,TO,SEGDUP \
            -h <(echo '##INFO=<ID=SEGDUP,Number=1,Type=String,Description="If variant is in a segmental duplication region">') \
            ~{vcf} > ~{sample_id}.segdup.vcf
        
        bcftools annotate \
            -a ~{repeatmasker_bed} \
            -c CHROM,FROM,TO,RM \
            -h <(echo '##INFO=<ID=RM,Number=1,Type=String,Description="If variant is in a Repeat Masker region">') \
            ~{sample_id}.segdup.vcf > ~{sample_id}.segdup.rm.vcf

        bgzip ~{sample_id}.segdup.rm.vcf

    }

    runtime {
        docker: docker
        bootDiskSizeGb: boot_disk_size
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} SSD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    output {
        File mask_annotated = "~{sample_id}.segdup.rm.vcf.gz"
    }
}
