version 1.0

workflow run_guide_mutation {
    input {
        File vcf
        String sample_id
    }

    call guide_mutation {
        input:
            vcf=vcf,
            sample_id=sample_id
    }

    output {
        File avana=guide_mutation.sample_avana_bed
    }
}

task guide_mutation {
    input {
        File vcf
        String sample_id

        File avana_bed = "gs://ccleparams/avana_guides.bed"
        File humagne_bed = "gs://ccleparams/humagne_guides.bed"
        File ky_bed = "gs://ccleparams/ky_score_guides.bed"
        String bcftools_format = "'%CHROM\\t%POS\\t%END\\t%ALT{0}\n'"
    
        Int memory = 4
        Int boot_disk_size = 10
        Int num_threads = 1
        Int num_preempt = 5
        Int disk_space = 40
        String docker = "dceoy/bcftools"
    }

    File vcf_index = vcf + ".tbi"

    command {
        set -euo pipefail

        bcftools query \
            --exclude FILTER!='PASS'&GT!='mis'&GT!~'\.' \
            --regions-file ${avana_bed} \
            --format ${bcftools_format} \
            ${vcf} > avana_${sample_id}.bed

        bcftools query\
            --exclude FILTER!='PASS'&GT!='mis'&GT!~'\.'\
            --regions-file ${humagne_bed} \
            --format ${bcftools_format} \
            ${vcf} > humagne_${sample_id}.bed

        bcftools query\
            --exclude FILTER!='PASS'&GT!='mis'&GT!~'\.'\
            --regions-file ${ky_bed} \
            --format ${bcftools_format} \
            ${vcf} > ky_${sample_id}.bed
        
    }

    output {
        File sample_avana_bed="avana_${sample_id}.bed"
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
        author: "Simone Zhang"
    }
}
