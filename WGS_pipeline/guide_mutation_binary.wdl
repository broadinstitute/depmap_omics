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
}

task guide_mutation {
    input {
        File vcf
        String sample_id

        File avana_bed = "gs://ccleparams/avana_guides.bed"
        File humagne_bed = "gs://ccleparams/humagne_guides.bed"
        File ky_bed = "gs://ccleparams/ky_score_guides.bed"
    
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

        bcftools +setGT ${vcf} -- -t q -i'INFO/DP>8 & AF>0.9' -n c:'m|m' > ${sample_id}.vcf
        bcftools +setGT ${sample_id}.vcf -- -t q -i'AD[*:0]=0 & INFO/DP>8 & GT="0/1/2"' -n c:'1/2' > ${sample_id}.vcf.1
        bcftools +setGT ${sample_id}.vcf.1 -- -t q -i'AD[*:0]=0 & INFO/DP>8 & GT="0/1/2/3"' -n c:'1/2/3' > ${sample_id}.vcf
        bcftools +setGT ${sample_id}.vcf -- -t q -i'AD[*:0]=0 & INFO/DP>8 & GT="0/1/2/3/4"' -n c:'1/2/3/4' > ${sample_id}.vcf.1
        bcftools +setGT ${sample_id}.vcf.1 -- -t q -i'AD[*:0]=0 & INFO/DP>8 & GT="0/1/2/3/4/5"' -n c:'1/2/3/4/5' > ${sample_id}.vcf
        bcftools +setGT ${sample_id}.vcf -- -t q -i'AD[*:0]=0 & INFO/DP>8 & GT="0/1/2/3/4/5/6"' -n c:'1/2/3/4/5/6' > ${sample_id}_fixedploidy.vcf

        bcftools query\
            --exclude \"FILTER!='PASS'&GT!='mis'&GT!~'\.'\"\
            --regions-file ${avana_bed}" \
            --format '%CHROM\\t%POS\\t%END\\t%ALT{0}\n' \
            ${vcf} > avana_${sample_id}.bed

        bcftools query\
            --exclude \"FILTER!='PASS'&GT!='mis'&GT!~'\.'\"\
            --regions-file ${humagne_bed}" \
            --format '%CHROM\\t%POS\\t%END\\t%ALT{0}\n' \
            ${vcf} > humagne_${sample_id}.bed

        bcftools query\
            --exclude \"FILTER!='PASS'&GT!='mis'&GT!~'\.'\"\
            --regions-file ${ky_bed}" \
            --format '%CHROM\\t%POS\\t%END\\t%ALT{0}\n' \
            ${vcf} > ky_${sample_id}.bed
        
        gzip ${sample_id}_fixedploidy.vcf
    }

    output {
        File vcf_fixedploid="${sample_id}_fixedploidy.vcf.gz"
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
        author: "Jeremie Kalfon"
    }
}
