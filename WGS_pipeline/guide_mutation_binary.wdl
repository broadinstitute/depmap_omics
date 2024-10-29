version 1.0

workflow run_guide_mutation {
    input {
        String sample_id
        File vcf
        String docker="us-docker.pkg.dev/depmap-omics/public/depmapomics:bcftools"
    }

    call guide_mutation {
        input:
            vcf=vcf,
            sample_id=sample_id,
            docker=docker,
    }

    output {
        File avana_binary_mut=guide_mutation.avana_binary_mut
        File humagne_binary_mut=guide_mutation.humagne_binary_mut
        File ky_binary_mut=guide_mutation.ky_binary_mut
        File brunello_binary_mut=guide_mutation.brunello_binary_mut
        File tkov3_binary_mut=guide_mutation.tkov3_binary_mut
    }
}

task guide_mutation {
    input {
        File vcf
        String sample_id
        String docker

        File avana_bed = "gs://ccleparams/avana_guides.bed"
        File humagne_bed = "gs://ccleparams/humagne_guides.bed"
        File ky_bed = "gs://ccleparams/ky_score_guides.bed"
        File brunello_bed = "gs://ccleparams/brunello_guides.bed"
        File TKOv3_bed = "gs://ccleparams/tkov3_guides.bed"
        String bcftools_format = '"%CHROM\\t%POS\\t%END\\t%ALT{0}\n"'
    
        Int memory = 4
        Int boot_disk_size = 10
        Int num_threads = 1
        Int num_preempt = 5
        Int disk_space = 40
    }

    command <<<
        set -euo pipefail

        bcftools index -t ~{vcf}

        bcftools query \
            --exclude "FILTER!='PASS'&GT!='mis'&GT!~'\.'" \
            --regions-file ~{avana_bed} \
            --format ~{bcftools_format} \
            ~{vcf} > avana_~{sample_id}.bed

        bcftools query\
            --exclude "FILTER!='PASS'&GT!='mis'&GT!~'\.'"\
            --regions-file ~{humagne_bed} \
            --format ~{bcftools_format} \
            ~{vcf} > humagne_~{sample_id}.bed

        bcftools query\
            --exclude "FILTER!='PASS'&GT!='mis'&GT!~'\.'"\
            --regions-file ~{ky_bed} \
            --format ~{bcftools_format} \
            ~{vcf} > ky_~{sample_id}.bed

        bcftools query\
            --exclude "FILTER!='PASS'&GT!='mis'&GT!~'\.'"\
            --regions-file ~{brunello_bed} \
            --format ~{bcftools_format} \
            ~{vcf} > brunello_~{sample_id}.bed

        bcftools query\
            --exclude "FILTER!='PASS'&GT!='mis'&GT!~'\.'"\
            --regions-file ~{TKOv3_bed} \
            --format ~{bcftools_format} \
            ~{vcf} > tkov3_~{sample_id}.bed

        touch "~{sample_id}_avana_mut_binary.csv" "~{sample_id}_humagne_mut_binary.csv" "~{sample_id}_ky_mut_binary.csv" "~{sample_id}_brunello_mut_binary.csv" "~{sample_id}_tkov3_mut_binary.csv"

        python -u /install/depmapomics/tasks/map_to_guides.py \
              --sample_id ~{sample_id} \
              --bed_filenames 'avana_~{sample_id}.bed,humagne_~{sample_id}.bed,ky_~{sample_id}.bed,brunello_~{sample_id}.bed,tkov3_~{sample_id}.bed'\
              --libraries 'avana,humagne,ky,brunello,tkov3' \
              --guides '~{avana_bed},~{humagne_bed},~{ky_bed},~{brunello_bed},~{TKOv3_bed}'
        
    >>>

    output {
        File avana_binary_mut="~{sample_id}_avana_mut_binary.csv"
        File humagne_binary_mut="~{sample_id}_humagne_mut_binary.csv"
        File ky_binary_mut="~{sample_id}_ky_mut_binary.csv"
        File brunello_binary_mut="~{sample_id}_brunello_mut_binary.csv"
        File tkov3_binary_mut="~{sample_id}_tkov3_mut_binary.csv"
    }

    runtime {
        docker: docker
        bootDiskSizeGb: "~{boot_disk_size}"
        memory: "~{memory}GB"
        disks: "local-disk ~{disk_space} HDD"
        cpu: "~{num_threads}"
        preemptible: "~{num_preempt}"
    }

    meta {
        author: "Simone Zhang"
    }
}
