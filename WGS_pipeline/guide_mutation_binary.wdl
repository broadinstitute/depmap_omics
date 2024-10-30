version 1.0

workflow run_guide_mutation {
    input {
        String sample_id
        File vcf
        String docker="us-docker.pkg.dev/depmap-omics/public/depmapomics:bcftools"
    }

    call guide_mutation_in_vcf {
        input:
            vcf=vcf,
            sample_id=sample_id,
            docker=docker,
    }

    call guide_mutation_intersect {
        input:
            sample_id=sample_id,
            avana_in_vcf_bed = guide_mutation_in_vcf.avana_mut,
            humagne_in_vcf_bed = guide_mutation_in_vcf.humagne_mut,
            ky_in_vcf_bed = guide_mutation_in_vcf.ky_mut,
            brunello_in_vcf_bed = guide_mutation_in_vcf.brunello_mut,
            TKOv3_in_vcf_bed = guide_mutation_in_vcf.tkov3_mut,
    }

    output {
        File avana_binary_mut=guide_mutation_intersect.avana_binary_mut
        File humagne_binary_mut=guide_mutation_intersect.humagne_binary_mut
        File ky_binary_mut=guide_mutation_intersect.ky_binary_mut
        File brunello_binary_mut=guide_mutation_intersect.brunello_binary_mut
        File tkov3_binary_mut=guide_mutation_intersect.tkov3_binary_mut
    }
}

task guide_mutation_in_vcf {
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
        
    >>>

    output {
        File avana_mut="avana_~{sample_id}.bed"
        File humagne_mut="humagne_~{sample_id}.bed"
        File ky_mut="ky_~{sample_id}.bed"
        File brunello_mut="brunello_~{sample_id}.bed"
        File tkov3_mut="tkov3_~{sample_id}.bed"
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

task guide_mutation_intersect {
    input {
        String sample_id
        String docker="biocontainers/bedtools:v2.28.0_cv2"

        File avana_in_vcf_bed
        File humagne_in_vcf_bed
        File ky_in_vcf_bed
        File brunello_in_vcf_bed
        File TKOv3_in_vcf_bed

        File avana_bed = "gs://ccleparams/avana_guides.bed"
        File humagne_bed = "gs://ccleparams/humagne_guides.bed"
        File ky_bed = "gs://ccleparams/ky_score_guides.bed"
        File brunello_bed = "gs://ccleparams/brunello_guides.bed"
        File TKOv3_bed = "gs://ccleparams/tkov3_guides.bed"
    
        Int memory = 4
        Int boot_disk_size = 10
        Int num_threads = 1
        Int num_preempt = 5
        Int disk_space = 40
    }

    command <<<
        set -euo pipefail

        bedtools intersect -a ~{avana_bed} -b ~{avana_in_vcf_bed} -c -sorted > ~{sample_id}_avana_mut_binary.bed
        bedtools intersect -a ~{humagne_bed} -b ~{humagne_in_vcf_bed} -c -sorted > ~{sample_id}_humagne_mut_binary.bed
        bedtools intersect -a ~{ky_bed} -b ~{ky_in_vcf_bed} -c -sorted > ~{sample_id}_ky_mut_binary.bed
        bedtools intersect -a ~{brunello_bed} -b ~{brunello_in_vcf_bed} -c -sorted > ~{sample_id}_brunello_mut_binary.bed
        bedtools intersect -a ~{TKOv3_bed} -b ~{TKOv3_in_vcf_bed} -c -sorted > ~{sample_id}_tkov3_mut_binary.bed
        
    >>>

    output {
        File avana_binary_mut="~{sample_id}_avana_mut_binary.bed"
        File humagne_binary_mut="~{sample_id}_humagne_mut_binary.bed"
        File ky_binary_mut="~{sample_id}_ky_mut_binary.bed"
        File brunello_binary_mut="~{sample_id}_brunello_mut_binary.bed"
        File tkov3_binary_mut="~{sample_id}_tkov3_mut_binary.bed"
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