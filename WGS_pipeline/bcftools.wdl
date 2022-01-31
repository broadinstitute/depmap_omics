# Given a set of samples, combine segment files into a single file
# more information available at https://open-cravat.readthedocs.io/en/latest/2.-Command-line-usage.html
workflow bcftools_fix_ploidy {
    call run_fix_ploidy
}

task run_fix_ploidy {
    File vcf
    String sample_id
 
    Int memory = 2
    Int boot_disk_size = 10
    Int num_threads = 1
    Int num_preempt = 5
    Int disk_space = 10
    String docker = "dceoy/bcftools"

    command {
      set -euo pipefail

      bcftools +setGT ${vcf} -- -t q -i'INFO/DP>8 & AF>0.9' -n c:'m|m' > ${sample_id}.vcf
      bcftools +setGT ${sample_id}.vcf -- -t q -i'AD[*:0]=0 & INFO/DP>8 & GT="0/1/2"' -n c:'1/2' > ${sample_id}.vcf.1
      bcftools +setGT ${sample_id}.vcf.1 -- -t q -i'AD[*:0]=0 & INFO/DP>8 & GT="0/1/2/3"' -n c:'1/2/3' > ${sample_id}.vcf
      bcftools +setGT ${sample_id}.vcf -- -t q -i'AD[*:0]=0 & INFO/DP>8 & GT="0/1/2/3/4"' -n c:'1/2/3/4' > ${sample_id}.vcf.1
      bcftools +setGT ${sample_id}.vcf.1 -- -t q -i'AD[*:0]=0 & INFO/DP>8 & GT="0/1/2/3/4/5"' -n c:'1/2/3/4/5' > ${sample_id}.vcf
      bcftools +setGT ${sample_id}.vcf -- -t q -i'AD[*:0]=0 & INFO/DP>8 & GT="0/1/2/3/4/5/6"' -n c:'1/2/3/4/5/6' > ${sample_id}_fixedploidy.vcf
      
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
