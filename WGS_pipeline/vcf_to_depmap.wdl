version 1.0

workflow run_vcf_to_depmap {
    input {
        String sample_id
        File input_vcf
    }

    call vcf_to_depmap {
        input:
            input_vcf=input_vcf,
            sample_id=sample_id
    }

    output {
        Array[File] full_file=vcf_to_depmap.full_file
        File depmap_maf=vcf_to_depmap.depmap_maf
    }
}

task vcf_to_depmap {
    input {
        File input_vcf
        String sample_id
        Boolean use_multi=false
        Boolean onco_kb=false
        Array[String] force_keep = ['clinvar_vcf_ssr', 'clinvar_vcf_clndisdbincl', 'clinvar_vcf_clnsigincl', 'clinvar_vcf_clndnincl', 'oc_oncokb__all', 'oc_oncokb__highestdiagnosticimplicationlevel', 'pid']
        Int n_rows=500000

        String docker_image="python"
        Int preemptible=3
        Int boot_disk_size=10
        Int cpu = 1
        String mem = "16 GB"
    }

    command {
        pip install broad-genepy
        git clone https://github.com/broadinstitute/depmap_omics.git
        cd depmap_omics && git checkout dev && pip install -e . && cd ..
        python depmap_omics/WGS_pipeline/vcf_to_depmap.py \
              ~{input_vcf} \
              ~{sample_id} \
              ~{n_rows} \
              ~{use_multi} \
              ~{onco_kb} \
              '~{sep="," force_keep}'
    }

    runtime {
        disks: "local-disk 10 HDD"
        memory: mem
        cpu: cpu
        preemptible: preemptible
        bootDiskSizeGb: boot_disk_size
        docker: docker_image
    }

    output {
        Array[File] full_file = glob("~{sample_id}-maf-full.parquet/*")
        File depmap_maf = "~{sample_id}-maf-coding_somatic-subset.csv.gz"
        
    }
}