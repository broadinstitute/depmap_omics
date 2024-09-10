version 1.0


workflow run_vcf_to_depmap {
    input {
        String sample_id
        File input_vcf
        String version=""
        String docker_image="us-docker.pkg.dev/depmap-omics/public/vcf_to_depmap:test"
    }

    call vcf_to_depmap {
        input:
            input_vcf=input_vcf,
            sample_id=sample_id,
            version=version,
            docker_image=docker_image,
    }

    output {
        Array[File] full_file=vcf_to_depmap.full_file
        File depmap_maf=vcf_to_depmap.depmap_maf
    }
}

# transforms a vcf file (if possible annotated with opencravat) into:
# 1. a set of parquet files (representing one full file)
#     with a much cleaner structure / annotation / logic than the original vcf file
#     with additional columns from parsing the existing ones (see improve method)
# 2. a maf file that represent what depmap releases. it is a tiny subset of the original vcf file
#     as we remove non coding and germline variants (with some whitelisting) (see to_depmap_maf() method)
task vcf_to_depmap {
    input {
        File input_vcf
        String sample_id
        Boolean use_multi=false
        Boolean whitelist=true
        Array[String] force_keep=['oc_brca1_func_assay__class', 'oc_brca1_func_assay__score', 'clinvar_vcf_ssr', 'clinvar_vcf_clndisdbincl', 'clinvar_vcf_clnsigincl', 'clinvar_vcf_clndnincl', 'oc_oncokb__all', 'oc_oncokb__highestdiagnosticimplicationlevel', 'cgc_other_germline_mut', 'cgc_other_syndrome/disease', 'clinvar_vcf_clnsigconf', 'center', 'cosmicfusion_fusion_id', 'tumor_barcode', 'source"', 'clinvar_vcf_filter', 'clinvar_vcf_dbvarid', 'normal_barcode', 'gencode_34_ncbibuild', 'oc_oncokb_dm__highestprognosticimplicationlevel', 'strandq', 'oc_oncokb_dm__highestsensitivelevel', 'ocm', 'oc_oncokb_dm__all', 'seqq', 'nlod', 'contq', 'nalod', 'oc_base__note', 'oc_cancer_hotspots__samples', 'oc_oncokb_dm__highestresistancelevel', 'oc_oncokb_dm__tumorsummary', 'oc_oncokb_dm__highestdiagnosticimplicationlevel', 'oc_hess_drivers__signature', 'oc_hess_drivers__is_driver']
        Int n_rows=100000
        String version=""

        String docker_image="depmapomics:test"
        Int preemptible=3
        Int boot_disk_size=10
        Int disk_space=40
        Int cpu = 4
        Int mem = 32
    }

    command {
        python -u /install/depmapomics/tasks/vcf_to_depmap.py \
              ~{input_vcf} \
              ~{sample_id} \
              --n_rows ~{n_rows} \
              --use_multi ~{use_multi} \
              --force_keep '~{sep="," force_keep}' \
              --whitelist ~{whitelist} \
              --version ~{version} \
    }

    runtime {
        disks: "local-disk ~{disk_space} HDD"
        memory: "~{mem} GB"
        cpu: cpu
        preemptible: preemptible
        bootDiskSizeGb: boot_disk_size
        docker: docker_image
    }

    output {
        Array[File] full_file = glob("~{sample_id}-maf-full.parquet/*.parquet")
        File depmap_maf = "~{sample_id}-maf-coding_somatic-subset.csv.gz"
    }
}
