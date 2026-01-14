version 1.0

# Given a set of samples, combine segment files into a single file
# more information available at https://open-cravat.readthedocs.io/en/latest/2.-Command-line-usage.html
# more information about the modules used: 
# https://docs.google.com/document/d/1CdeeprU3oNE9j9-rwjRX_CommeRqwW78Lu2HEHOLdK0

workflow run_opencravat {
    input {
        File vcf
    }
    call opencravat {
        input:
            vcf=vcf
    }
    output {
        File vcf_out=opencravat.vcf_out
    }
}

# only the new version of opencravat actually works and it is not in this docker
task opencravat {
    input {
        File vcf
        File? oc_modules # a tar ball of the entie oc module folder (must start with the module folder in the path)
        String format = "vcf"
        File cosmic_annotation = "gs://cds-cosmic/cosmic_cmc_20230509_tier123.csv"
        Array[String] annotators_to_use = []
        #Int stripfolder = 0 
        String genome = "hg38"
        String modules_options = "vcfreporter.type=separate"

        Int memory = 16
        Int boot_disk_size = 100
        Int retries=1
        Int disk_space = 20
        Int num_threads = 4
        Int num_preempt = 2
        String docker = "karchinlab/opencravat:2.2.6"
    }
    String oc_install = "oc module install-base && oc module install -y vcfreporter hg19 cscape_coding civic brca1_func_assay sift provean dann_coding revel spliceai gtex funseq2 pharmgkb dida gwas_catalog mavedb alfa ccre_screen"
    
    command {
        set -euo pipefail
        
        # only the new version of opencravat actually works and it is not in this docker
        pip install open-cravat --upgrade

        mkdir /usr/local/lib/python3.6/site-packages/cravat/modules/
        mkdir /usr/local/lib/python3.6/site-packages/cravat/modules/annotators/
        
        oc new annotator hess_drivers

        oc new annotator cosmic_sig

        oc module install -y hg38 cravat-converter vcf-converter tagsampler varmeta vcfinfo casecontrol vcfreporter

        git clone https://github.com/broadinstitute/depmap_omics.git
        cd depmap_omics && git checkout add-cosmic-to-oc && git pull && cd ..
        cp -r depmap_omics/WGS_pipeline/hess_drivers /usr/local/lib/python3.6/site-packages/cravat/modules/annotators/
        
        mkdir depmap_omics/WGS_pipeline/cosmic_sig/data && cp ${cosmic_annotation} depmap_omics/WGS_pipeline/cosmic_sig/data/cosmic.csv
        cp -r depmap_omics/WGS_pipeline/cosmic_sig /usr/local/lib/python3.6/site-packages/cravat/modules/annotators/

        pip install bgzip pytabix scipy

        oc run ${vcf} \
            -l hg38 \
            -t ${format} \
            --mp ${num_threads} \
            ${"--module-option "+modules_options} \
            -d out \
            -a cosmic_sig
    }

    output {
        File vcf_out = "out/${basename(vcf)}.${format}"
    }

    runtime {
        docker: docker
        bootDiskSizeGb: boot_disk_size
        memory: "${memory} GB"
        disks: "local-disk ${disk_space} SSD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
        maxRetries: "${retries}"
    }

    meta {
        author: "Jeremie Kalfon"
    }
}
