version 1.0

# Given a set of samples, combine segment files into a single file
# more information available at https://open-cravat.readthedocs.io/en/latest/2.-Command-line-usage.html
# more information about the modules used: 
# https://docs.google.com/document/d/1CdeeprU3oNE9j9-rwjRX_CommeRqwW78Lu2HEHOLdK0

workflow run_opencravat {
    input {
        File vcf
        Int boot_disk_size = 100
        Int disk_space = 100
        Int memory = 32
    }
    call opencravat {
        input:
            vcf=vcf,
            boot_disk_size=boot_disk_size,
            disk_space=disk_space,
            memory=memory
    }
    output {
        File oc_error_file=opencravat.oc_error_file
        File oc_log_file=opencravat.oc_log_file
        #File oc_sql_file=opencravat.oc_sql_file
        File oc_main_file=opencravat.oc_main_file
    }
}

# only the new version of opencravat actually works and it is not in this docker
task opencravat {
    input {
        File vcf
        File oc_modules = "gs://ccleparams/oc_modules_trimmed.tar.gz"# a tar ball of the entie oc module folder (must start with the module folder in the path)
        String format = "vcf"
        File cosmic_annotation = "gs://cds-cosmic/cosmic_cmc_20230509_tier123.csv"
        File oncokb_annotation = "gs://cds-oncokb-data/OncoKB_Annotated_Final_2024-03-19_08-12-54.csv"
        Array[String] annotators_to_use = ["brca1_func_assay", "provean", "revel", "spliceai", "gtex", "pharmgkb", "dida", "gwas_catalog", "ccre_screen", "alfa"]
        #Int stripfolder = 0 
        String genome = "hg38"
        String modules_options = "vcfreporter.type=separate"

        Int memory = 32
        Int boot_disk_size = 100
        Int retries=1
        Int disk_space = 100
        Int num_threads = 4
        Int num_preempt = 2
        String docker = "karchinlab/opencravat:2.2.6"
    }
    String oc_install = "oc module install-base && oc module install -y vcfreporter hg19 cscape_coding civic brca1_func_assay sift provean dann_coding revel spliceai gtex funseq2 pharmgkb dida gwas_catalog mavedb alfa ccre_screen"
    
    command {
        set -euo pipefail
        
        # only the new version of opencravat actually works and it is not in this docker
        pip install open-cravat --upgrade
        
        # only the new version of opencravat actually works and it is not in this docker
        ${if defined(oc_modules) then "cd /usr/local/lib/python3.6/site-packages/cravat/ && tar -xzvf "+oc_modules+" && cd -" else oc_install}
        
        oc new annotator hess_drivers
        rm -r /usr/local/lib/python3.6/site-packages/cravat/modules/annotators/hess_drivers

        oc new annotator cosmic_sig

        oc new annotator oncokb

        git clone https://github.com/broadinstitute/depmap_omics.git
        cd depmap_omics && git checkout add-oncokb-to-oc && git pull && cd ..
        cp -r depmap_omics/WGS_pipeline/hess_drivers /usr/local/lib/python3.6/site-packages/cravat/modules/annotators/
        
        mkdir depmap_omics/WGS_pipeline/cosmic_sig/data && cp ${cosmic_annotation} depmap_omics/WGS_pipeline/cosmic_sig/data/cosmic.csv
        cp -r depmap_omics/WGS_pipeline/cosmic_sig /usr/local/lib/python3.6/site-packages/cravat/modules/annotators/ 

        mkdir depmap_omics/WGS_pipeline/oncokb/data && cp ${oncokb_annotation} depmap_omics/WGS_pipeline/oncokb/data/oncokb.csv
        cp -r depmap_omics/WGS_pipeline/oncokb /usr/local/lib/python3.6/site-packages/cravat/modules/annotators/ 

        pip install bgzip pytabix scipy
        
        echo """
import re
import sys
import gzip
import shutil
print(sys.argv)

done=False
with open(sys.argv[1],'rb') as f:
    with gzip.open(sys.argv[2],'wb') as fout:
        for i, line in enumerate(f):
            original_string = line.decode('utf-8')
            if original_string[0] == '#':
                if original_string.startswith('##INFO=<ID=OC_provean__prediction'):
                    original_string = original_string.replace('\"D(amaging)\"', 'D(amaging)').replace('\"N(eutral)\"', 'N(eutral)')
                    done = True
                fout.write(original_string.encode())
            else:
                done = True
            if done:
                break
        shutil.copyfileobj(f, fout)
""" > fix_name.py

        oc run ${vcf} \
            -l ${genome} \
            -t ${format} \
            --mp ${num_threads} \
            ${"--module-option "+modules_options} \
            -d out \
            -a hess_drivers cosmic_sig oncokb ~{sep=" " annotators_to_use}
    
        python fix_name.py out/${basename(vcf)}.${format} out/${basename(vcf, '.vcf.gz')}.${format}.gz
    }

    output {
        File oc_error_file="out/${basename(vcf)}.err"
        File oc_log_file="out/${basename(vcf)}.log"
        #File oc_sql_file="out/${basename(vcf)}.sqlite"
        File oc_main_file="out/${basename(vcf, '.vcf.gz')}.${format}.gz"
    }

    runtime {
        docker: docker
        bootDiskSizeGb: "${boot_disk_size}"
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
