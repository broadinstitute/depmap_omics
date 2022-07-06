version 1.0

# Given a set of samples, combine segment files into a single file
# more information available at https://open-cravat.readthedocs.io/en/latest/2.-Command-line-usage.html
workflow run_opencravat {
    input {
        File vcf
    }
    call opencravat {
        input:
            vcf=vcf
    }
    output {
        File oc_error_file=opencravat.oc_error_file
        File oc_log_file=opencravat.oc_log_file
        #File oc_sql_file=opencravat.oc_sql_file
        File oc_main_file=opencravat.oc_main_file
    }
}

task opencravat {
    input {
        File vcf
        String format = "vcf"
        Array[String] annotators_to_use = []
        #Int stripfolder = 0 
        String genome = "hg38"
        String modules_options = "vcfreporter.type=separate"
        # see https://github.com/rkimoakbioinformatics/oak-cravat-modules/tree/master/annotators/oncokb
        File? oncokb_api_key

        Int memory = 16
        Int boot_disk_size = 100
        Int retries=1
        Int disk_space = 20
        Int num_threads = 4
        Int num_preempt = 2
        String docker = "karchinlab/opencravat"
    }
    
    command {
        set -euo pipefail
            
        pip install open-cravat --upgrade

        oc module install-base
        oc module install -y vcfreporter hg19 ~{sep=" " annotators_to_use}

        ${if defined(oncokb_api_key) then "mv "+oncokb_api_key+" /usr/local/lib/python3.6/site-packages/cravat/modules/annotators/oncokb_dm/data/token.txt" else ""}
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
            -a ~{sep=" " annotators_to_use}
    
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
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
        maxRetries: "${retries}"
    }

    meta {
        author: "Jeremie Kalfon"
    }
}