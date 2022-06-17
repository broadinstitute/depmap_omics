version 1.0

# Given a set of samples, combine segment files into a single file
# more information available at https://open-cravat.readthedocs.io/en/latest/2.-Command-line-usage.html
workflow run_opencravat {
    input {
        String sample_id
        File vcf
        String? format
        String? annotators_to_use
        String? genome
        String? modules_options
        String? docker
        # see https://github.com/rkimoakbioinformatics/oak-cravat-modules/tree/master/annotators/oncokb
        File? oncokb_api_key
    }
    call opencravat {
        input:
            sample_id=sample_id,
            vcf=vcf,
            format=format,
            annotators_to_use=annotators_to_use,
            genome=genome,
            modules_options=modules_options,
            docker=docker,
            oncokb_api_key=oncokb_api_key
    }
    output {
        File oc_error_files=opencravat.oc_error_files
        File oc_log_files=opencravat.oc_log_files
        File oc_sql_files=opencravat.oc_sql_files
        File oc_main_files=opencravat.oc_main_files
    }
}

task opencravat {
    input {
        String sample_id
        File vcf
        String format = "vcf"
        String annotators_to_use = ""
        #Int stripfolder = 0 
        String genome = "hg38"
        String modules_options = "vcfreporter.type=separate"
        File? oncokb_api_key

        Int memory = 16
        Int boot_disk_size = 20
        Int disk_space=20
        Int num_threads = 10
        Int num_preempt = 2
        String docker = "karchinlab/opencravat"
    }
    
    command {
      set -euo pipefail
        
      # regular version
      # ---------------
      pip install open-cravat --upgrade
      oc module install-base
      oc module install -y ${annotators_to_use} vcfreporter hg19
      ${if defined(oncokb_api_key) then "mkdir /usr/local/lib/python3.6/site-packages/cravat/modules/annotators/oncokb/data/" else ""}
      ${if defined(oncokb_api_key) then "mv "+oncokb_api_key+" /usr/local/lib/python3.6/site-packages/cravat/modules/annotators/oncokb/data/token.txt" else ""}
      
      # fast version
      # ------------
      # to make it faster we should use a copy from a bucket using gsutil cp
      # we install everything on a machine, then get path to data using `oc config md`
      # we then cp it on a bucket.
      # now we can add two additional commands here:
      # 1. to copy the content of the bucket here: gsutil cp gs://path/to/modules.gz .
      # 2. to use this location as the md location: oc config md LOCATION
      # gsutil cp [modules] modules.tar
      # tar -tvf modules.tar --strip-components=[stripfolder]
      # oc config md ./modules
      oc run ${vcf} \
        -l ${genome} \
        -t ${format} \
        --mp ${num_threads} \
        ${"--module-option "+modules_options} \
        -d out \
        -a ${annotators_to_use}
      
      pip install bgzip
      
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
      python fix_name.py out/${basename(vcf)}.${format} out/${basename(vcf)}.${format}.gz
    }

    output {
        File oc_error_files="out/${basename(vcf)}.err"
        File oc_log_files="out/${basename(vcf)}.log"
        File oc_sql_files="out/${basename(vcf)}.sqlite"
        File oc_main_files="out/${basename(vcf)}.${format}.gz"
    }

    runtime {
        docker: docker
        bootDiskSizeGb: "${boot_disk_size}"
        memory: "${memory} GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    meta {
        author: "Jeremie Kalfon"
    }
}
