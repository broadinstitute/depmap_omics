# Given a set of samples, combine segment files into a single file
# more information available at https://open-cravat.readthedocs.io/en/latest/2.-Command-line-usage.html
workflow opencravat {
    call run_opencravat
}

task run_opencravat {
    String sample_id
    Array[File] vcfs
    String format = "vcf"
    String annotators_to_use = ""
    Int stripfolder = 0 
    String genome = "hg38"
    String module_option = "vcfreporter.type=separate"
    
    Int memory = 16
    Int boot_disk_size = 20
    Int disk_space
    Int num_threads = 1
    Int num_preempt = 5
    String docker = "karchinlab/opencravat"

    command {
      set -euo pipefail
        
      # regular version
      # ---------------
      oc module install-base
      oc module install -y ${annotators_to_use} vcfreporter
      
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
      oc run ${sep=" " vcfs} \
        -l ${genome} \
        -t ${format} \
        --mp ${num_threads} \
        --separatesample \
        ${"--module-option " + module_options}

      gzip ${sep(" ", suffix(${format}, vcfs))}
    }

    output {
        Array[File] oc_error_files="*.err"
        Array[File] oc_log_files="*.log"
        Array[File] oc_sql_files="*.sqlite"
        Array[File] oc_main_files="*.${format}.gz"
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
