# Given a set of samples, combine segment files into a single file
# more information available at https://open-cravat.readthedocs.io/en/latest/2.-Command-line-usage.html
workflow opencravat {
    call run_opencravat
}

task run_opencravat {
    String sample_id
    File vcf
    String format
    String annotators_to_use
    
    Int memory
    Int boot_disk_size
    Int disk_space
    Int num_threads
    Int num_preempt

    command {
      set -euo pipefail
        
      oc module install-base
      oc module install -y ${annotators_to_use}
      oc run ${vcf} -l hg38 -t ${format} 

      mv ${vcf}.err ${sample_id}.variant_annotations.err
      mv ${vcf}.log ${sample_id}.variant_annotations.log
      mv ${vcf}.sqlite ${sample_id}.variant_annotations.sqlite
      mv ${vcf}.tsv ${sample_id}.variant_annotations.tsv

      gzip ${sample_id}.*
    }

    output {
        File oc_error_file="${sample_id}.variant_annotations.err.gz"
        File oc_log_file="${sample_id}.variant_annotations.log.gz"
        File oc_sqlite_file="${sample_id}.variant_annotations.sqlite.gz"
        File oc_tsv_file="${sample_id}.variant_annotations.tsv.gz"
    }

    runtime {
        docker: "karchinlab/opencravat"
        bootDiskSizeGb: "${boot_disk_size}"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    meta {
        author: "David Wu"
    }
}
