version 1.0

workflow run_bedpe_to_depmap {
    input {
        String sample_id
        File input_bedpe
    }

    call bedpe_to_depmap {
        input:
            input_bedpe=input_bedpe,
            sample_id=sample_id
    }

    output {
        File expanded_sv_bedpe=bedpe_to_depmap.expanded_bedpe
        File expanded_filtered_sv_bedpe=bedpe_to_depmap.expanded_filtered_bedpe
    }
}

task bedpe_to_depmap {
    input {
        File input_bedpe
        String sample_id

        Int mem_gb = 8
        Int cpu = 1
        Int preemptible = 3
        Int max_retries = 0
        Int additional_disk_gb = 0
    }

    command {
        python -u /home/bedpe_to_depmap.py \
              ~{input_bedpe} \
              ~{sample_id} 
    }

    runtime {
        docker: "us-docker.pkg.dev/depmap-omics/public/bedpe_to_depmap:production"
        memory: "~{mem_gb} GB"
        disks: "local-disk ~{disk_space} SSD"
        preemptible: preemptible
        maxRetries: max_retries
        cpu: cpu
    }

    output {
        File expanded_bedpe = "~{sample_id}.svs.expanded.bedpe"
        File expanded_filtered_bedpe = "~{sample_id}.svs.expanded.filtered.bedpe"
    }

    meta {
        allowNestedInputs: true
    }
}
