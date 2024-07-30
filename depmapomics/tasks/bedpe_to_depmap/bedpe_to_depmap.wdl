version 1.0

workflow run_bedpe_to_depmap {
    input {
        String sample_id
        File input_bedpe
        String docker_image="us-docker.pkg.dev/depmap-omics/public/bedpe_to_depmap:test"
    }

    call bedpe_to_depmap {
        input:
            input_bedpe=input_bedpe,
            sample_id=sample_id,
            docker_image=docker_image,
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

        String docker_image="us-docker.pkg.dev/depmap-omics/public/bedpe_to_depmap:test"
        Int preemptible=3
        Int boot_disk_size=10
        Int disk_space=40
        Int cpu = 4
        Int mem = 32
    }

    command {
        python -u /home/bedpe_to_depmap.py \
              ~{input_bedpe} \
              ~{sample_id} 
    }

    runtime {
        disks: "local-disk ~{disk_space} SSD"
        memory: "~{mem} GB"
        cpu: cpu
        preemptible: preemptible
        bootDiskSizeGb: boot_disk_size
        docker: docker_image
    }

    output {
        File expanded_bedpe = "~{sample_id}.svs.expanded.bedpe"
        File expanded_filtered_bedpe = "~{sample_id}.svs.expanded.filtered.bedpe"
    }
}
