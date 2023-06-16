version 1.0


workflow run_echtvar_cosmic {
    input {
        String sample_id
        File input_vcf
    }

    call echtvar_cosmic {
        input:
            input_vcf=input_vcf,
            sample_id=sample_id
    }

    output {
        File cosmic_annotated_vcf=echtvar_cosmic.annotated_vcf
    }
}

task echtvar_cosmic {
    input {
        File input_vcf
        String sample_id

        String docker_image="depmapomics:test"
        Int preemptible=3
        Int boot_disk_size=10
        Int disk_space=40
        Int cpu = 4
        Int mem = 32
    }

    command {
        
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
        File annotated_vcf = ""        
    }
}
