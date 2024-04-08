version 1.0


workflow RNAStrandWorkflow {
    input {
        File input_bam
        File refseq = "gs://cds-vep-data/ucsc_ncbi_refseq_curated_hg38_20240303.bed"

        String sample_id
        String docker_image = "qianqin/strand:latest"

        Int boot_disk_size=60
        Int disk_space=60
        Int cpu = 10
        Int mem = 80
    }

    call strand_check_task as strand {
        input:
            input_bam=input_bam,
            refseq=refseq,
            sample_id=sample_id,

            docker_image=docker_image,
            boot_disk_size=boot_disk_size,
            disk_space=disk_space,
            cpu = cpu,
            mem = mem
    }

    output {
        File rna_strand_output=strand.strand_info
    }
}


# Standard interface to run vcf to maf
task strand_check_task {
    input {
        File input_bam
        File refseq
        String sample_id

        String docker_image
        Int preemptible=2
        Int boot_disk_size=60
        Int disk_space=60
        Int cpu = 10
        Int mem = 80
    }

    command {
        infer_experiment.py -r ~{refseq} -i ~{input_bam}
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
        File strand_info = stdout()
    }
}
