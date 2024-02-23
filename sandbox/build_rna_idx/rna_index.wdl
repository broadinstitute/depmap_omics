version 1.0


workflow RNAIndex {
    input {
        File fasta
        File gtf
        String sample_id
    }

    call star_build_index {
        input:
            File fasta,
            File gtf,
            String sample_id
    }

    call rsem_build_index {
        input:
            File fasta,
            File gtf,
            String sample_id
    }

    output {
        File star_tar = star_build_index.output_tar
        File rsem_tar = rsem_build_index.output_tar
    }
}


# Standard interface to run vcf to maf
task star_build_index {
    input {
        File fasta
        File gtf
        String sample_id

        String docker_image="gcr.io/broad-cga-francois-gtex/gtex_rnaseq:V10"
        String assembly="GRCh38"
        Int preemptible=2
        Int boot_disk_size=60
        Int disk_space=60
        Int cpu = 4
        Int mem = 80
    }

    command {
        STAR --runThreadN ~{cpu} --runMode genomeGenerate --genomeDir ~{sample_id} --genomeFastaFiles ~{fasta} --sjdbGTFfile ~{gtf} --sjdbOverhang 100
        tar cvf ~{sample_id}
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
        File output_tar = "~{sample_id}"
    }
}

# Standard interface to run vcf to maf
task rsem_build_index {
    input {
        File fasta
        File gtf
        String sample_id

        String docker_image="gcr.io/broad-cga-francois-gtex/gtex_rnaseq:V10"
        String assembly="GRCh38"
        Int preemptible=2
        Int boot_disk_size=60
        Int disk_space=60
        Int cpu = 4
        Int mem = 80
    }

    command {
        rsem-prepare-reference --gtf ~{gtf} ~{fasta} ~{sample_id}
        tar cvf ~{sample_id}
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
        File output_tar = "~{sample_id}"
    }
}
