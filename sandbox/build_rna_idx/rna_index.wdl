version 1.0


workflow RNAIndex {
    input {
        File fasta
        File gtf
        String sample_id
    }

    call star_build_index {
        input:
            fasta=fasta,
            gtf=gtf,
            sample_id=sample_id
    }

    call rsem_build_index {
        input:
            fasta=fasta,
            gtf=gtf,
            sample_id=sample_id
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
        Int boot_disk_size=100
        Int disk_space=200
        Int cpu = 12
        Int mem = 80
    }

    command {
        STAR --runThreadN ~{cpu} --runMode genomeGenerate --genomeDir ~{sample_id} --genomeFastaFiles ~{fasta} --sjdbGTFfile ~{gtf} --sjdbOverhang 100
        tar cvfz ~{sample_id}_star.tar ~{sample_id}
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
        File output_tar = "~{sample_id}_star.tar.gz"
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
        Int cpu = 12
        Int mem = 80
    }

    command {
        mkdir -p ~{sample_id}_rsem_references
        rsem-prepare-reference --gtf ~{gtf} ~{fasta} ~{sample_id}
        cp ~{sample_id}.chrlist ~{sample_id}.seq ~{sample_id}.idx.fa ~{sample_id}.ti ~{sample_id}.n2g.idx.fa ~{sample_id}.transcripts.fa ~{sample_id}.grp ~{sample_id}_rsem_references
        tar cvfz ~{sample_id}_rsem.tar.gz ~{sample_id}_rsem_references
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
        File output_tar = "~{sample_id}_rsem.tar.gz"
    }
}
