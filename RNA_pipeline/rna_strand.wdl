version 1.0


workflow RNAStrandWorkflow {
    input {
        File input_bam
        File refseq = "gs://cds-vep-data/ucsc_ncbi_refseq_curated_hg38_20240303.bed"

        String sample_id
        String docker_image = "us-docker.pkg.dev/depmap-omics/public/strand:24q2"

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

    call parse_line as parse_1pp1mm {
        input:
            line = strand.line_1pp1mm
    }

    call parse_line as parse_1pm1mp {
        input:
            line = strand.line_1pm1mp
    }

    output {
        Float percentage_1pp1mm = parse_1pp1mm.percentage
        Float percentage_1pm1mp = parse_1pm1mp.percentage
        # If fraction explained by 1+-,1-+,2++,2-- is higher than 0.7, it is stranded 
        Boolean inferred_stranded = percentage_1pm1mp > 0.7
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
        String line_1pp1mm = read_lines(stdout())[4]
        String line_1pm1mp = read_lines(stdout())[5]
    }
}

task parse_line {
    input {
        String line
    }

    command {
        echo ~{line} | grep -oE "[^:]+$"
    }

    output {
        Float percentage = read_float(stdout())
    }
}