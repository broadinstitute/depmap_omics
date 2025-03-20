version 1.0


workflow RNAStrandWorkflow {
    input {
        File input_bam
        File refseq = "gs://cds-vep-data/ucsc_ncbi_refseq_curated_hg38_20240303.bed"

        String sample_id
        String docker_image = "us-docker.pkg.dev/depmap-omics/public/strand:24q2"
    }

    call strand_check_task as strand {
        input:
            input_bam=input_bam,
            refseq=refseq,
            sample_id=sample_id,
            docker_image=docker_image,
    }

    output {
        Float percentage_1pp1mm = strand.line_1pp1mm
        Float percentage_1pm1mp = strand.line_1pm1mp
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
        Int preemptible = 2
        Int cpu = 1
        Int mem_gb = 8
        Int additional_disk_gb = 0
    }

    Int disk_space = ceil(size(input_bam, "GiB")) + 10 + additional_disk_gb

    command <<<
        infer_experiment.py -r ~{refseq} -i ~{input_bam} > results.txt
        
        sed -n '5p' results.txt | grep -oE "[^:]+$" > 1pp1mm
        sed -n '6p' results.txt | grep -oE "[^:]+$" > 1pm1mp
    >>>

    runtime {
        disks: "local-disk ~{disk_space} SSD"
        memory: "~{mem_gb} GB"
        cpu: cpu
        preemptible: preemptible
        docker: docker_image
    }

    output {
        Float line_1pp1mm = read_float("1pp1mm")
        Float line_1pm1mp = read_float("1pm1mp")
    }

    meta {
        allowNestedInputs: true
    }
}
