version 1.0

task rsem {

    input {
        File transcriptome_bam
        File rsem_reference
        String prefix

        Int memory = 128
        Int disk_space = 500
        Int num_threads = 32
        Int num_preempt = 1

        Int max_frag_len = 1000
        String estimate_rspd = "true"
        String? is_stranded
        String? paired_end
        String? calc_ci
        Int? ci_memory
    }

    command {
        set -euo pipefail
        mkdir rsem_reference
        tar -xvvf ${rsem_reference} -C rsem_reference --strip-components=1

        git clone https://github.com/broadinstitute/ccle_processing.git

        chmod +x ccle_processing/RNA_pipeline/run_RSEM_david.py

        ccle_processing/RNA_pipeline/run_RSEM_david.py \
            ${"--max_frag_len " + max_frag_len} \
            ${"--estimate_rspd " + estimate_rspd} \
            ${"--is_stranded " + is_stranded} \
            ${"--paired_end " + paired_end} \
            ${"--calc_ci " + calc_ci} \
            ${"--ci_memory " + ci_memory} \
            --threads ${num_threads} \
            rsem_reference ${transcriptome_bam} ${prefix}
        gzip *.results
    }

    output {
        File genes="${prefix}.rsem.genes.results.gz"
        File isoforms="${prefix}.rsem.isoforms.results.gz"
    }

    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/gtex_rnaseq:V10"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    meta {
        author: "David Wu"
    }
}


workflow rsem_workflow {
    input {
        File transcriptome_bam
        File rsem_reference
        String prefix
    }

    call rsem {
        input:
            transcriptome_bam = transcriptome_bam,
            rsem_reference = rsem_reference,
            prefix = prefix
    }

    output {   
        File genes=rsem.genes
        File isoforms=rsem.isoforms
    }
}