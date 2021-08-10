task cram_to_bam {

    File cram_file
    File cram_index
    File reference_fasta
    File reference_fasta_index
    String? region

    Int memory
    Int disk_space
    Int num_threads
    Int num_preempt

    command {
        set -euo pipefail
        echo $(date +"[%b %d %H:%M:%S] Converting CRAM to BAM")
        prefix=$(basename ${cram_file} ".cram")
        bam_file=$prefix".bam"
        samtools view -b -T ${reference_fasta} -o $bam_file ${cram_file} ${region}
        samtools index $bam_file
        echo $(date +"[%b %d %H:%M:%S] done")
    }

    output {
        File bam_file=glob("*.bam")[0]
        File bam_index=glob("*.bam.bai")[0]
    }

    runtime {
        docker: "us-docker.pkg.dev/depmap-omics/public/gtex-rnaseq:V8"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    meta {
        author: "Francois Aguet"
    }
}


workflow cram_to_bam_workflow {
    call cram_to_bam
}
