version 1.0

task rnaseqc2 {

    input {
        File bam_file
        File genes_gtf
        String sample_id
        String? strandedness 
        File? intervals_bed
        String? flags

        Int memory = 8
        Int disk_space = 75
        Int num_threads = 4
        Int num_preempt = 1
    }

    command {
        set -euo pipefail
        echo $(date +"[%b %d %H:%M:%S] Running RNA-SeQC 2")
        touch ${sample_id}.fragmentSizes.txt
        rnaseqc ${genes_gtf} ${bam_file} . -s ${sample_id} ${"--bed " + intervals_bed} ${"--stranded " + strandedness} -vv ${flags}
        echo "  * compressing outputs"
        gzip *.gct
        echo $(date +"[%b %d %H:%M:%S] done")
    }

    output {
        File gene_tpm = "${sample_id}.gene_tpm.gct.gz"
        File gene_counts = "${sample_id}.gene_reads.gct.gz"
        File exon_counts = "${sample_id}.exon_reads.gct.gz"
        File metrics = "${sample_id}.metrics.tsv"
        File insertsize_distr = "${sample_id}.fragmentSizes.txt"
    }

    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/gtex_rnaseq:V10"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    meta {
        author: "Francois Aguet"
    }
}


workflow rnaseqc2_workflow {
    input {
        File bam_file
        File genes_gtf
        String sample_id
        String? strandedness 
        File? intervals_bed
        String? flags
    }

    call rnaseqc2 {
        input:
            bam_file = bam_file,
            genes_gtf = genes_gtf,
            sample_id = sample_id,
            strandedness = strandedness,
            intervals_bed = intervals_bed,
            flags = flags
    }

    output {
        File gene_tpm = rnaseqc2.gene_tpm
        File gene_counts = rnaseqc2.gene_counts
        File exon_counts = rnaseqc2.exon_counts
        File metrics = rnaseqc2.metrics
        File insertsize_distr = rnaseqc2.insertsize_distr       
    }
}
