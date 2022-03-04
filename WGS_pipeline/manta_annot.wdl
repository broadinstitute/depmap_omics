version 1.0

workflow run_manta_annotator {
    # more info at https://github.com/acranej/MantaSVAnnotator
    input {
        File sv
        File gene_annot="gs://ccle_default_params/ensembl_102_gene_exons.filt.tsv"
    }

    call manta_annotator {
        input:
            sv=sv,
            gene_annot=gene_annot
    }
}

task manta_annotator {
    input {
        File sv
        File gene_annot

        Int mem_size = 2
        Int disk_size = 20
        Int cores = 8
        String docker_image="jkobject/manta_annot"
        Int preemptible_tries = 3
    }

    command {
        git clone https://github.com/acranej/MantaSVAnnotator.git
        mkdir out
        
        Rscript MantaSVAnnotator/MANTA_vcf2bedpe.R \
        -i ${sv} \
        -o ./out/

        Rscript MantaSVAnnotator/Manta_SV_Annotator_2.R \
        -i out/${sv}.bedpe \
        -r MantaSVAnnotator/gencode_hg38_annotations_table.txt
        ~{"-e " + gene_annot} \
        -g MantaSVAnnotator/gnomad_germline_hg38all.txt \
        -o ./out/ \
        -c ${cores}
    }
    # drops germline, <1kb, indels, IMPRECISE, alt/mito/Y chromosomes, failed QC
    runtime {
        preemptible: "${preemptible_tries}"
        docker: docker_image
        memory: "{mem_size}GB"
        cpu: "${cores}"
        disks: "local-disk ${disk_size} HDD"
    }

    output {
        File somatic_annotated_sv = "out/${sv}.somatic_only_sv.annotated.bedpe"
        File filtered_annotated_sv = "out/${sv}.sv.annotated.bedpe"
        File dropped= "out/${sv}_removed_calls"
    }
}