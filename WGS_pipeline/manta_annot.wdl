version 1.0

workflow Manta_Annotator {
    # more info at https://github.com/acranej/MantaSVAnnotator
    input {
        Int cores
        File sv
    }

    call run_Manta_Annotator {
        input:
            sv=sv,
            cores=cores
    }
}

task run_Manta_Annotator {
    input {
        File sv
        
        Int mem_size = 2
        Int disk_size = 20
        Int cores = 4
        String docker_image="jkobject/manta_annot"
        Int preemptible_tries = 3
    }

    command {
        git clone https://github.com/acranej/MantaSVAnnotator.git
        Rscript MantaSVAnnotator/MANTA_vcf2bedpe.R \
        -i ${sv} \
        -o ./

        Rscript MantaSVAnnotator/Manta_SV_Annotator_2.R \
        -i out/${sv}.bedpe \
        -r MantaSVAnnotator/gencode_hg38_annotations_table.txt \
        -g MantaSVAnnotator/gnomad_germline_hg38all.txt \
        -o ./out/ \
        -c {cores}
    }
    # drops germline, <1kb, indels, IMPRECISE, alt/mito/Y chromosomes, failed QC
    runtime {
        preemptible: preemptible_tries
        docker: docker_image
        memory: mem_size
        cpu: cores
        disks: "local-disk " + disk_size + " HDD"
    }

    output {
        File somatic_annotated_sv = "out/${sv}.bedpe.sv_annotation"
        File filtered_annotated_sv = "out/${sv}.bedpe.sv_annotation"
        File indel_vcf = "out/.indels"
        File dropped= ""
    }
}