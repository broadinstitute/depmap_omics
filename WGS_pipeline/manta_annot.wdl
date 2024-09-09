version 1.0

workflow run_manta_annotator {
    # more info at https://github.com/acranej/MantaSVAnnotator
    input {
        File sv
        File exon_annot="gs://ccleparams/hg38_ensembl_exonlocations_formatted.txt"

        Int mem_size = 4
        Int disk_size = 30
        Int cores = 8
        String docker_image="jkobject/manta_annot"
        Int preemptible_tries = 3
    }

    call manta_annotator {
        input:
            sv=sv,
            exon_annot=exon_annot,
            mem_size=mem_size,
            disk_size=disk_size,
            cores=cores,
            docker_image=docker_image,
            preemptible_tries=preemptible_tries
    }
    output {
        File somatic_annotated_sv = manta_annotator.somatic_annotated_sv 
        File filtered_annotated_sv = manta_annotator.filtered_annotated_sv 
        File dropped = manta_annotator.dropped 
    }
}

task manta_annotator {
    input {
        File sv
        File exon_annot

        Int mem_size = 4
        Int disk_size = 30
        Int cores = 8
        String docker_image="jkobject/manta_annot"
        Int preemptible_tries = 3
    }

    String newname = sub(sub(sv, "\\.gz$", ""), "gs://","/cromwell_root/")

    command {
        git clone https://github.com/acranej/MantaSVAnnotator.git
        mkdir out

        gunzip ${sv}
        
        Rscript MantaSVAnnotator/MANTA_vcf2bedpe.R \
        -i ${newname} \
        -o ./out/

        Rscript MantaSVAnnotator/Manta_SV_Annotator_2.R \
        -i out/${basename(newname)}.bedpe \
        -r MantaSVAnnotator/hg38_ensembl_genelocatins_formatted.txt \
        ~{"-e " + exon_annot} \
        -g MantaSVAnnotator/gnomad_germline_hg38all.txt \
        -o ./out/ \
        -c ${cores}
    }
    # drops germline, <1kb, indels, IMPRECISE, alt/mito/Y chromosomes, failed QC
    runtime {
        preemptible: "${preemptible_tries}"
        docker: docker_image
        memory: "${mem_size} GB"
        cpu: "${cores}"
        disks: "local-disk ${disk_size} SSD"
    }

    output {
        File somatic_annotated_sv = "out/${basename(newname)}.bedpe_somatic_only_sv_annotated.bedpe"
        File filtered_annotated_sv = "out/${basename(newname)}.bedpe_sv_annotated.bedpe"
        File dropped= "out/${basename(newname)}.removed_calls"
    }
}
