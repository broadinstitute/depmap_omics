# This workflow runs STAR Fusion after STAR has already been run
task StarFusion {
    File star_chimeric_junctions
    Array[File] ctat_genome_lib_build_dir_files
    Array[File] ref_genome_fa_star_idx_files

    String prefix

    Int memory
    Int disk_space
    Int num_threads
    Int num_preempt
    String docker_source

    command {
        # We need to recreate the directory structure for STAR-Fusion
        mkdir ctat_genome_lib_build_dir
        for path in "${sep='" "' ctat_genome_lib_build_dir_files}"; do
            mv $path ctat_genome_lib_build_dir
        done
        mkdir ctat_genome_lib_build_dir/ref_genome.fa.star.idx
        for path in "${sep='" "' ref_genome_fa_star_idx_files}"; do
            mv $path ctat_genome_lib_build_dir/ref_genome.fa.star.idx
        done

        ls ctat_genome_lib_build_dir | xargs echo

        # Unzip the junctions file, otherwise STAR-Fusion complains
        gunzip -k -c ${star_chimeric_junctions} > chimeric_junctions

        # Run STAR fusion
        /usr/local/src/STAR-Fusion/STAR-Fusion --genome_lib_dir ctat_genome_lib_build_dir \
             -J chimeric_junctions \
             --no_annotation_filter \
             --min_FFPM 0 \
             --output_dir star_fusion_outdir
            
        cp star_fusion_outdir/star-fusion.fusion_predictions.abridged.tsv "${prefix}.star-fusion.fusion_predictions.abridged.tsv"
        cp star_fusion_outdir/star-fusion.fusion_predictions.tsv "${prefix}.star-fusion.fusion_predictions.tsv"
    }

    output {
        File fusion_predictions = "${prefix}.star-fusion.fusion_predictions.tsv"
        File fusion_predictions_abridged = "${prefix}.star-fusion.fusion_predictions.abridged.tsv"
    }

    runtime {
        docker: "${docker_source}"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }
}

#FusionFilter is what was used to generate the reference used in this pipeline.
workflow STAR_Fusion_standalone {
    call StarFusion {

    }
}
