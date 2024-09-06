version 1.0

task StarFusion {
    
    input {
        File left_fastq
        File right_fastq
        Array[File] ctat_genome_lib_build_dir_files
        Array[File] ref_genome_fa_star_idx_files

        String prefix

        Int memory=64
        Int disk_space=500
        Int num_threads=8
        Int num_preempt=2
        String docker="trinityctat/starfusion:1.7.0"
    }

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

        # Run STAR fusion
        /usr/local/src/STAR-Fusion/STAR-Fusion --genome_lib_dir ctat_genome_lib_build_dir \
             --left_fq ${left_fastq} \
             --right_fq ${right_fastq} \
             --min_FFPM 0 \
             --no_annotation_filter \
             --output_dir star_fusion_outdir
            
        cp star_fusion_outdir/star-fusion.fusion_predictions.abridged.tsv "${prefix}.star-fusion.fusion_predictions.abridged.tsv"
        cp star_fusion_outdir/star-fusion.fusion_predictions.tsv "${prefix}.star-fusion.fusion_predictions.tsv"
    }

    output {
        File fusion_predictions = "${prefix}.star-fusion.fusion_predictions.tsv"
        File fusion_predictions_abridged = "${prefix}.star-fusion.fusion_predictions.abridged.tsv"
    }

    runtime {
        docker: "${docker}"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} SSD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }
}

#FusionFilter is what was used to generate the reference used in this pipeline.
workflow trinity_cleaned {
    
    input {
        String prefix
        File fastq1
        File fastq2
        Array[File] ctat_genome_lib_build_dir_files
        Array[File] ref_genome_fa_star_idx_files
    }
    call StarFusion {
        input: 
            prefix = prefix,
            left_fastq=fastq1,
            right_fastq=fastq2,
            ctat_genome_lib_build_dir_files=ctat_genome_lib_build_dir_files,
            ref_genome_fa_star_idx_files=ref_genome_fa_star_idx_files
    }
    output {
        File fusion_predictions = StarFusion.fusion_predictions
        File fusion_predictions_abridged = StarFusion.fusion_predictions_abridged
    }
}
