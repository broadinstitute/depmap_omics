# Given a set of samples, combine segment files into a single file
workflow aggregate_CN_segments_wrkflw {
    call aggregate_CN_segments
}

task aggregate_CN_segments {
    Array[File] sample_seg_files
    String sample_set_id
    
    Int memory
    Int disk_space
    Int num_preempt

    
    command {
        git clone https://github.com/broadinstitute/ccle_processing.git
        Rscript ccle_processing/WGS_pipeline/generate_single_seg_file.R "${sample_set_id}.called.seg" ${write_lines(sample_seg_files)} 
    }

    output {
        File combined_cn_file="${sample_set_id}.called.seg"
    }

    runtime {
        docker: "us-docker.pkg.dev/depmap-omics/public/ccle_rnaseq:latest:latest"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        preemptible: "${num_preempt}"
    }

    meta {
        author: "Guillaume Kugener"
    }
}