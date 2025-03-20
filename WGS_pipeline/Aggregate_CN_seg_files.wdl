# Given a set of samples, combine segment files into a single file
workflow aggregate_CN_segments_wrkflw {
    call aggregate_CN_segments
}

task aggregate_CN_segments {
    Array[File] sample_seg_files
    String sample_set_id
    File aggregate_seg_files_script
    
    Int memory
    Int disk_space
    Int num_preempt

    
    command {
        Rscript ${aggregate_seg_files_script} "${sample_set_id}.called.seg" ${write_lines(sample_seg_files)} 
    }

    output {
        File combined_cn_file="${sample_set_id}.called.seg"
    }

    runtime {
        docker: "flyingrobin/cds_shiny"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} SSD"
        preemptible: "${num_preempt}"
    }

    meta {
        author: "Guillaume Kugener"
    }
}
