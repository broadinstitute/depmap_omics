# Given a set of samples, combine segment files into a single file
workflow aggregate_set_files_workflow {
    call aggregate_set_files
}

task aggregate_set_files {
    Array[File] sample_files
    String output_file_name
    File aggregate_files_script
    
    Int memory
    Int disk_space
    Int num_preempt

    
    command {
        Rscript ${aggregate_files_script} "${output_file_name}" ${write_lines(sample_files)} 
    }

    output {
        File output_merged_file="${output_file_name}"
    }

    runtime {
        docker: "flyingrobin/cds_shiny"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        preemptible: "${num_preempt}"
    }

    meta {
        author: "Guillaume Kugener"
    }
}