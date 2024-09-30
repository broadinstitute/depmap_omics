task aggregate_csvs {


    Array[File] in_files
    String out_file_name
    File aggregate_csvs_script
    
    Int memory
    Int disk_space
    Int num_preempt

    
    command {
        Rscript ${aggregate_csvs_script} ${write_lines(in_files)} "${out_file_name}"
    }

    output {
        File outfile="${out_file_name}"
    }

    runtime {
        docker: "colganwi/tidyverse_plus"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        preemptible: "${num_preempt}"
    }

    meta {
        author: "William Colgan"
    }
}


workflow Aggregate_CSVs {
    call aggregate_csvs
}
