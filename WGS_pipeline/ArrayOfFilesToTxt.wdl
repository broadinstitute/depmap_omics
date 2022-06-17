version 1.0
workflow ArrayToTxt_workflow {
    input {
        Array[String] array_of_files
        String list_name
    }

	call CreateTxt {
        input:
            array_of_files=array_of_files,
            list_name=list_name
	}

    output {
        File file_list_name = CreateTxt.file_list_name
    }
}

# Create the txt file from the list of files
task CreateTxt {
    input {
        # Command parameters
        Array[String] array_of_files
        String list_name
    }
 
    command {
        mv ${write_lines(array_of_files)} ${list_name}.txt
    }
    output {
        File file_list_name = "${list_name}.txt"
    }
    runtime {
        docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1500064817"
        preemptible: 3
    }
}
	