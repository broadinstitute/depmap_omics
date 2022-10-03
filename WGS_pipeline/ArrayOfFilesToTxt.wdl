workflow ArrayToTxt_workflow {
	call CreateTxt {

	}
}

# Create the txt file from the list of files
task CreateTxt {
  # Command parameters
  Array[String] array_of_files
  String list_name
 
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
	