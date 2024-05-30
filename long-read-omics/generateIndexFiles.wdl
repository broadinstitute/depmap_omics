version 1.0

workflow generateIndexFilesWorkflow {
	
	call generateIndexFiles 

	output {
		File rawBamIndex = generateIndexFiles.rawBamIndex		
	}
}

task generateIndexFiles {
	
	input {
		File inputFile
		String sampleName
		Int memory
    		Int disk_size
		String docker
	}

	command {
		#samtools view -bS \
		#	      -o ${sampleName}.merged.sorted.bam \
		#	      ${inputFile}
		
		cp ${inputFile} ./${sampleName}.merged.sorted.bam 
		
		samtools index ${sampleName}.merged.sorted.bam 
	}

	output {
		File rawBamIndex = "${sampleName}.merged.sorted.bam.bai"
	}

	runtime {
		docker: "${docker}"
		memory: "${memory}GB"
		max_retries: 3
		disks: "local-disk " + disk_size + " HDD"
	}
}
