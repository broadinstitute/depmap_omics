version 1.0

workflow helloLongReadMapper {
  call longReadMapper
}

task longReadMapper {
  
  input {

    File RefFasta
    String sampleName
    String readType
    File inputFile

  }

  command {
	
    minimap2 -ax ${readType} \
          ${RefFasta} \
          ${inputFile} \
          > ${sampleName}.sam 

    samtools view \
          -Sb \
          -o ${sampleName}.bam \
          ${sampleName}.sam
	
    samtools sort \
          ${sampleName}.bam \
          -o ${sampleName}.sorted.bam

    samtools index \
          ${sampleName}.sorted.bam 

  }

  output {
      File rawBam = "${sampleName}.sorted.bam"
      File rawBamIndex = "${sampleName}.sorted.bam.bai"
  }

  runtime {
    docker: "us.gcr.io/cds-docker-containers/gulatide/long-read-align"

  }
}

