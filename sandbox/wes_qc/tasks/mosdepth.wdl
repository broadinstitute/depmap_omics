version 1.0


task Depth {
    input {
        File inputBam
        File bamIndex
	File targetBed
        String outputPath

        Int? MAPQthreshold
        Int threads = 3

        String memory = "2GiB"
        Int timeMinutes = 1 + ceil(size(inputBam, "GiB") * 5)
        String dockerImage = "quay.io/biocontainers/mosdepth:0.3.3--h37c5b7d_2"
    }

    # Always output to bam and output header.
    command <<<
        mosdepth -n --fast-mode -t ~{threads} -b ~{targetBed} ~{outputPath} ~{inputBam}
    >>>

    output {
        File result = outputPath + ".mosdepth.region.dist.txt"
        File result2 = outputPath + ".mosdepth.global.dist.txt"
        File result_sum = outputPath + ".mosdepth.summary.txt"
    }

    runtime {
        cpu: threads
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
        bootDiskSizeGb: 50
        disks: "local-disk 100 HDD"
    }

    parameter_meta {
        # inputs
        inputBam: {description: "A BAM, SAM or CRAM file.", category: "required"}
    }
}
