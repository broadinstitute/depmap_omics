version 1.0

# Copyright (c) 2017 Leiden University Medical Center
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

task BgzipAndIndex {
    input {
        File inputFile
        String outputDir
        String type = "vcf"

        String memory = "2GiB"
        Int timeMinutes = 1 + ceil(size(inputFile, "GiB"))
        String dockerImage = "quay.io/biocontainers/tabix:0.2.6--ha92aebf_0"
    }

    String outputGz = outputDir + "/" + basename(inputFile) + ".gz"

    command {
        set -e
        mkdir -p "$(dirname ~{outputGz})"
        bgzip -c ~{inputFile} > ~{outputGz}
        tabix ~{outputGz} -p ~{type}
    }

    output {
        File compressed = outputGz
        File index = outputGz + ".tbi"
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        inputFile: {description: "The file to be compressed and indexed.", category: "required"}
        outputDir: {description: "The directory in which the output will be placed.", category: "required"}
        type: {description: "The type of file (eg. vcf or bed) to be compressed and indexed.", category: "common"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        compressed: {description: "Compressed input file."}
        index: {description: "Index of the compressed input file."}
    }
}

task DictAndFaidx {
    input {
        File inputFile
        String javaXmx = "2G"
        String memory = "3GiB"
        Int timeMinutes = 5 + ceil(size(inputFile, "GiB") * 5)
        String dockerImage = "quay.io/biocontainers/samtools:1.16.1--h6899075_1"
    }

    String outputFile = basename(inputFile)
    # Capture .fa¸ .fna and .fasta
    String outputDict = sub(outputFile, "\.fn?as?t?a?$", "") + ".dict"
    # This executes both dict and faidx, so indexes are co-located in the same folder.
    command <<<
        set -e
        cp ~{inputFile} ~{outputFile}
        samtools dict -o ~{outputDict}  ~{outputFile}
        samtools faidx ~{outputFile} --fai-idx ~{outputFile}.fai
    >>>

    output {
        File outputFasta = outputFile
        File outputFastaDict = outputDict
        File outputFastaFai = outputFile + ".fai"
    }

    runtime {
        memory: memory
        docker: dockerImage
        time_minutes: timeMinutes
        cpu: 1
    }

    parameter_meta {
        # inputs
        inputFile: {description: "The input fasta file.", category: "required"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.", category: "advanced"}
        memory: {description: "The amount of memory available to the job.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}
        # outputs
        outputFasta: {description: "Fasta file that is co-located with the indexes"}
        outputFastaFai: {description: "Fasta index file for the outputFasta file."}
        outputFastaDict: {description: "Sequence dictionary for the outputFasta file."}
    }
}

task Faidx {
    input {
        File inputFile
        String outputDir

        String memory = "2GiB"
        String dockerImage = "quay.io/biocontainers/samtools:1.16.1--h6899075_1"
    }

    command {
        set -e
        mkdir -p "~{outputDir}"
        ln -s ~{inputFile} "~{outputDir}/$(basename ~{inputFile})"
        samtools faidx \
        "~{outputDir}/$(basename ~{inputFile})"
    }

    output {
        File outputIndex = outputDir + "/" + basename(inputFile) + ".fai"
    }

    runtime {
        memory: memory
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        inputFile: {description: "The input fasta file.", category: "required"}
        outputDir: {description: "Output directory path.", category: "required"}
        memory: {description: "The amount of memory available to the job.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        outputIndex: {description: "Index of the input fasta file."}
    }
}

task Fastq {
    input {
        File inputBam
        String outputRead1
        String? outputRead2
        String? outputRead0
        Boolean appendReadNumber = false
        Boolean outputQuality = false

        Int? includeFilter
        Int? excludeFilter
        Int? excludeSpecificFilter
        Int? compressionLevel

        Int threads = 1
        String memory = "1GiB"
        Int timeMinutes = 1 + ceil(size(inputBam) * 2)
        String dockerImage = "quay.io/biocontainers/samtools:1.16.1--h6899075_1"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputRead1})"
        samtools fastq \
        ~{true="-1" false="-s" defined(outputRead2)} ~{outputRead1} \
        ~{"-2 " + outputRead2} \
        ~{"-0 " + outputRead0} \
        ~{"-f " + includeFilter} \
        ~{"-F " + excludeFilter} \
        ~{"-G " + excludeSpecificFilter} \
        ~{true="-N" false="-n" appendReadNumber} \
        ~{true="-O" false="" outputQuality} \
        ~{"-c " + compressionLevel} \
        ~{"--threads " + threads} \
        ~{inputBam}
    }

    output {
        File read1 = outputRead1
        File? read2 = outputRead2
        File? read0 = outputRead0
    }

    runtime {
        cpu: threads
        memory: memory
        docker: dockerImage
        time_minutes: timeMinutes
    }

    parameter_meta {
        # inputs
        inputBam: {description: "The bam file to process.", category: "required"}
        outputRead1: {description: "The location the reads (first reads for pairs, in case of paired-end sequencing) should be written to.", category: "required"}
        outputRead2: {description: "The location the second reads from pairs should be written to.", category: "common"}
        outputRead0: {description: "The location the unpaired reads should be written to (in case of paired-end sequenicng).", category: "advanced"}
        appendReadNumber: {description: "Append /1 and /2 to the read name, or don't. Corresponds to `-n/N`.", category: "advanced"}
        outputQuality: {description: "Equivalent to samtools fastq's `-O` flag.", category: "advanced"}
        includeFilter: {description: "Include reads with ALL of these flags. Corresponds to `-f`.", category: "advanced"}
        excludeFilter: {description: "Exclude reads with ONE OR MORE of these flags. Corresponds to `-F`.", category: "advanced"}
        excludeSpecificFilter: {description: "Exclude reads with ALL of these flags. Corresponds to `-G`.", category: "advanced"}
        compressionLevel: {description: "Set compression level when writing gz or bgzf fastq files.", category: "advanced"}
        threads: {description: "The number of threads to use.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        read1: {description: "Reads with the READ1 FLAG set."}
        read2: {description: "Reads with the READ2 FLAG set."}
        read0: {description: "Reads with either READ1 FLAG or READ2 flag set."}
    }
}

task FilterShortReadsBam {
    input {
        File bamFile
        String outputPathBam

        String memory = "1GiB"
        Int timeMinutes = 1 + ceil(size(bamFile, "GiB") * 8)
        String dockerImage = "quay.io/biocontainers/samtools:1.16.1--h6899075_1"
    }

    String outputPathBamIndex = sub(outputPathBam, "\.bam$", ".bai")

    command {
        set -e
        mkdir -p "$(dirname ~{outputPathBam})"
        samtools view -h ~{bamFile} | \
        awk 'length($10) > 30 || $1 ~/^@/' | \
        samtools view -bS -> ~{outputPathBam}
        samtools index ~{outputPathBam} ~{outputPathBamIndex}
    }

    output {
        File filteredBam = outputPathBam
        File filteredBamIndex = outputPathBamIndex
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        bamFile: {description: "The bam file to process.", category: "required"}
        outputPathBam: {description: "The filtered bam file.", category: "common"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        filteredBam: {description: "BAM file filtered for short reads."}
        filteredBamIndex: {description: "Index of filtered BAM file."}
    }
}

task Flagstat {
    input {
        File inputBam
        String outputPath

        String memory = "256MiB"  # Only 40.5 MiB used for 150G bam file.
        Int timeMinutes = 1 + ceil(size(inputBam, "G"))
        String dockerImage = "quay.io/biocontainers/samtools:1.16.1--h6899075_1"
	# biocontainers/samtools
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputPath})"
        samtools flagstat ~{inputBam} > ~{outputPath}
    }

    output {
        File flagstat = outputPath
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
	bootDiskSizeGb: 50
        disks: "local-disk 100 HDD"
    }

    parameter_meta {
        # inputs
        inputBam: {description: "The BAM file for which statistics should be retrieved.", category: "required"}
        outputPath: {description: "The location the ouput should be written to.", category: "required"}
        memory: {description: "The amount of memory needed for the job.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        flagstat: {description: "The number of alignments for each FLAG type."}
    }
}

task Index {
    input {
        File bamFile

        String? outputBamPath

        String memory = "2GiB"
        Int timeMinutes = 1 + ceil(size(bamFile, "GiB") * 4)
        String dockerImage = "quay.io/biocontainers/samtools:1.16.1--h6899075_1"
    }

    # Select_first is needed, otherwise womtool validate fails.
    String outputPath = select_first([outputBamPath, basename(bamFile)])
    String bamIndexPath = sub(outputPath, "\.bam$", ".bai")

    command {
        bash -c '
        set -e
        # Make sure outputBamPath does not exist.
        if [ ! -f ~{outputPath} ]
        then
            mkdir -p "$(dirname ~{outputPath})"
            ln ~{bamFile} ~{outputPath} || cp ~{bamFile} ~{outputPath}
        fi
        samtools index ~{outputPath} ~{bamIndexPath}
        '
    }

    output {
        File indexedBam = outputPath
        File index =  bamIndexPath
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        bamFile: {description: "The BAM file for which an index should be made.", category: "required"}
        outputBamPath: {description: "The location where the BAM file should be written to. The index will appear alongside this link to the BAM file.", category: "common"}
        memory: {description: "The amount of memory needed for the job.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        indexedBam: {description: "BAM file that was indexed."}
        index: {description: "Index of the input BAM file."}
    }
}

task Markdup {
    input {
        File inputBam
        String outputBamPath

        Int timeMinutes = 1 + ceil(size(inputBam, "GiB") * 2)
        String dockerImage = "quay.io/biocontainers/samtools:1.16.1--h6899075_1"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputBamPath})"
        samtools markdup ~{inputBam} ~{outputBamPath}
    }

    output {
        File outputBam = outputBamPath
    }

    runtime {
        docker: dockerImage
        time_minutes: timeMinutes
    }

    parameter_meta {
        # inputs
        inputBam: {description: "The BAM file to be processed.", category: "required"}
        outputBamPath: {description: "The location of the output BAM file.", category: "required"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        outputBam: {description: "BAM file with duplicate alignments marked."}
    }
}

task Merge {
    input {
        Array[File]+ bamFiles
        String outputBamPath = "merged.bam"
        Boolean force = true

        Int threads = 1
        String memory = "4GiB"
        Int timeMinutes = 1 + ceil(size(bamFiles, "GiB") * 2)
        String dockerImage = "quay.io/biocontainers/samtools:1.16.1--h6899075_1"
    }

    String indexPath = sub(outputBamPath, "\.bam$",".bai")

    # Samtools uses additional threads for merge.
    command {
        set -e
        mkdir -p "$(dirname ~{outputBamPath})"
        samtools merge \
        --threads ~{threads - 1} \
        ~{true="-f" false="" force} \
        ~{outputBamPath} ~{sep=' ' bamFiles}
        samtools index ~{outputBamPath} ~{indexPath}
    }

    output {
        File outputBam = outputBamPath
        File outputBamIndex = indexPath
    }

    runtime {
        cpu: threads
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        bamFiles: {description: "The BAM files to merge.", category: "required"}
        outputBamPath: {description: "The location the merged BAM file should be written to.", category: "common"}
        force: {description: "Equivalent to samtools merge's `-f` flag.", category: "advanced"}
        threads: {description: "Number of threads to use.", category: "common"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        outputBam: {description: "Multiple BAM files merged into one."}
        outputBamIndex: {description: "Index of the merged BAM file."}
    }
}

task Sort {
    input {
        File inputBam
        String outputPath = basename(inputBam, "\.bam") + ".sorted.bam"
        Boolean sortByName = false
        Int compressionLevel = 1

        Int memoryPerThreadGb = 4
        Int threads = 1
        Int memoryGb = 1 + threads * memoryPerThreadGb
        Int timeMinutes = 1 + ceil(size(inputBam, "GiB") * 3)
        String dockerImage = "quay.io/biocontainers/samtools:1.16.1--h6899075_1"
    }

    # Select first needed as outputPath is optional input (bug in cromwell).
    String bamIndexPath = sub(select_first([outputPath]), "\.bam$", ".bai")

    command {
        set -e
        mkdir -p "$(dirname ~{outputPath})"
        samtools sort \
        -l ~{compressionLevel} \
        ~{true="-n" false="" sortByName} \
        ~{"--threads " + threads} \
        -m ~{memoryPerThreadGb}G \
        -o ~{outputPath} \
        ~{inputBam}
        samtools index \
        -@ ~{threads} \
        ~{outputPath} ~{bamIndexPath}
    }

    output {
        File outputBam = outputPath
        File outputBamIndex = bamIndexPath
    }

    runtime {
        cpu: threads
        memory: "~{memoryGb}GiB"
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        inputBam: {description: "The input SAM file.", category: "required"}
        outputPath: {description: "Output directory path + output file.", category: "required"}
        sortByName: {description: "Sort the inputBam by read name instead of position.", category: "advanced"}
        compressionLevel: {description: "Compression level from 0 (uncompressed) to 9 (best).", category: "advanced"}
        memoryPerThreadGb: {description: "The amount of memory used per sort thread in gigabytes.", category: "advanced"}
        threads: {description: "The number of additional threads that will be used for this task.", category: "advanced"}
        memoryGb: {description: "The amount of memory available to the job in gigabytes.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        outputBam: {description: "Sorted BAM file."}
        outputBamIndex: {description: "Sorted BAM file index."}
    }
}

task Tabix {
    input {
        File inputFile
        String outputFilePath = basename(inputFile)
        String type = "vcf"

        Int timeMinutes = 1 + ceil(size(inputFile, "GiB") * 2)
        String dockerImage = "quay.io/biocontainers/tabix:0.2.6--ha92aebf_0"
    }

    # FIXME: It is better to do the indexing on VCF creation.
    # Not in a separate task. With file localization this gets hairy fast.
    command {
        set -e
        mkdir -p "$(dirname ~{outputFilePath})"
        if [ ! -f ~{outputFilePath} ]
        then
            ln ~{inputFile} ~{outputFilePath} || cp ~{inputFile} ~{outputFilePath}
        fi
        tabix ~{outputFilePath} -p ~{type}
    }

    output {
        File indexedFile = outputFilePath
        File index = outputFilePath + ".tbi"
    }

    runtime {
        memory: "2GiB"
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        inputFile: {description: "The file to be indexed.", category: "required"}
        outputFilePath: {description: "The location where the file should be written to. The index will appear alongside this link to the file.", category: "common"}
        type: {description: "The type of file (eg. vcf or bed) to be indexed.", category: "common"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        indexedFile: {description: "Indexed input file."}
        index: {description: "Index of the input file."}
    }
}

task View {
    input {
        File inFile
        String outputFileName = "view.bam"
        Boolean uncompressedBamOutput = false

        File? referenceFasta
        Int? includeFilter
        Int? excludeFilter
        Int? excludeSpecificFilter
        Int? MAPQthreshold
        File? targetFile

        Int threads = 1
        String memory = "1GiB"
        Int timeMinutes = 1 + ceil(size(inFile, "GiB") * 5)
        String dockerImage = "quay.io/biocontainers/samtools:1.16.1--h6899075_1"
    }

    String outputIndexPath = basename(outputFileName) + ".bai"

    # Always output to bam and output header.
    command {
        set -e
        mkdir -p "$(dirname ~{outputFileName})"
        samtools view -b \
        ~{"-T " + referenceFasta} \
        ~{"-o " + outputFileName} \
        ~{true="-u " false="" uncompressedBamOutput} \
        ~{"-f " + includeFilter} \
        ~{"-F " + excludeFilter} \
        ~{"-G " + excludeSpecificFilter} \
        ~{"-q " + MAPQthreshold} \
        ~{"--threads " + (threads - 1)} \
        ~{"--target-file " + targetFile} \
        ~{inFile}
        samtools index ~{outputFileName} ~{outputIndexPath}
    }

    output {
        File outputBam = outputFileName
        File outputBamIndex = outputIndexPath
    }

    runtime {
        cpu: threads
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        inFile: {description: "A BAM, SAM or CRAM file.", category: "required"}
        outputFileName: {description: "The location the output BAM file should be written.", category: "common"}
        uncompressedBamOutput: {description: "Equivalent to samtools view's `-u` flag.", category: "advanced"}
        referenceFasta: {description: "The reference fasta file also used for mapping.", category: "advanced"}
        includeFilter: {description: "Equivalent to samtools view's `-f` option.", category: "advanced"}
        excludeFilter: {description: "Equivalent to samtools view's `-F` option.", category: "advanced"}
        excludeSpecificFilter: {description: "Equivalent to samtools view's `-G` option.", category: "advanced"}
        MAPQthreshold: {description: "Equivalent to samtools view's `-q` option.", category: "advanced"}
        targetFile: {description: "A BED file with regions to include", caegory: "advanced"}
        threads: {description: "The number of threads to use.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        outputBam: {description: "Processed input file."}
        outputBamIndex: {description: "Index of the processed input file."}
    }
}

task ViewCount {
    input {
        File inputBam

        Int? MAPQthreshold

        String memory = "1GiB"
        Int timeMinutes = 1 + ceil(size(inputBam, "GiB") * 5)
        String dockerImage = "quay.io/biocontainers/samtools:1.16.1--h6899075_1"
    }

    # Always output to bam and output header.
    command <<<
        samtools view -c \
        ~{"-q " + MAPQthreshold} \
        ~{inputBam}
    >>>

    output {
        File count_stat = stdout()
    }

    runtime {
        cpu: 1
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

task Depth {
    input {
        File inputBam
        File bamIndex

        String memory = "1GiB"
        Int timeMinutes = 1 + ceil(size(inputBam, "GiB") * 5)
        String dockerImage = "quay.io/biocontainers/samtools:1.16.1--h6899075_1"
    }

    # Always output to bam and output header.
    command <<<
        samtools depth ~{inputBam} |  awk '{sum+=$3} END { print "Average = ",sum/NR}'
    >>>
    output {
        File result = stdout()
    }

    runtime {
        cpu: 1
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

