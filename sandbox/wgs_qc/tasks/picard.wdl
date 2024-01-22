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

task BedToIntervalList {
    input {
        File bedFile
        File dict
        String outputPath = "regions.interval_list"

        String javaXmx = "3G"
        String memory = "4GiB"
        Int timeMinutes = 5
        String dockerImage = "quay.io/biocontainers/picard:2.26.10--hdfd78af_0"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputPath})"
        picard -Xmx~{javaXmx} -XX:ParallelGCThreads=1 \
        BedToIntervalList \
        I=~{bedFile} \
        O=~{outputPath} \
        SD=~{dict}
    }

    output {
        File intervalList = outputPath
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        bedFile: {description: "A bed file.", category: "required"}
        dict: {description: "A sequence dict file.", category: "required"}
        outputPath: {description: "The location the output interval list should be written to.", category: "advanced"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        intervalList: {description: "Picard Interval List from a BED file."}
    }
}

task CollectHsMetrics {
    input {
        File inputBam
        File inputBamIndex
        File referenceFasta
        File referenceFastaDict
        File referenceFastaFai
        File targets
        String basename

        File? baits

        # Use the targets file as baits as a fallback, since often the baits
        # for a certain capture kit are not available.
        File baitsFile = select_first([baits, targets])
        File targetsFile = targets

        Int javaXmxMb = 3072
        Int memoryMb = javaXmxMb + 512
        # Additional * 2 because picard multiple metrics reads the
        # reference fasta twice.
        Int timeMinutes = 1 + ceil(size(referenceFasta, "GiB") * 3 * 2) + ceil(size(inputBam, "GiB") * 6)
        String dockerImage = "quay.io/biocontainers/picard:2.26.10--hdfd78af_0"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{basename})"

        picard -Xmx~{javaXmxMb}M -XX:ParallelGCThreads=1 \
        BedToIntervalList \
        I=~{targetsFile} \
        O=~{basename}.intervals \
        SD=~{referenceFastaDict}

        picard -Xmx~{javaXmxMb}M -XX:ParallelGCThreads=1 \
        CollectHsMetrics \
        I=~{inputBam} \
        R=~{referenceFasta} \
        BAIT_INTERVALS=~{basename}.intervals \
        TARGET_INTERVALS=~{basename}.intervals \
        O="~{basename}.hs_metrics.txt"
    }

    output {
        File HsMetrics = basename + ".hs_metrics.txt"
    }

    runtime {
        memory: "~{memoryMb}MiB"
        time_minutes: timeMinutes
        docker: dockerImage
	bootDiskSizeGb: 100
        disks: "local-disk 100 HDD"
    }

    parameter_meta {
        # inputs
        inputBam: {description: "The input BAM file for which metrics will be collected.", category: "required"}
        inputBamIndex: {description: "The index of the input BAM file.", category: "required"}
        referenceFasta: {description: "The reference fasta file which was also used for mapping.", category: "required"}
        referenceFastaDict: {description: "The sequence dictionary associated with the reference fasta file.", category: "required"}
        referenceFastaFai: {description: "The index for the reference fasta file.", category: "required"}
        targets: {description: "Picard interval file of the capture targets.", category: "required"}
        targetsFile: {description: "Picard interval file of the capture targets, the same as targets.", category: "advanced"}
        basename: {description: "The basename/prefix of the output files (may include directories).", category: "required"}
        baits: {description: "Picard interval file of the capture bait set.", category: "advanced"}
        baitsFile: {description: "Picard interval file of the bait set. Uses targets as a fallback when baits is not set.", category: "advanced"}
        javaXmxMb: {description: "The maximum memory available to the program in megabytes. Should be lower than `memoryMb` to accommodate JVM overhead.", category: "advanced"}
        memoryMb: {description: "The amount of memory this job will use in megabytes.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        HsMetrics: {description: "Hybrid-selection (HS) metrics for the input BAM file."}
    }
}

task CollectInsertSizeMetrics {
    input {
        File inputBam
        File inputBamIndex

        Float? minimumPercentage
        String basename = "./insertSize_metrics"

        String memory = "5GiB"
        String javaXmx = "4G"
        Int timeMinutes = 1 + ceil(size(inputBam, "GiB") * 6)
        String dockerImage = "quay.io/biocontainers/picard:2.23.2--0"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{basename})"
        picard -Xmx~{javaXmx} -XX:ParallelGCThreads=1 \
        CollectInsertSizeMetrics \
        I=~{inputBam} \
        O=~{basename}.txt \
        H=~{basename}.pdf \
        ~{"M=" + minimumPercentage}
    }

    output {
        File metricsTxt = "~{basename}.txt"
        File metricsPdf = "~{basename}.pdf"
    }

    runtime {
        docker: dockerImage
        time_minutes: timeMinutes
        memory: memory
    }
    
    parameter_meta {
        # inputs
        inputBam: {description: "The input BAM file for which metrics will be collected.", category: "required"}
        inputBamIndex: {description: "The index of the input BAM file.", category: "required"}
        minimumPercentage: {description: "Equivalent to picard CollectInsertSizeMetrics' `M` option.", category: "advanced"}
        basename: {description: "The basename for the output files.", category: "common"}

        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.",
                  category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}

task CollectMultipleMetrics {
    input {
        File inputBam
        File inputBamIndex
        File referenceFasta
        File referenceFastaDict
        File referenceFastaFai
        String basename
        Boolean collectAlignmentSummaryMetrics = true
        Boolean collectInsertSizeMetrics = true
        Boolean qualityScoreDistribution = true
        Boolean meanQualityByCycle = true
        Boolean collectBaseDistributionByCycle = true
        Boolean collectGcBiasMetrics = true
        #FIXME: Boolean rnaSeqMetrics = false # There is a bug in picard https://github.com/broadinstitute/picard/issues/999
        Boolean collectSequencingArtifactMetrics = true
        Boolean collectQualityYieldMetrics = true

        Int javaXmxMb = 3072
        Int memoryMb = javaXmxMb + 512
        # Additional * 2 because picard multiple metrics reads the reference fasta twice.
        Int timeMinutes = 1 + ceil(size(referenceFasta, "GiB") * 3 * 2) + ceil(size(inputBam, "GiB") * 6)
        String dockerImage = "quay.io/biocontainers/picard:2.26.10--hdfd78af_0"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{basename})"
        picard -Xmx~{javaXmxMb}M -XX:ParallelGCThreads=1 \
        CollectMultipleMetrics \
        I=~{inputBam} \
        R=~{referenceFasta} \
        O=~{basename} \
        PROGRAM=null \
        ~{true="PROGRAM=CollectAlignmentSummaryMetrics" false="" collectAlignmentSummaryMetrics} \
        ~{true="PROGRAM=CollectInsertSizeMetrics" false="" collectInsertSizeMetrics} \
        ~{true="PROGRAM=QualityScoreDistribution" false="" qualityScoreDistribution} \
        ~{true="PROGRAM=MeanQualityByCycle" false="" meanQualityByCycle} \
        ~{true="PROGRAM=CollectBaseDistributionByCycle" false="" collectBaseDistributionByCycle} \
        ~{true="PROGRAM=CollectGcBiasMetrics" false="" collectGcBiasMetrics} \
        ~{true="PROGRAM=CollectSequencingArtifactMetrics" false="" collectSequencingArtifactMetrics} \
        ~{true="PROGRAM=CollectQualityYieldMetrics" false="" collectQualityYieldMetrics}
    }

    output {
        File? alignmentSummary = basename + ".alignment_summary_metrics"
        File? baitBiasDetail = basename + ".bait_bias_detail_metrics"
        File? baitBiasSummary = basename + ".bait_bias_summary_metrics"
        File? baseDistributionByCycle = basename + ".base_distribution_by_cycle_metrics"
        File? baseDistributionByCyclePdf = basename + ".base_distribution_by_cycle.pdf"
        File? errorSummary = basename + ".error_summary_metrics"
        File? gcBiasDetail = basename + ".gc_bias.detail_metrics"
        File? gcBiasPdf = basename + ".gc_bias.pdf"
        File? gcBiasSummary = basename + ".gc_bias.summary_metrics"
        File? insertSizeHistogramPdf = basename + ".insert_size_histogram.pdf"
        File? insertSize = basename + ".insert_size_metrics"
        File? preAdapterDetail = basename + ".pre_adapter_detail_metrics"
        File? preAdapterSummary = basename + ".pre_adapter_summary_metrics"
        File? qualityByCycle = basename + ".quality_by_cycle_metrics"
        File? qualityByCyclePdf = basename + ".quality_by_cycle.pdf"
        File? qualityDistribution = basename + ".quality_distribution_metrics"
        File? qualityDistributionPdf = basename + ".quality_distribution.pdf"
        File? qualityYield = basename + ".quality_yield_metrics"
        # Using a glob is easier. But will lead to very ugly output directories.
        Array[File] allStats = select_all([
            alignmentSummary,
            baitBiasDetail,
            baitBiasSummary,
            baseDistributionByCycle,
            baseDistributionByCyclePdf,
            errorSummary,
            gcBiasDetail,
            gcBiasPdf,
            gcBiasSummary,
            insertSizeHistogramPdf,
            insertSize,
            preAdapterDetail,
            preAdapterSummary,
            qualityByCycle,
            qualityByCyclePdf,
            qualityDistribution,
            qualityDistributionPdf,
            qualityYield
        ])
    }

    runtime {
        memory: "~{memoryMb}MiB"
        time_minutes: timeMinutes
        docker: dockerImage
	bootDiskSizeGb: 100
        disks: "local-disk 100 HDD"
    }

    parameter_meta {
        # inputs
        inputBam: {description: "The input BAM file for which metrics will be collected.", category: "required"}
        inputBamIndex: {description: "The index of the input BAM file.", category: "required"}
        referenceFasta: {description: "The reference fasta file which was also used for mapping.", category: "required"}
        referenceFastaDict: {description: "The sequence dictionary associated with the reference fasta file.", category: "required"}
        referenceFastaFai: {description: "The index for the reference fasta file.", category: "required"}
        basename: {description: "The basename/prefix of the output files (may include directories).", category: "required"}
        collectAlignmentSummaryMetrics: {description: "Equivalent to the `PROGRAM=CollectAlignmentSummaryMetrics` argument.", category: "advanced"}
        collectInsertSizeMetrics: {description: "Equivalent to the `PROGRAM=CollectInsertSizeMetrics` argument.", category: "advanced"}
        qualityScoreDistribution: {description: "Equivalent to the `PROGRAM=QualityScoreDistribution` argument.", category: "advanced"}
        meanQualityByCycle: {description: "Equivalent to the `PROGRAM=MeanQualityByCycle` argument.", category: "advanced"}
        collectBaseDistributionByCycle: {description: "Equivalent to the `PROGRAM=CollectBaseDistributionByCycle` argument.", category: "advanced"}
        collectGcBiasMetrics: {description: "Equivalent to the `PROGRAM=CollectGcBiasMetrics` argument.", category: "advanced"}
        collectSequencingArtifactMetrics: {description: "Equivalent to the `PROGRAM=CollectSequencingArtifactMetrics` argument.", category: "advanced"}
        collectQualityYieldMetrics: {description: "Equivalent to the `PROGRAM=CollectQualityYieldMetrics` argument.", category: "advanced"}
        javaXmxMb: {description: "The maximum memory available to the program in megabytes. Should be lower than `memoryMb` to accommodate JVM overhead.", category: "advanced"}
        memoryMb: {description: "The amount of memory this job will use in megabytes.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        alignmentSummary: {description: ""}
        baitBiasDetail: {description: ""}
        baitBiasSummary: {description: ""}
        baseDistributionByCycle: {description: ""}
        baseDistributionByCyclePdf: {description: ""}
        errorSummary: {description: ""}
        gcBiasDetail: {description: ""}
        gcBiasPdf: {description: ""}
        gcBiasSummary: {description: ""}
        insertSizeHistogramPdf: {description: ""}
        insertSize: {description: ""}
        preAdapterDetail: {description: ""}
        preAdapterSummary: {description: ""}
        qualityByCycle: {description: ""}
        qualityByCyclePdf: {description: ""}
        qualityDistribution: {description: ""}
        qualityDistributionPdf: {description: ""}
        qualityYield: {description: ""}
        allStats: {description: ""}
    }
}

task CollectRnaSeqMetrics {
    input {
        File inputBam
        File inputBamIndex
        File refRefflat
        String basename
        String strandSpecificity = "NONE"

        String javaXmx =  "8G"
        String memory = "9GiB"
        # With 6 minutes per G there were several timeouts.
        Int timeMinutes = 1 + ceil(size(inputBam, "GiB") * 12)
        String dockerImage = "quay.io/biocontainers/picard:2.26.10--hdfd78af_0"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{basename})"
        picard -Xmx~{javaXmx} \
        CollectRnaSeqMetrics -XX:ParallelGCThreads=1 \
        I=~{inputBam} \
        O=~{basename}.RNA_Metrics \
        CHART_OUTPUT=~{basename}.RNA_Metrics.pdf \
        STRAND_SPECIFICITY=~{strandSpecificity} \
        REF_FLAT=~{refRefflat}
    }

    output {
        File metrics = basename + ".RNA_Metrics"
        File? chart = basename + ".RNA_Metrics.pdf"
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        inputBam: {description: "The input BAM file for which metrics will be collected.", category: "required"}
        inputBamIndex: {description: "The index of the input BAM file.", category: "required"}
        refRefflat: {description: "A refflat file containing gene annotations.", catehory: "required"}
        basename: {description: "The basename/prefix of the output files (may include directories).", category: "required"}
        strandSpecificity: {description: "Equivalent to the `STRAND_SPECIFICITY` option of picard's CollectRnaSeqMetrics.", category: "common"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        metrics: {description: "Metrics describing the distribution of bases within the transcripts."}
        chart: {description: "Plot of normalized position vs. coverage."}
    }
}

task CollectTargetedPcrMetrics {
    input {
        File inputBam
        File inputBamIndex
        File referenceFasta
        File referenceFastaDict
        File referenceFastaFai
        File ampliconIntervals
        Array[File]+ targetIntervals
        String basename

        String javaXmx = "3G"
        String memory = "4GiB"
        Int timeMinutes = 1 + ceil(size(inputBam, "GiB") * 6)
        String dockerImage = "quay.io/biocontainers/picard:2.26.10--hdfd78af_0"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{basename})"
        picard -Xmx~{javaXmx} -XX:ParallelGCThreads=1 \
        CollectTargetedPcrMetrics \
        I=~{inputBam} \
        R=~{referenceFasta} \
        AMPLICON_INTERVALS=~{ampliconIntervals} \
        TARGET_INTERVALS=~{sep=" TARGET_INTERVALS=" targetIntervals} \
        O=~{basename}.targetPcrMetrics \
        PER_BASE_COVERAGE=~{basename}.targetPcrPerBaseCoverage \
        PER_TARGET_COVERAGE=~{basename}.targetPcrPerTargetCoverage
    }

    output {
        File perTargetCoverage = basename + ".targetPcrPerTargetCoverage"
        File perBaseCoverage = basename + ".targetPcrPerBaseCoverage"
        File metrics = basename + ".targetPcrMetrics"
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        inputBam: {description: "The input BAM file for which metrics will be collected.", category: "required"}
        inputBamIndex: {description: "The index of the input BAM file.", category: "required"}
        referenceFasta: {description: "The reference fasta file which was also used for mapping.", category: "required"}
        referenceFastaDict: {description: "The sequence dictionary associated with the reference fasta file.", category: "required"}
        referenceFastaFai: {description: "The index for the reference fasta file.", category: "required"}
        ampliconIntervals: {description: "An interval list describinig the coordinates of the amplicons sequenced.", category: "required"}
        targetIntervals: {description: "An interval list describing the coordinates of the targets sequenced.", category: "required"}
        basename: {description: "The basename/prefix of the output files (may include directories).", category: "required"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        perTargetCoverage: {description: "Per target coverage information."}
        perBaseCoverage: {description: "Per base coverage information to."}
        metrics: {description: "File containing metrics."}
    }
}

task CollectVariantCallingMetrics {
    input {
        File dbsnp
        File dbsnpIndex
        File inputVCF
        File inputVCFIndex
        String basename

        String javaXmx =  "8G"
        String memory = "9GiB"
        Int timeMinutes = 1440
        String dockerImage = "quay.io/biocontainers/picard:2.26.10--hdfd78af_0"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{basename})"
        picard -Xmx~{javaXmx} \
        CollectVariantCallingMetrics -XX:ParallelGCThreads=1 \
        DBSNP=~{dbsnp} \
        INPUT=~{inputVCF} \
        OUTPUT=~{basename}
    }

    output {
        File details = basename + ".variant_calling_detail_metrics"
        File summary = basename + ".variant_calling_summary_metrics"
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        dbsnp: {description: "DBSNP vcf file to use with CollectVariantCallingMetrics.", category: "required"}
        dbsnpIndex: {description: "Index file for the DBSNP VCF.", category: "required"}
        inputVCF: {description: "Input VCF file.", category: "required"}
        inputVCFIndex: {description: "Index file for the input VCF.", category: "required"}
        basename: {description: "The basename/prefix of the output files (may include directories).", category: "required"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        details: {description: ""}
        summary: {description: ""}
    }
}

task CollectWgsMetrics {
    input {
        File inputBam
        File inputBamIndex
        File referenceFasta
        File referenceFastaDict
        File referenceFastaFai
        String outputPath = "./wgs_metrics.txt"
        
        Int? minimumMappingQuality
        Int? minimumBaseQuality
        Int? coverageCap

        String memory = "25GiB"
        String javaXmx = "8G"
        Int timeMinutes = 1 + ceil(size(inputBam, "GiB") * 10)
        String dockerImage = "quay.io/biocontainers/picard:2.23.2--0"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputPath})"

        picard -Xmx~{javaXmx} -XX:ParallelGCThreads=1 \
        CollectWgsMetrics \
        REFERENCE_SEQUENCE=~{referenceFasta} \
        INPUT=~{inputBam} \
        OUTPUT=~{outputPath} \
        ~{"MINIMUM_MAPPING_QUALITY=" + minimumMappingQuality} \
        ~{"MINIMUM_BASE_QUALITY=" + minimumBaseQuality} \
        ~{"COVERAGE_CAP=" + coverageCap}
    }

    output {
        File metrics = outputPath
    }

    runtime {
        docker: dockerImage
        time_minutes: timeMinutes
        memory: memory
    }
    
    parameter_meta {
        # inputs
        inputBam: {description: "The input BAM file for which metrics will be collected.", category: "required"}
        inputBamIndex: {description: "The index of the input BAM file.", category: "required"}
        referenceFasta: {description: "The reference fasta file which was also used for mapping.", category: "required"}
        referenceFastaDict: {description: "The sequence dictionary associated with the reference fasta file.",
                             category: "required"}
        referenceFastaFai: {description: "The index for the reference fasta file.", category: "required"}
        outputPath: {description: "The path picard CollectWgsMetrics' output should be written to.", category: "common"}
        minimumMappingQuality: {description: "Equivalent to picard CollectWgsMetrics' MINIMUM_MAPPING_QUALITY option.", category: "advanced"}
        minimumBaseQuality: {description: "Equivalent to picard CollectWgsMetrics' MINIMUM_BASE_QUALITY option.", category: "advanced"}
        coverageCap: {description: "Equivalent to picard CollectWgsMetrics' OVERAGE_CAP option.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.",
                  category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}

task CreateSequenceDictionary {
    input {
        File inputFile
        String outputDir

        String javaXmx = "2G"
        String memory = "3GiB"
        String dockerImage = "quay.io/biocontainers/picard:2.26.10--hdfd78af_0"
    }

    command {
        set -e
        mkdir -p "~{outputDir}"
        picard -Xmx~{javaXmx} \
        -XX:ParallelGCThreads=1 \
        CreateSequenceDictionary \
        REFERENCE=~{inputFile} \
        OUTPUT="~{outputDir}/$(basename ~{inputFile}).dict"
    }

    output {
        File outputDict = outputDir + "/" + basename(inputFile) + ".dict"
    }

    runtime {
        memory: memory
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        inputFile: {description: "The input fasta file.", category: "required"}
        outputDir: {description: "Output directory path.", category: "required"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.", category: "advanced"}
        memory: {description: "The amount of memory available to the job.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        outputDict: {description: "Dictionary of the input fasta file."}
    }
}

# Combine multiple recalibrated BAM files from scattered
# ApplyRecalibration runs.
task GatherBamFiles {
    input {
        Array[File]+ inputBams
        Array[File]+ inputBamsIndex
        String outputBamPath
        Boolean createMd5File = false

        Int compressionLevel = 1
        Boolean useJdkInflater = false
        Boolean useJdkDeflater = true  # Achieves much better compression rates than the intel deflater

        Int javaXmxMb = 1024
        Int memoryMb = javaXmxMb + 512
        # One minute per input gigabyte.
        Int timeMinutes = 1 + ceil(size(inputBams, "GiB") * 1)
        String dockerImage = "quay.io/biocontainers/picard:2.26.10--hdfd78af_0"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputBamPath})"
        picard -Xmx~{javaXmxMb}M -XX:ParallelGCThreads=1 \
        GatherBamFiles \
        INPUT=~{sep=' INPUT=' inputBams} \
        OUTPUT=~{outputBamPath} \
        COMPRESSION_LEVEL=~{compressionLevel} \
        USE_JDK_INFLATER=~{true="true" false="false" useJdkInflater} \
        USE_JDK_DEFLATER=~{true="true" false="false" useJdkDeflater} \
        CREATE_INDEX=true \
        CREATE_MD5_FILE=~{true="true" false="false" createMd5File}
    }

    output {
        File outputBam = outputBamPath
        File outputBamIndex = sub(outputBamPath, "\.bam$", ".bai")
        File? outputBamMd5 = outputBamPath + ".md5"
    }

    runtime {
        memory: "~{memoryMb}MiB"
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        inputBams: {description: "The BAM files to be merged together.", category: "required"}
        inputBamsIndex: {description: "The indexes of the input BAM files.", category: "required"}
        outputBamPath: {description: "The path where the merged BAM file will be written.", caregory: "required"}
        createMd5File: {decription: "Whether to create an md5 file of the output BAM.", category: "advanced"}
        compressionLevel: {description: "The compression level at which the BAM files are written.", category: "advanced"}
        useJdkInflater: {description: "True, uses the java inflater. False, uses the optimized intel inflater.", category: "advanced"}
        useJdkDeflater: {description: "True, uses the java deflator to compress the BAM files. False uses the optimized intel deflater.", category: "advanced"}
        javaXmxMb: {description: "The maximum memory available to the program in megabytes. Should be lower than `memoryMb` to accommodate JVM overhead.", category: "advanced"}
        memoryMb: {description: "The amount of memory this job will use in megabytes.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        outputBam: {description: "Concatenated BAM files."}
        outputBamIndex: {description: "Index of the output `outputBam`."}
        outputBamMd5: {description: "MD5 of the output `outputBam`."}
    }
}

task GatherVcfs {
    input {
        Array[File]+ inputVcfs
        Array[File]+ inputVcfIndexes
        String outputVcfPath = "out.vcf.gz"

        Int compressionLevel = 1
        Boolean useJdkInflater = false
        Boolean useJdkDeflater = true  # Achieves much better compression rates than the intel deflater

        String javaXmx = "4G"
        String memory = "5GiB"
        Int timeMinutes = 1 + ceil(size(inputVcfs, "GiB") * 2)
        String dockerImage = "quay.io/biocontainers/picard:2.26.10--hdfd78af_0"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputVcfPath})"
        picard -Xmx~{javaXmx} -XX:ParallelGCThreads=1 \
        GatherVcfs \
        COMPRESSION_LEVEL=~{compressionLevel} \
        USE_JDK_INFLATER=~{true="true" false="false" useJdkInflater} \
        USE_JDK_DEFLATER=~{true="true" false="false" useJdkDeflater} \
        CREATE_INDEX=true \
        INPUT=~{sep=' INPUT=' inputVcfs} \
        OUTPUT=~{outputVcfPath}
    }

    output {
        File outputVcf = outputVcfPath
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        inputVcfs: {description: "The VCF files to be merged together.", category: "required"}
        inputVcfIndexes: {description: "The indexes of the input VCF files.", category: "required"}
        outputVcfPath: {description: "The path where the merged VCF file will be written.", caregory: "required"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        compressionLevel: {description: "The compression level at which the BAM files are written.", category: "advanced"}
        useJdkInflater: {description: "True, uses the java inflater. False, uses the optimized intel inflater.", category: "advanced"}
        useJdkDeflater: {description: "True, uses the java deflator to compress the BAM files. False uses the optimized intel deflater.", category: "advanced"}

        # outputs
        outputVcf: {description: "Multiple VCF files gathered into one file."}
    }
}

# Mark duplicate reads to avoid counting non-independent observations.
task MarkDuplicates {
    input {
        Array[File]+ inputBams
        String outputBamPath
        String metricsPath
        Boolean createMd5File = false
        
        Int compressionLevel = 1
        Boolean useJdkInflater = false
        Boolean useJdkDeflater = true  # Achieves much better compression rates than the intel deflater

        # The program default for READ_NAME_REGEX is appropriate in nearly every case.
        # Sometimes we wish to supply "null" in order to turn off optical duplicate detection.
        # This can be desirable if you don't mind the estimated library size
        # being wrong and optical duplicate detection is taking >7 days and failing.
        String? read_name_regex

        # In GATK Best practices pipeline MarkDuplicates is given a 7G VM.
        # https://github.com/gatk-workflows/broad-prod-wgs-germline-snps-indels/blob/d2934ed656ade44801f9cfe1c0e78d4f80684b7b/PairedEndSingleSampleWf-fc-hg38.wdl#L1040
        Int javaXmxMb =  6656  # 6.5G
        String memoryMb = javaXmxMb + 512

        Int timeMinutes = 1 + ceil(size(inputBams, "GiB") * 8)
        String dockerImage = "quay.io/biocontainers/picard:2.26.10--hdfd78af_0"
    }

    # Task is assuming query-sorted input so that the Secondary and Supplementary reads get
    # marked correctly. This works because the output of BWA is query-grouped and therefore,
    # so is the output of MergeBamAlignment. While query-grouped isn't actually query-sorted,
    # it's good enough for MarkDuplicates with ASSUME_SORT_ORDER="queryname".
    command {
        set -e
        mkdir -p "$(dirname ~{outputBamPath})"
        picard -Xmx~{javaXmxMb}M -XX:ParallelGCThreads=1 \
        MarkDuplicates \
        INPUT=~{sep=' INPUT=' inputBams} \
        OUTPUT=~{outputBamPath} \
        METRICS_FILE=~{metricsPath} \
        COMPRESSION_LEVEL=~{compressionLevel} \
        USE_JDK_INFLATER=~{true="true" false="false" useJdkInflater} \
        USE_JDK_DEFLATER=~{true="true" false="false" useJdkDeflater} \
        VALIDATION_STRINGENCY=SILENT \
        ~{"READ_NAME_REGEX=" + read_name_regex} \
        OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
        CLEAR_DT="false" \
        CREATE_INDEX=true \
        ADD_PG_TAG_TO_READS=false \
        CREATE_MD5_FILE=~{true="true" false="false" createMd5File} \
    }

    output {
        File outputBam = outputBamPath
        File outputBamIndex = sub(outputBamPath, "\.bam$", ".bai")
        File? outputBamMd5 = outputBamPath + ".md5"
        File metricsFile = metricsPath
    }

    runtime {
        memory: "~{memoryMb}MiB"
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        inputBams: {description: "The BAM files for which the duplicate reads should be marked.", category: "required"}
        outputBamPath: {description: "The location where the ouptut BAM file should be written.", category: "required"}
        metricsPath: {description: "The location where the output metrics file should be written.", category: "required"}
        compressionLevel: {description: "The compression level at which the BAM files are written.", category: "advanced"}
        useJdkInflater: {description: "True, uses the java inflater. False, uses the optimized intel inflater.", category: "advanced"}
        useJdkDeflater: {description: "True, uses the java deflator to compress the BAM files. False uses the optimized intel deflater.", category: "advanced"}
        createMd5File: {description: "Whether to create a md5 file for the created BAM file.", category: "advanced"}
        read_name_regex: {description: "Equivalent to the `READ_NAME_REGEX` option of MarkDuplicates.", category: "advanced"}
        javaXmxMb: {description: "The maximum memory available to the program in megabytes. Should be lower than `memoryMb` to accommodate JVM overhead.", category: "advanced"}
        memoryMb: {description: "The amount of memory this job will use in megabytes.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        outputBam: {description: ""}
        outputBamIndex: {description: ""}
        outputBamMd5: {description: ""}
        metricsFile: {description: ""}
    }
}

# Combine multiple VCFs or GVCFs from scattered HaplotypeCaller runs.
task MergeVCFs {
    input {
        Array[File]+ inputVCFs
        Array[File]+ inputVCFsIndexes
        String outputVcfPath
        Int compressionLevel = 1
        Boolean useJdkInflater = false
        # Better results for compression level 1 (much smaller).
        # Higher compression levels similar to intel deflater.
        # NOTE: this might change in the future when the intel deflater is updated!
        # Second NOTE: No it did not change. Only the fastest algorithm with
        # worse compression is wrapped in the intel GKL. Instead of using
        # one of the slightly slower but better compressing alternatives from ISA-L. 
        # (Which are also faster than zlib.)
        Boolean useJdkDeflater = true  # Achieves much better compression rates than the intel deflater

        String javaXmx = "4G"
        String memory = "5GiB"
        Int timeMinutes = 1 + ceil(size(inputVCFs, "GiB")) * 2
        String dockerImage = "quay.io/biocontainers/picard:2.26.10--hdfd78af_0"
    }

    # Using MergeVcfs instead of GatherVcfs so we can create indices.
    # See https://github.com/broadinstitute/picard/issues/789 for relevant GatherVcfs ticket.
    command {
        set -e
        mkdir -p "$(dirname ~{outputVcfPath})"
        picard -Xmx~{javaXmx} -XX:ParallelGCThreads=1 \
        MergeVcfs \
        INPUT=~{sep=' INPUT=' inputVCFs} \
        OUTPUT=~{outputVcfPath} \
        COMPRESSION_LEVEL=~{compressionLevel} \
        USE_JDK_INFLATER=~{true="true" false="false" useJdkInflater} \
        USE_JDK_DEFLATER=~{true="true" false="false" useJdkDeflater}
    }

    output {
        File outputVcf = outputVcfPath
        File outputVcfIndex = outputVcfPath + ".tbi"
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        inputVCFs: {description: "The VCF files to be merged.", category: "required"}
        inputVCFsIndexes: {description: "The indexes of the VCF files.", category: "required"}
        outputVcfPath: {description: "The location the output VCF file should be written to.", category: "required"}
        compressionLevel: {description: "The compression level at which the BAM files are written.", category: "advanced"}
        useJdkInflater: {description: "True, uses the java inflater. False, uses the optimized intel inflater.", category: "advanced"}
        useJdkDeflater: {description: "True, uses the java deflator to compress the BAM files. False uses the optimized intel deflater.", category: "advanced"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        outputVcf: {description: "Multiple variant files combined into a single variant file."}
        outputVcfIndex: {description: "Index of `outputVcf`."}
    }
}

task SamToFastq {
    input {
        File inputBam
        File inputBamIndex
        Boolean paired = true

        String javaXmx = "16G" # High memory default to avoid crashes.
        String memory = "17GiB"
        Int timeMinutes = 30
        String dockerImage = "quay.io/biocontainers/picard:2.26.10--hdfd78af_0"

        File? noneFile
    }

    String outputRead1 = basename(inputBam, "\.[bs]am") + "_R1.fastq.gz"
    String outputRead2 = basename(inputBam, "\.[bs]am") + "_R2.fastq.gz"
    String outputUnpaired = basename(inputBam, "\.[bs]am") + "_unpaired.fastq.gz"

    command {
        set -e
        picard -Xmx~{javaXmx} -XX:ParallelGCThreads=1 \
        SamToFastq \
        I=~{inputBam} \
        ~{"FASTQ=" + outputRead1} \
        ~{if paired then "SECOND_END_FASTQ=" + outputRead2 else ""} \
        ~{if paired then "UNPAIRED_FASTQ=" + outputUnpaired else ""}
    }

    output {
        File read1 = outputRead1
        File? read2 = if paired then outputRead2 else noneFile
        File? unpairedRead = if paired then outputUnpaired else noneFile
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        inputBam: {description: "Input BAM file to extract reads from.", category: "required"}
        inputBamIndex: {description: "Input BAM index file.", category: "required"}
        paired: {description: "Set to false when input data is single-end.", category: "common"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        read1: {description: "Fastq file containing reads from the first pair."}
        read2: {description: "Fastq file containing reads from the second pair."}
        unpairedRead: {description: "Fastq file containing unpaired reads."}
    }

    meta {
        WDL_AID: {
            exclude: ["noneFile"]
        }
    }
}

task ScatterIntervalList {
    input {
        File interval_list
        Int scatter_count

        String javaXmx = "3G"
        String memory = "4GiB"
        String dockerImage = "quay.io/biocontainers/picard:2.26.10--hdfd78af_0"
    }

    command {
        set -e
        mkdir scatter_list
        picard -Xmx~{javaXmx} -XX:ParallelGCThreads=1 \
        IntervalListTools \
        SCATTER_COUNT=~{scatter_count} \
        SUBDIVISION_MODE=BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \
        UNIQUE=true \
        SORT=true \
        INPUT=~{interval_list} \
        OUTPUT=scatter_list
    }

    output {
        Array[File] out = glob("scatter_list/*/*.interval_list")
        Int interval_count = read_int(stdout())
    }

    runtime {
        memory: memory
        docker: dockerImage
    }
}

task SortSam {
    input {
        File inputBam
        String outputPath
        Boolean sortByName = false
        Boolean createMd5File = false
        Int maxRecordsInRam = 500000
        Int compressionLevel = 1
        Boolean useJdkInflater = false
        Boolean useJdkDeflater = true  # Achieves much better compression rates than the intel deflater

        # Default ram of 4 GB. Using 125001.0  to prevent an answer of
        # 4.000000001 which gets rounded to 5.
        # GATK Best practices uses 75000 here: https://github.com/gatk-workflows/broad-prod-wgs-germline-snps-indels/blob/d2934ed656ade44801f9cfe1c0e78d4f80684b7b/PairedEndSingleSampleWf-fc-hg38.wdl#L778
        Int XmxGb = ceil(maxRecordsInRam / 125001.0)
        Int timeMinutes = 1 + ceil(size(inputBam, "GiB") * 3)
        String dockerImage = "quay.io/biocontainers/picard:2.26.10--hdfd78af_0"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputPath})"
        picard -Xmx~{XmxGb}G -XX:ParallelGCThreads=1 SortSam \
        INPUT=~{inputBam} \
        OUTPUT=~{outputPath} \
        MAX_RECORDS_IN_RAM=~{maxRecordsInRam} \
        SORT_ORDER=~{true="queryname" false="coordinate" sortByName} \
        CREATE_INDEX=true \
        COMPRESSION_LEVEL=~{compressionLevel} \
        USE_JDK_INFLATER=~{true="true" false="false" useJdkInflater} \
        USE_JDK_DEFLATER=~{true="true" false="false" useJdkDeflater} \
        VALIDATION_STRINGENCY=SILENT \
        CREATE_MD5_FILE=~{true="true" false="false" createMd5File}

    }

    output {
        File outputBam = outputPath
        File outputBamIndex = sub(outputPath, "\.bam$", ".bai")
    }

    runtime {
        cpu: 1
        memory: "~{1 + XmxGb}GiB"
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        inputBam: {description: "The unsorted input BAM file.", category: "required"}
        outputPath: {description: "The location the output BAM file should be written to.", category: "required"}
        sortByName: {description: "Sort the output file by name, default is position.", category: "advanced"}
        createMd5File: {description: "Whether to create an MD5 digest for any BAM or FASTQ files created.", category: "advanced"}
        maxRecordsInRam: {description: "This will specify the number of records stored in RAM before spilling to disk.", category: "advanced"}
        compressionLevel: {description: "The compression level at which the BAM files are written.", category: "advanced"}
        useJdkInflater: {description: "True, uses the java inflater. False, uses the optimized intel inflater.", category: "advanced"}
        useJdkDeflater: {description: "True, uses the java deflator to compress the BAM files. False uses the optimized intel deflater.", category: "advanced"}
        XmxGb: {description: "The maximum memory available to picard SortSam. Should be lower than `memory` to accommodate JVM overhead and BWA mem's memory usage.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        outputBam: {description: "Sorted BAM file."}
        outputBamIndex: {description: "Index of sorted BAM file."}
    }
}

task SortVcf {
    input {
        Array[File]+ vcfFiles
        String outputVcfPath

        File? dict

        String javaXmx = "8G"
        String memory = "9GiB"
        Int timeMinutes = 1 + ceil(size(vcfFiles, "GiB") * 5)
        String dockerImage = "quay.io/biocontainers/picard:2.26.10--hdfd78af_0"
    }


    command {
        set -e
        mkdir -p "$(dirname ~{outputVcfPath})"
        picard -Xmx~{javaXmx} -XX:ParallelGCThreads=1 \
        SortVcf \
        I=~{sep=" I=" vcfFiles} \
        ~{"SEQUENCE_DICTIONARY=" + dict} \
        O=~{outputVcfPath}
    }

    output {
        File outputVcf = outputVcfPath
        File outputVcfIndex = outputVcfPath + ".tbi"
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        vcfFiles: {description: "The VCF files to merge and sort.", category: "required"}
        outputVcfPath: {description: "The location the sorted VCF files should be written to.", category: "required"}
        dict: {description: "A sequence dictionary matching the VCF files.", category: "advanced"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        outputVcf: {description: "Sorted VCF file(s)."}
        outputVcfIndex: {description: "Index(es) of sort(ed) VCF file(s)."}
    }
}

task RenameSample {
    input {
        File inputVcf
        String outputPath = "./picard/renamed.vcf"
        String newSampleName

        String javaXmx = "8G"
        String memory = "9GiB"
        Int timeMinutes = 1 + ceil(size(inputVcf, "GiB") * 2)
        String dockerImage = "quay.io/biocontainers/picard:2.26.10--hdfd78af_0"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputPath})"
        picard -Xmx~{javaXmx} -XX:ParallelGCThreads=1 \
        RenameSampleInVcf \
        I=~{inputVcf} \
        O=~{outputPath} \
        NEW_SAMPLE_NAME=~{newSampleName}
    }

    output {
        File renamedVcf = outputPath
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        inputVcf: {description: "The VCF file to process.", category: "required"}
        outputPath: {description: "The location the output VCF file should be written.", category: "common"}
        newSampleName: {description: "A string to replace the old sample name.", category: "required"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.", category: "advanced"}
        memory: {description: "The memory required to run the programs.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        renamedVcf: {description: "New VCF with renamed sample."}
    }
}

task UmiAwareMarkDuplicatesWithMateCigar {
    input {
        Array[File] inputBams
        String outputPath
        String outputPathMetrics = outputPath + ".metrics"
        String outputPathUmiMetrics = outputPath + ".umi-metrics"
        Int maxRecordsInRam = 1500000  # Default is 500_000 but that will lead to very small files on disk.
        String? assumeSortOrder
        String tempdir = "temp"
        Boolean removeDuplicates = true
        String umiTagName = "RX"
        Int compressionLevel = 1
        Boolean useJdkInflater = false
        Boolean useJdkDeflater = true  # Achieves much better compression rates than the intel deflater
        String javaXmx = "8G"
        String memory = "9GiB"
        Int timeMinutes = 360
        String dockerImage = "quay.io/biocontainers/picard:2.26.10--hdfd78af_0"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputPath})" ~{tempdir}
        picard -Xmx~{javaXmx} -XX:ParallelGCThreads=1 \
        UmiAwareMarkDuplicatesWithMateCigar \
        INPUT=~{sep=' INPUT=' inputBams} \
        O=~{outputPath} \
        M=~{outputPathMetrics} \
        UMI_TAG_NAME=~{umiTagName} \
        UMI_METRICS_FILE=~{outputPathUmiMetrics} \
        TMP_DIR=~{tempdir} \
        REMOVE_DUPLICATES=~{removeDuplicates} \
        MAX_RECORDS_IN_RAM=~{maxRecordsInRam} \
        CREATE_INDEX=true \
        COMPRESSION_LEVEL=~{compressionLevel} \
        USE_JDK_INFLATER=~{true="true" false="false" useJdkInflater} \
        USE_JDK_DEFLATER=~{true="true" false="false" useJdkDeflater} \
        ~{"ASSUME_SORT_ORDER=" + assumeSortOrder}
    }

    output {
        File outputBam = outputPath
        File outputBamIndex = sub(outputPath, "\.bam$", ".bai")
        File outputMetrics = outputPathMetrics
        File outputUmiMetrics = outputPathUmiMetrics
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        inputBams: {description: "The BAM files for which the duplicate reads should be marked.", category: "required"}
        outputPath: {description: "The location the output BAM file should be written to.", category: "required"}
        outputPathMetrics: {description: "The location the output metrics file should be written to.", category: "required"}
        outputPathUmiMetrics: {description: "The location the output UMI metrics file should be written to.", category: "required"}
        removeDuplicates: {description: "Whether the duplicate reads should be removed instead of marked.", category: "common"}
        umiTagName: {description: "Which tag in the BAM file holds the UMI.", category: "common"}
        assumeSortOrder: {description: "Assume a certain sort order even though the header might say otherwise.", category: "common"}
        tempdir: {description: "Temporary directory.", category: "advanced"}
        compressionLevel: {description: "The compression level at which the BAM files are written.", category: "advanced"}
        maxRecordsInRam: {description: "This will specify the number of records stored in RAM before spilling to disk.", category: "advanced"}
        useJdkInflater: {description: "True, uses the java inflater. False, uses the optimized intel inflater.", category: "advanced"}
        useJdkDeflater: {description: "True, uses the java deflator to compress the BAM files. False uses the optimized intel deflater.", category: "advanced"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}
        
    }
}
#
# Mark duplicate reads to avoid counting non-independent observations.
task EstimateComplexity {
    input {
        File inputBam
        String metricsPath
        File referenceFasta="gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta"

        Int javaXmxMb =  6656  # 6.5G
        String memoryMb = javaXmxMb + 512
        Int? hardware_memory_GB = 16

        Int timeMinutes = 1 + ceil(size(inputBam, "GiB") * 8)
        String dockerImage = "broadinstitute/picard:2.27.5"
    }

    command {
        java "-Xmx${hardware_memory_GB}g" -jar /usr/picard/picard.jar \
        EstimateLibraryComplexity \
        -I ~{inputBam} \
        -R ~{referenceFasta} \
        -O ~{metricsPath}
    }

    output {
        File metricsFile = metricsPath
    }

    runtime {
        memory: "~{hardware_memory_GB}G"
        time_minutes: timeMinutes
        docker: dockerImage
	bootDiskSizeGb: 32
        disks: "local-disk 100 HDD"
        preemptible: 2
        maxRetries: 0
    }

    parameter_meta {
        # inputs
        inputBam: {description: "The BAM file for which the duplicate reads should be marked.", category: "required"}
        metricsPath: {description: "The location where the output metrics file should be written.", category: "required"}
        memoryMb: {description: "The amount of memory this job will use in megabytes.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        metricsFile: {description: "library complexity text file"}
    }
}
