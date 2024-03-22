version 1.0

import "tasks/common.wdl" as common
import "tasks/picard.wdl" as picard
import "tasks/samtools.wdl" as samtools
import "tasks/mosdepth.wdl" as mosdepth

workflow BamMetrics {
    input {
        File bam
        File bamIndex
        File referenceFasta
        File referenceFastaDict
        File referenceFastaFai

        String outputDir = "."
        String strandedness = "None"
        Boolean collectAlignmentSummaryMetrics = true
        Boolean meanQualityByCycle = true

        Map[String, String] dockerImages = {
            "samtools":"quay.io/biocontainers/samtools:1.11--h6270b1f_0",
            "picard":"quay.io/biocontainers/picard:2.23.8--0",
        }
    }

    meta {
        allowNestedInputs: true
    }

    Array[Int] mapq = [20]
    String prefix = outputDir + "/" + basename(bam, ".bam")

    call samtools.Flagstat as Flagstat {
        input:
            inputBam = bam,
            outputPath = prefix + ".flagstats",
            dockerImage = dockerImages["samtools"]
    }
    scatter (q in mapq) {
        call samtools.ViewCount as Count {
            input:
                inputBam = bam,
                MAPQthreshold=q
        }
    }
    call picard.EstimateComplexity as Libcomplex {
        input:
            inputBam = bam,
            metricsPath = prefix + ".libcomplex"
    }
    call picard.CollectWgsMetrics as WgsMetrics {
        input:
            inputBam = bam,
            inputBamIndex = bamIndex,
            referenceFasta = referenceFasta,
            referenceFastaDict = referenceFastaDict,
            referenceFastaFai = referenceFastaFai
    }
    # call picard.CollectHsMetrics as CollectHsMetrics {
    #     input:
    #     	inputBam = bams[bam_and_index],
    #     	inputBamIndex = bamIndexes[bam_and_index],
    #         referenceFasta = referenceFasta,
    #         referenceFastaDict = referenceFastaDict,
    #         referenceFastaFai = referenceFastaFai,
    # targets = targetBed[bam_and_index],
    #     	basename = prefix
    # }
#Not suitable for WES due to capture of 2% genome
    #call samtools.Depth as Depth {
    #    input: inputBam = bams[bam_and_index], bamIndex = bamIndexes[bam_and_index]
    #}
    # call mosdepth.Depth as Depth2 {
    #     input:
    #     	inputBam = bams[bam_and_index],
    #     	bamIndex = bamIndexes[bam_and_index],
    # targetBed = targetBed[bam_and_index],
    #     	outputPath = prefix + "_depth"
    # }
    call picard.CollectMultipleMetrics as picardMetrics {
        input:
            inputBam = bam,
            inputBamIndex = bamIndex,
            referenceFasta = referenceFasta,
            referenceFastaDict = referenceFastaDict,
            referenceFastaFai = referenceFastaFai,
            basename = prefix,
            collectAlignmentSummaryMetrics = collectAlignmentSummaryMetrics,
            meanQualityByCycle = meanQualityByCycle,
            dockerImage = dockerImages["picard"]
    }

    output {
        File flagstats = Flagstat.flagstat
        # Array[File] depth_summary = Depth.result
        # Array[File] depth_summary2 = Depth2.result2
        # Array[File] depth_summary3 = Depth2.result_sum

        File libcomp = Libcomplex.metricsFile
        # Array[File] HsMetrics = CollectHsMetrics.HsMetrics
        Array[File] count = Count.count_stat
        Array[File] picardMetricsFiles = picardMetrics.allStats
        File wgs_metrics = WgsMetrics.metrics
    }

    parameter_meta {
        # inputs
        bam: {description: "The BAM file for which metrics will be collected.", category: "required"}
        bamIndex: {description: "The index for the bam file.", category: "required"}
        outputDir: {description: "The directory to which the outputs will be written.", category: "common"}
        referenceFasta: {description: "The reference fasta file.", category: "required"}
        referenceFastaDict: {description: "The sequence dictionary associated with the reference fasta file.", category: "required"}
        referenceFastaFai: {description: "The index for the reference fasta file.", category: "required"}
        strandedness: {description: "The strandedness of the RNA sequencing library preparation. One of \"None\" (unstranded), \"FR\" (forward-reverse: first read equal transcript) or \"RF\" (reverse-forward: second read equals transcript).", category: "common"}
        collectAlignmentSummaryMetrics: {description: "Equivalent to the `PROGRAM=CollectAlignmentSummaryMetrics` argument in Picard.", category: "advanced"}
        meanQualityByCycle: {description: "Equivalent to the `PROGRAM=MeanQualityByCycle` argument in Picard.", category: "advanced"}
        refRefflat: {description: "A refflat file containing gene annotations. If defined RNA sequencing metrics will be collected.", category: "common"}
        targetIntervals: {description: "An interval list describing the coordinates of the targets sequenced. This should only be used for targeted sequencing or WES. If defined targeted PCR metrics will be collected. Requires `ampliconIntervals` to also be defined.", category: "common"}
        ampliconIntervals: {description: "An interval list describinig the coordinates of the amplicons sequenced. This should only be used for targeted sequencing or WES. Required if `ampliconIntervals` is defined.", category: "common"}
        dockerImages: {description: "The docker images used. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        flagstats: {description: "Statistics output from flagstat."}
        reports: {description: "All reports from this pipeline gathered into one array."}
    }
}
