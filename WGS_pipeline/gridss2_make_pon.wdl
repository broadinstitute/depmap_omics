version 1.0

# modified "https://raw.githubusercontent.com/biowdl/tasks/develop/gridss.wdl"

import "gridss2_workflow.wdl" as gridss

struct BwaIndex {
    File fastaFile
    Array[File] indexFiles
}

workflow generate_pon_workflow {
    input {
        Array[File] bams
        Array[File] bais
        Array[String] sample_ids

        File ref_fasta
        File ref_fasta_index
        File ref_dict
        File ref_alt
        File ref_sa
        File ref_ann
        File ref_bwt
        File ref_pac
        File ref_amb

        String dockerImage = "quay.io/biowdl/gridss:2.12.2"

        BwaIndex ref_struct = object {
            fastaFile: ref_fasta,
            indexFiles: [ref_fasta_index, ref_dict, ref_alt, ref_sa, ref_ann, ref_bwt, ref_pac, ref_amb]
        }
    }

    scatter (idx in range(length(sample_ids)) ) {
        call gridss.GRIDSS as GRIDSS{
            input:
                tumorBam = [bams[idx]],
                tumorBai = [bais[idx]],
                tumorLabel = [sample_ids[idx]],
                reference = ref_struct,
                dockerImage = dockerImage,
                outputPrefix = sample_ids[idx]
        }
    }

    call GeneratePonBedpe {
        input:
            vcfFiles=GRIDSS.vcf,
            vcfIndexes=GRIDSS.vcfIndex,
            referenceFasta=ref_fasta,
            referenceFastaFai=ref_fasta_index,
    }

    output {
        File pon_bedpe = GeneratePonBedpe.bedpe
        File pon_bed = GeneratePonBedpe.bed
    }
}



task GeneratePonBedpe {
    input {
        Array[File]+ vcfFiles
        Array[File]+ vcfIndexes
        File referenceFasta
        File referenceFastaFai
        String outputDir = "."

        Int threads = 8
        String javaXmx = "8G"
        String memory = "9GiB"
        String dockerImage = "quay.io/biowdl/gridss:2.12.2"
        Int timeMinutes = 120
    }

    command {
        set -e
        mkdir -p ~{outputDir}
        java -Xmx~{javaXmx} \
        -cp /usr/local/share/gridss-2.12.2-0/gridss.jar \
        gridss.GeneratePonBedpe \
        INPUT=~{sep=" INPUT=" vcfFiles} \
        NO=0 \
        O=~{outputDir}/gridss_pon_breakpoint.bedpe \
        SBO=~{outputDir}/gridss_pon_single_breakend.bed \
        REFERENCE_SEQUENCE=~{referenceFasta} \
        THREADS=~{threads}
    }

    output {
        File bedpe = "~{outputDir}/gridss_pon_breakpoint.bedpe"
        File bed = "~{outputDir}/gridss_pon_single_breakend.bed"
    }

    runtime {
        cpu: threads
        memory: memory
        time_minutes: timeMinutes # !UnknownRuntimeKey
        docker: dockerImage
    }

    parameter_meta {
        vcfFiles: {description: "The vcf files with the normals as the first sample.", category: "required"}
        referenceFasta: {description: "The fasta of the reference genome.", category: "required"}
        referenceFastaFai: {description: "The index for the reference genome fasta.", category: "required"}
        outputDir: {description: "The directory the output will be written to.", category: "common"}
        threads: {description: "The number of the threads to use.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.",
                  category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}
