version 1.0

# modified "https://raw.githubusercontent.com/biowdl/tasks/develop/gridss.wdl"


struct BwaIndex {
    File fastaFile
    Array[File] indexFiles
}

workflow run_gridss2 {
    input {
        File bam
        File bai
        String sample_id

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

    call GRIDSS{
        input:
            tumorBam = [bam],
            tumorBai = [bai],
            tumorLabel = [sample_id],
            reference = ref_struct,
            dockerImage = dockerImage
    }

    output {
        File vcf = GRIDSS.vcf
        File vcfIndex = GRIDSS.vcfIndex
        File assembly = GRIDSS.assembly
        File assemblyIndex = GRIDSS.assemblyIndex
    }
}


task GRIDSS {
    input {
        Array[File]+ tumorBam
        Array[File]+ tumorBai
        Array[String]+ tumorLabel
        BwaIndex reference
        String outputPrefix = "gridss"

        File? normalBam
        File? normalBai
        String? normalLabel
        File? blacklistBed
        File? gridssProperties

        Int jvmHeapSizeGb = 64
        Int nonJvmMemoryGb = 10
        Int threads = 12
        Int disk_size = 150
        String dockerImage = "quay.io/biowdl/gridss:2.12.2"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputPrefix})"
        gridss \
        -w . \
        --reference ~{reference.fastaFile} \
        --output ~{outputPrefix}.vcf.gz \
        --assembly ~{outputPrefix}_assembly.bam \
        ~{"-c " + gridssProperties} \
        ~{"-t " + threads} \
        ~{"--jvmheap " + jvmHeapSizeGb + "G"} \
        --labels ~{normalLabel}~{true="," false="" defined(normalLabel)}~{sep="," tumorLabel} \
        ~{"--blacklist " + blacklistBed} \
        ~{normalBam} \
        ~{sep=" " tumorBam}
        samtools index ~{outputPrefix}_assembly.bam ~{outputPrefix}_assembly.bai

        # For some reason the VCF index is sometimes missing
        if [ ! -e ~{outputPrefix}.vcf.gz.tbi ]
          then
            tabix ~{outputPrefix}.vcf.gz
        fi
    }

    output {
        File vcf = outputPrefix + ".vcf.gz"
        File vcfIndex = outputPrefix + ".vcf.gz.tbi"
        File assembly = outputPrefix + "_assembly.bam"
        File assemblyIndex = outputPrefix + "_assembly.bai"
    }

    runtime {
        cpu: select_first([threads, 32])
        disks: "local-disk " + select_first([disk_size, 200]) + " HDD"
        memory: "~{jvmHeapSizeGb + nonJvmMemoryGb}GiB"
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        tumorBam: {description: "The input BAM file. This should be the tumor/case sample in case of a paired analysis.", category: "required"}
        tumorBai: {description: "The index for tumorBam.", category: "required"}
        tumorLabel: {description: "The name of the (tumor) sample.", category: "required"}
        reference: {description: "A BWA index, this should also include the fasta index file (.fai).", category: "required"}
        outputPrefix: {description: "The prefix for the output files. This may include parent directories.", category: "common"}
        normalBam: {description: "The BAM file for the normal/control sample.", category: "advanced"}
        normalBai: {description: "The index for normalBam.", category: "advanced"}
        normalLabel: {description: "The name of the normal sample.", category: "advanced"}
        blacklistBed: {description: "A bed file with blaclisted regins.", category: "advanced"}
        gridssProperties: {description: "A properties file for gridss.", category: "advanced"}

        threads: {description: "The number of the threads to use.", category: "advanced"}
        disk_size: {description: "Size of boot disk", category: "advanced"}
        jvmHeapSizeGb: {description: "The size of JVM heap for assembly and variant calling", category: "advanced"}
        nonJvmMemoryGb: {description: "The amount of memory in Gb to be requested besides JVM memory.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        vcf: {description: "VCF file including variant allele fractions."}
        vcfIndex: {description: "Index of output VCF."}
        assembly: {description: "The GRIDSS assembly BAM."}
        assemblyIndex: {description: "Index of output BAM file."}
    }
}
