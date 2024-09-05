version 1.0

task star {

    input {
        File fastq1
        File? fastq2
        String prefix
        File star_index

        # STAR options
        Int outFilterMultimapNmax = 20
        Int alignSJoverhangMin = 8
        Int alignSJDBoverhangMin = 1
        Int outFilterMismatchNmax = 999
        Float outFilterMismatchNoverLmax = 0.1
        Int alignIntronMin = 20
        Int alignIntronMax = 1000000
        Int alignMatesGapMax = 1000000
        String outFilterType = "BySJout"
        Float outFilterScoreMinOverLread = 0.33
        Float outFilterMatchNminOverLread = 0.33
        Int limitSjdbInsertNsj = 1200000
        String outSAMstrandField = "intronMotif"
        String outFilterIntronMotifs = "None"
        String alignSoftClipAtReferenceEnds = "Yes"
        String quantMode = "TranscriptomeSAM GeneCounts"
        String outSAMattrRGline = "ID:rg1 SM:sm1"
        String outSAMattributes = "NH HI AS nM NM ch"
        File? varVCFfile
        String? waspOutputMode
        Int chimSegmentMin = 15
        Int chimJunctionOverhangMin = 15
        String chimOutType = "Junctions WithinBAM SoftClip"
        Int chimMainSegmentMultNmax = 1
        Int chimOutJunctionFormat = 1
        File? sjdbFileChrStartEnd

        Int memory = 64
        Int disk_space = 500
        Int num_threads = 16
        Int num_preempt = 1
    }

    command {
        set -euo pipefail

        if [[ ${fastq1} == *".tar" || ${fastq1} == *".tar.gz" ]]; then
            tar -xvvf ${fastq1}
            fastq1_abs=$(for f in *_1.fastq*; do echo "$(pwd)/$f"; done | paste -s -d ',')
            fastq2_abs=$(for f in *_2.fastq*; do echo "$(pwd)/$f"; done | paste -s -d ',')
            if [[ $fastq1_abs == *"*_1.fastq*" ]]; then  # no paired-end FASTQs found; check for single-end FASTQ
                fastq1_abs=$(for f in *.fastq*; do echo "$(pwd)/$f"; done | paste -s -d ',')
                fastq2_abs=''
            fi
        else
            # make sure paths are absolute
            fastq1_abs=${fastq1}
            fastq2_abs=${fastq2}
            if [[ $fastq1_abs != /* ]]; then
                fastq1_abs=$PWD/$fastq1_abs
                fastq2_abs=$PWD/$fastq2_abs
            fi
        fi

        echo "FASTQs:"
        echo $fastq1_abs
        echo $fastq2_abs

        # extract index
        echo $(date +"[%b %d %H:%M:%S] Extracting STAR index")
        mkdir star_index
        tar -xvvf ${star_index} -C star_index --strip-components=1

        mkdir star_out
        # placeholders for optional outputs
        touch star_out/${prefix}.Aligned.toTranscriptome.out.bam
        touch star_out/${prefix}.Chimeric.out.sorted.bam
        touch star_out/${prefix}.Chimeric.out.sorted.bam.bai
        touch star_out/${prefix}.ReadsPerGene.out.tab  # run_STAR.py will gzip

        /src/run_STAR.py \
            star_index $fastq1_abs $fastq2_abs ${prefix} \
            --output_dir star_out \
            ${"--outFilterMultimapNmax " + outFilterMultimapNmax} \
            ${"--alignSJoverhangMin " + alignSJoverhangMin} \
            ${"--alignSJDBoverhangMin " + alignSJDBoverhangMin} \
            ${"--outFilterMismatchNmax " + outFilterMismatchNmax} \
            ${"--outFilterMismatchNoverLmax " + outFilterMismatchNoverLmax} \
            ${"--alignIntronMin " + alignIntronMin} \
            ${"--alignIntronMax " + alignIntronMax} \
            ${"--alignMatesGapMax " + alignMatesGapMax} \
            ${"--outFilterType " + outFilterType} \
            ${"--outFilterScoreMinOverLread " + outFilterScoreMinOverLread} \
            ${"--outFilterMatchNminOverLread " + outFilterMatchNminOverLread} \
            ${"--limitSjdbInsertNsj " + limitSjdbInsertNsj} \
            ${"--outSAMstrandField " + outSAMstrandField} \
            ${"--outFilterIntronMotifs " + outFilterIntronMotifs} \
            ${"--alignSoftClipAtReferenceEnds " + alignSoftClipAtReferenceEnds} \
            ${"--quantMode " + quantMode} \
            ${"--outSAMattrRGline " + outSAMattrRGline} \
            ${"--outSAMattributes " + outSAMattributes} \
            ${"--varVCFfile " + varVCFfile} \
            ${"--waspOutputMode " + waspOutputMode} \
            ${"--chimSegmentMin " + chimSegmentMin} \
            ${"--chimJunctionOverhangMin " + chimJunctionOverhangMin} \
            ${"--chimOutType " + chimOutType} \
            ${"--chimMainSegmentMultNmax " + chimMainSegmentMultNmax} \
            ${"--chimOutJunctionFormat " + chimOutJunctionFormat} \
            ${"--sjdbFileChrStartEnd " + sjdbFileChrStartEnd} \
            --threads ${num_threads}
    }

    output {
        File bam_file = "star_out/${prefix}.Aligned.sortedByCoord.out.bam"
        File bam_index = "star_out/${prefix}.Aligned.sortedByCoord.out.bam.bai"
        File transcriptome_bam = "star_out/${prefix}.Aligned.toTranscriptome.out.bam"
        File chimeric_junctions = "star_out/${prefix}.Chimeric.out.junction.gz"
        File chimeric_bam_file = "star_out/${prefix}.Chimeric.out.sorted.bam"
        File chimeric_bam_index = "star_out/${prefix}.Chimeric.out.sorted.bam.bai"
        File read_counts = "star_out/${prefix}.ReadsPerGene.out.tab.gz"
        File junctions = "star_out/${prefix}.SJ.out.tab.gz"
        File junctions_pass1 = "star_out/${prefix}._STARpass1/${prefix}.SJ.pass1.out.tab.gz"
        Array[File] logs = ["star_out/${prefix}.Log.final.out", "star_out/${prefix}.Log.out", "star_out/${prefix}.Log.progress.out"]
    }

    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/gtex_rnaseq:V10"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} SSD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    meta {
        author: "Francois Aguet"
    }
}


workflow star_workflow {
    input {
        File fastq1
        File? fastq2
        String prefix
        File star_index        
    }
    call star {
        input:
            fastq1 = fastq1,
            fastq2 = fastq2,
            prefix = prefix,
            star_index = star_index, 
    }
    output {
        File bam_file = star.bam_file
        File bam_index = star.bam_index
        File transcriptome_bam = star.transcriptome_bam
        File chimeric_junctions = star.chimeric_junctions
        File chimeric_bam_file = star.chimeric_bam_file
        File chimeric_bam_index = star.chimeric_bam_index
        File read_counts = star.read_counts
        File junctions = star.junctions
        File junctions_pass1 = star.junctions_pass1
        Array[File] logs = star.logs
    }    
}
