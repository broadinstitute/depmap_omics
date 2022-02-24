task star {

    File fastq1
    File? fastq2
    String prefix
    File star_index

    # STAR options
    Int? outFilterMultimapNmax
    Int? alignSJoverhangMin
    Int? outFilterMismatchNmax
    Float? outFilterMismatchNoverLmax
    Int? alignIntronMin
    String? outFilterType
    Float? outFilterScoreMinOverLread
    Float? outFilterMatchNminOverLread
    Int? limitSjdbInsertNsj
    String? outFilterIntronMotifs
    String? alignSoftClipAtReferenceEnds
    String? quantMode
    String? outSAMattributes
    File? varVCFfile
    String? waspOutputMode
    String? chimOutType
    Int? chimMainSegmentMultNmax
    File? sjdbFileChrStartEnd


    # parameters for STAR as part of STAR-fusion
    String outReadsUnmapped="None"
    String outSAMstrandField="IntronMotif"
    String outSAMunmapped="Within"
    Int chimSegmentMin=12
    Int chimJunctionOverhangMin=8
    Int chimOutJunctionFormat=1
    Int alignSJDBoverhangMin=10
    Int alignMatesGapMax=100000
    Int alignIntronMax=100000
    String lignSJstitchMismatchNmax="5 -1 5 5"
    String outSAMattrRGline="ID:GRPundef"
    Int chimMultimapScoreRange=3
    Int chimScoreJunctionNonGTAG=-4
    Int chimMultimapNmax=20
    Int chimNonchimScoreDropMin=10
    Int peOverlapNbasesMin=12
    Float peOverlapMMp=0.1
    String alignInsertionFlush="Right"
    Int alignSplicedMateMapLminOverLmate=0
    Int alignSplicedMateMapLmin=30



    Int memory
    Int disk_space
    Int num_threads
    Int num_preempt

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

        git clone -b test_star --single-branch https://github.com/broadinstitute/depmap_omics.git

        depmap_omics/RNA_pipeline/run_STAR_for_fusion.py \
            star_index $fastq1_abs $fastq2_abs ${prefix} \
            --output_dir star_out \
            ${"--outFilterMultimapNmax " + outFilterMultimapNmax} \
            ${"--alignSJoverhangMin " + alignSJoverhangMin} \
            ${"--outFilterMismatchNmax " + outFilterMismatchNmax} \
            ${"--outFilterMismatchNoverLmax " + outFilterMismatchNoverLmax} \
            ${"--alignIntronMin " + alignIntronMin} \
            ${"--outFilterType " + outFilterType} \
            ${"--outFilterScoreMinOverLread " + outFilterScoreMinOverLread} \
            ${"--outFilterMatchNminOverLread " + outFilterMatchNminOverLread} \
            ${"--limitSjdbInsertNsj " + limitSjdbInsertNsj} \
            ${"--outFilterIntronMotifs " + outFilterIntronMotifs} \
            ${"--alignSoftClipAtReferenceEnds " + alignSoftClipAtReferenceEnds} \
            ${"--quantMode " + quantMode} \
            ${"--outSAMattributes " + outSAMattributes} \
            ${"--varVCFfile " + varVCFfile} \
            ${"--waspOutputMode " + waspOutputMode} \
            ${"--chimOutType " + chimOutType} \
            ${"--chimMainSegmentMultNmax " + chimMainSegmentMultNmax} \
            ${"--sjdbFileChrStartEnd " + sjdbFileChrStartEnd} \
            ${"--outReadsUnmapped " + outReadsUnmapped} \
            ${"--outSAMstrandField " + outSAMstrandField} \  # include for potential use with StringTie for assembly
            ${"--outSAMunmapped " + outSAMunmapped} \
            ${"--chimSegmentMin " + chimSegmentMin} \  # ** essential to invoke chimeric read detection & reporting **
            ${"--chimJunctionOverhangMin " + chimJunctionOverhangMin} \
            ${"--chimOutJunctionFormat " + chimOutJunctionFormat} \   # **essential** includes required metadata in Chimeric.junction.out file.
            ${"--alignSJDBoverhangMin " + alignSJDBoverhangMin} \
            ${"--alignMatesGapMax " + alignMatesGapMax} \   # avoid readthru fusions within 100k
            ${"--alignIntronMax " + alignIntronMax} \
            ${"--alignSJstitchMismatchNmax " + alignSJstitchMismatchNmax} \   # settings improved certain chimera detections
            ${"--outSAMattrRGline " + outSAMattrRGline} \
            ${"--chimMultimapScoreRange " + chimMultimapScoreRange} \
            ${"--chimScoreJunctionNonGTAG " + chimScoreJunctionNonGTAG} \
            ${"--chimMultimapNmax " + chimMultimapNmax} \
            ${"--chimNonchimScoreDropMin " + chimNonchimScoreDropMin} \
            ${"--peOverlapNbasesMin " + peOverlapNbasesMin} \
            ${"--peOverlapMMp " + peOverlapMMp} \
            ${"--alignInsertionFlush " + alignInsertionFlush} \
            ${"--alignSplicedMateMapLminOverLmate " + alignSplicedMateMapLminOverLmate} \
            ${"--alignSplicedMateMapLmin " + alignSplicedMateMapLmin} \
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
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    meta {
        author: "Simone Zhang"
    }
}


workflow star_workflow {
    call star
}