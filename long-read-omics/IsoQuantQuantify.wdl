version 1.0

task isoquantQuantifyTask {
    input {
        String sampleName
        File inputBAM
        File inputBAMIndex
        File referenceFasta
        File ?referenceAnnotation
        Boolean ?isCompleteGeneDB
        String dataType
        String ?strandedness
        String transcriptQuantification = "with_ambiguous"
        String geneQuantification = "with_inconsistent"
        String modelConstructionStrategy
        Boolean noModelConstruction
        Int cpu = 16
        Int numThreads = 32
        Int memoryGB = 128
        Int diskSizeGB = 500
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lrtools-isoquant/lrtools-isoquant-plus@sha256:46554ca6d93d3b84d06faf493f07dd484e8300c6961ea1277483666b20e9cbcf"
        Int preemptible_tries
        File monitoringScript = "gs://mdl-refs/util/cromwell_monitoring_script2.sh"
    }

    # String file_name = basename("~{inputBAM}", ".bam")
    String model_reconstruction_arg = if noModelConstruction then "--no_model_construction" else ""

    String strandedness_present = if defined(strandedness) then select_first([strandedness]) else ""
    String stranded_arg = if defined(strandedness) then "--stranded ~{strandedness_present}" else ""

    Boolean is_complete_gene_db = if defined(isCompleteGeneDB) then select_first([isCompleteGeneDB]) else false
    String complete_gene_db_arg = if is_complete_gene_db then "--complete_genedb" else ""


    command <<<
        bash ~{monitoringScript} > monitoring.log &

        # Check if reference_annotation is provided
        ref_annotation_arg=""
        if [ -n "~{referenceAnnotation}" ]; then
            ref_annotation_arg="--genedb ~{referenceAnnotation}"
        fi



        /usr/local/src/IsoQuant-3.3.1/isoquant.py \
            --reference ~{referenceFasta} \
            ${ref_annotation_arg} ~{complete_gene_db_arg} \
            --bam ~{inputBAM} \
            --data_type ~{dataType} \
            ~{stranded_arg} \
            --transcript_quantification ~{transcriptQuantification} \
            --gene_quantification ~{geneQuantification} \
            --model_construction_strategy ~{modelConstructionStrategy} \
            --threads ~{numThreads} ~{model_reconstruction_arg} \
            --labels ~{sampleName} \
            --prefix ~{sampleName} \
            -o isoquant_output

            ls -ltrR
            echo "zipping"
            find isoquant_output/~{sampleName}/ -maxdepth 1 -type f -exec gzip {} +
            ls -ltrR
    >>>

    output {
        Array[File] isoquantOutputs = glob("isoquant_output/~{sampleName}/*.gz")
        File ?transcriptModelsGTF = "isoquant_output/~{sampleName}/~{sampleName}.transcript_models.gtf.gz"
        File ?readAssignmentsTSV = "isoquant_output/~{sampleName}/~{sampleName}.read_assignments.tsv.gz"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
        preemptible: preemptible_tries
    }
}


workflow isoquantQuantify {
    meta {
        description: "Run IsoQuant quantification (on an already gffutils preprocessed reference geneDB ideally)."
    }

    input {
        String sampleName
        File inputBAM
        File inputBAMIndex
        File referenceFasta
        File ?referenceAnnotation
        Boolean ?isCompleteGeneDB
        String dataType = "pacbio_ccs"
        String ?strandedness = "forward"
        String transcriptQuantification = "with_ambiguous"
        String geneQuantification = "with_inconsistent"
        String modelConstructionStrategy
        Boolean noModelConstruction = false
        Int preemptible_tries = 3
    }

    call isoquantQuantifyTask {
        input:
            sampleName = sampleName,
            inputBAM = inputBAM,
            inputBAMIndex = inputBAMIndex,
            referenceFasta = referenceFasta,
            referenceAnnotation = referenceAnnotation,
            isCompleteGeneDB = isCompleteGeneDB,
            dataType = dataType,
            strandedness = strandedness,
            transcriptQuantification = transcriptQuantification,
            geneQuantification = geneQuantification,
            modelConstructionStrategy = modelConstructionStrategy, 
            noModelConstruction = noModelConstruction,
            preemptible_tries = preemptible_tries
    }

    output {
        Array[File] isoquantOutputs = isoquantQuantifyTask.isoquantOutputs
        File ?transcriptModelsGTF = isoquantQuantifyTask.transcriptModelsGTF
        File ?readAssignmentsTSV = isoquantQuantifyTask.readAssignmentsTSV
        # File monitoringLog = isoquantQuantifyTask.monitoringLog
    }
}