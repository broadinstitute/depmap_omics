version 1.0

#import "subworkflows/annotate_variants.wdl" as VariantAnnotation
import "annotate_variants.wdl" as VariantAnnotation


workflow ctat_mutations {
    input {
        String sample_id

        # different entry points based on inputs
        File? left
        File? right
        File? bam
        File? bai
        File? vcf
        File? vcf_index

        File? extra_fasta
        Boolean merge_extra_fasta = true

        # resources - all resources derive from the ctat genome lib.
        File ref_dict
        File ref_fasta
        File ref_fasta_index
        File? gtf

        File? db_snp_vcf
        File? db_snp_vcf_index

        File? gnomad_vcf
        File? gnomad_vcf_index

        File? rna_editing_vcf
        File? rna_editing_vcf_index

        File? repeat_mask_bed

        File? ref_splice_adj_regions_bed

        File? cosmic_vcf
        File? cosmic_vcf_index

        File? ref_bed

        File? cravat_lib_tar_gz
        String? cravat_lib_dir

        String? genome_version

        File? star_reference
        String? star_reference_dir

        File? intervals
      
        Boolean annotate_variants = true


        # annotation options
        Boolean incl_snpEff = true
        Boolean incl_dbsnp = true
        Boolean incl_gnomad = true
        Boolean incl_rna_editing = true
        Boolean include_read_var_pos_annotations = true
        Boolean incl_repeats = true
        Boolean incl_homopolymers = true
        Boolean incl_splice_dist = true
        Boolean incl_blat_ED = true
        Boolean incl_cosmic = true
        Boolean incl_cravat = true


        Boolean filter_variants = true
        Boolean filter_cancer_variants = true
		
        Boolean variant_ready_bam = false
        Boolean filter_ready_vcf = false

        Boolean apply_bqsr = true
        Boolean mark_duplicates = true
        Boolean add_read_groups = true
        
        Int variant_filtration_cpu = 1
        Int variant_annotation_cpu = 5


        # boosting
        String boosting_alg_type = "classifier" #["classifier", "regressor"],
        String boosting_method = "XGBoost" #  ["none", "AdaBoost", "XGBoost", "LR", "NGBoost", "RF", "SGBoost", "SVM_RBF", "SVML"]
        Boolean seperate_snps_indels = true

        # variant attributes on which to perform boosting
        Array[String] boosting_attributes =
        ["AC","ALT","BaseQRankSum","DJ","DP","ED","Entropy","ExcessHet","FS","Homopolymer","LEN","MLEAF","MMF","QUAL","REF","RPT","RS","ReadPosRankSum","SAO","SOR","TCR","TDM","VAF","VMMF"]
        # minimum score threshold for boosted variant selection"
        Float boosting_score_threshold = 0.05

        String gatk_path = "gatk" # assume in path

        Float star_extra_disk_space = 30
        Float star_fastq_disk_space_multiplier = 10
        Boolean star_use_ssd = false
        Int star_cpu = 12
        Float star_memory = 43
        Boolean output_unmapped_reads = false

		String haplotype_caller_xtra_args = ""
        String haplotype_caller_args = "-dont-use-soft-clipped-bases --stand-call-conf 20 --recover-dangling-heads true " + haplotype_caller_xtra_args
        String haplotype_caller_args_for_extra_reads = "-dont-use-soft-clipped-bases --stand-call-conf 20 --recover-dangling-heads true " + haplotype_caller_xtra_args
        Float haplotype_caller_memory = 6.5
        String sequencing_platform = "ILLUMINA"
        Int preemptible = 2
        String docker = "trinityctat/ctat_mutations:3.0.0"
        Int variant_scatter_count = 6
        String plugins_path = "/usr/local/src/ctat-mutations/plugins"
        String scripts_path = "/usr/local/src/ctat-mutations/src"

        Boolean include_read_var_pos_annotations = true
        Float mark_duplicates_memory = 16
        Float split_n_cigar_reads_memory = 32

        Float filter_memory = 10
    }

    Boolean vcf_input = defined(vcf)

    parameter_meta {

        left:{help:"One of the two paired RNAseq samples"}
        right:{help:"One of the two paired RNAseq samples"}
        bam:{help:"Previously aligned bam file. When VCF is provided, the output from ApplyBQSR should be provided as the bam input."}
        bai:{help:"Previously aligned bam index file"}
        vcf:{help:"Previously generated vcf file to annotate and filter. When provided, the output from ApplyBQSR should be provided as the bam input."}
        sample_id:{help:"Sample id"}

        # resources
        ref_fasta:{help:"Path to the reference genome to use in the analysis pipeline."}
        ref_fasta_index:{help:"Index for ref_fasta"}
        ref_dict:{help:"Sequence dictionary for ref_fasta"}
        gtf:{help:"Annotations GTF."}

        intervals:{help:"Intervals file to restrict variant calling to. (eg. exome target list file)"}
        
        extra_fasta:{help:"Extra genome to use in alignment and variant calling."}
        merge_extra_fasta:{help:"Whether to merge extra genome fasta to use in variant calling. Set to false when extra fasta is already included in primary fasta, but you want to process the reads from extra_fasta differently."}

        db_snp_vcf:{help:"dbSNP VCF file for the reference genome."}
        db_snp_vcf_index:{help:"dbSNP vcf index"}

        gnomad_vcf:{help:"gnomad vcf file w/ allele frequencies"}
        gnomad_vcf_index:{help:"gnomad VCF index"}

        rna_editing_vcf:{help:"RNA editing VCF file"}
        rna_editing_vcf_index:{help:"RNA editing VCF index"}

        cosmic_vcf:{help:"Coding Cosmic Mutation VCF annotated with Phenotype Information"}
        cosmic_vcf_index:{help:"COSMIC VCF index"}

        repeat_mask_bed:{help:"bed file containing repetive element (repeatmasker) annotations (from UCSC genome browser download)"}

        ref_splice_adj_regions_bed:{help:"For annotating exon splice proximity"}

        ref_bed:{help:"Reference bed file for IGV cancer mutation report (refGene.sort.bed.gz)"}
        cravat_lib_tar_gz:{help:"CRAVAT resource archive"}
        cravat_lib_dir:{help:"CRAVAT resource directory (for non-Terra use)"}

        star_reference:{help:"STAR index archive"}
        star_reference_dir:{help:"STAR directory (for non-Terra use)"}

        genome_version:{help:"Genome version for annotating variants using Cravat and SnpEff", choices:["hg19", "hg38"]}

        add_read_groups : {help:"Whether to add read groups and sort the bam. Turn off for optimization with prealigned sorted bam with read groups."}
        mark_duplicates : {help:"Whether to mark duplicates"}
        filter_cancer_variants:{help:"Whether to generate cancer VCF file"}
        annotate_variants:{help:"Whether to annotate the vcf file (needed for boosting)"}
        filter_variants:{help:"Whether to filter VCF file"}
        apply_bqsr:{help:"Whether to apply base quality score recalibration"}
        #        recalibration_plot:{help:"Generate recalibration plot"}

        sequencing_platform:{help:"The sequencing platform used to generate the sample"}
        include_read_var_pos_annotations :{help: "Add vcf annotation that requires variant to be at least 6 bases from ends of reads."}

        boosting_method:{help:"Variant calling boosting method", choices:["none", "AdaBoost", "XGBoost", "LR", "NGBoost", "RF", "SGBoost", "SVM_RBF", "SVML"]}
        boosting_alg_type:{help:"Boosting algorithm type: classifier or regressor", choices:["classifier", "regressor"]}
        boosting_score_threshold:{help:"Minimum score threshold for boosted variant selection"}
        boosting_attributes:{help:"Variant attributes on which to perform boosting"}

        star_cpu:{help:"STAR aligner number of CPUs"}
        star_memory:{help:"STAR aligner memory"}
        output_unmapped_reads:{help:"Whether to output unmapped reads from STAR"}

        variant_scatter_count:{help:"Number of parallel variant caller jobs"}
        variant_filtration_cpu:{help:"Number of CPUs for variant filtration task"}
        variant_annotation_cpu:{help:"Number of CPUs for variant annotation task"}

        gatk_path:{help:"Path to GATK"}
        plugins_path:{help:"Path to plugins"}
        scripts_path:{help:"Path to scripts"}

        docker:{help:"Docker or singularity image"}
    }

    if(!vcf_input && !defined(bam) && (defined(star_reference_dir) || defined(star_reference))) {
        call StarAlign {
            input:
                star_reference = star_reference,
                star_reference_dir = star_reference_dir,
                fastq1 = left,
                fastq2 = right,
                output_unmapped_reads = output_unmapped_reads,
                genomeFastaFiles=extra_fasta,
                base_name = sample_id + '.star',
                extra_disk_space = star_extra_disk_space,
                fastq_disk_space_multiplier = star_fastq_disk_space_multiplier,
                memory = star_memory,
                use_ssd = star_use_ssd,
                cpu = star_cpu,
                docker = docker,
                preemptible = preemptible
        }
    }

    if(!vcf_input && !variant_ready_bam && add_read_groups) {
        call AddOrReplaceReadGroups {
            input:
                input_bam = select_first([StarAlign.bam, bam]),
                sample_id = sample_id,
                base_name = sample_id + '.sorted',
                sequencing_platform=sequencing_platform,
                gatk_path = gatk_path,
                docker = docker,
                preemptible = preemptible
        }
    }

    if(!vcf_input && !variant_ready_bam && mark_duplicates) {
        call MarkDuplicates {
            input:
                input_bam = select_first([AddOrReplaceReadGroups.bam, StarAlign.bam, bam]),
                base_name = sample_id + ".dedupped",
                gatk_path = gatk_path,
                memory = mark_duplicates_memory,
                docker = docker,
                preemptible = preemptible
        }
    }
    if(!vcf_input && !variant_ready_bam && defined(extra_fasta) && merge_extra_fasta) {
        call MergeFastas {
            input:
                name = "combined",
                ref_fasta = ref_fasta,
                extra_fasta = extra_fasta,
                gatk_path = gatk_path,
                docker = docker,
                preemptible = preemptible
        }
    }
    File fasta = select_first([MergeFastas.fasta, ref_fasta])
    File fasta_index = select_first([MergeFastas.fasta_index, ref_fasta_index])
    File sequence_dict = select_first([MergeFastas.sequence_dict, ref_dict])


    if( (!vcf_input) && (!variant_ready_bam) ) {
        call SplitNCigarReads {
            input:
                input_bam = select_first([MarkDuplicates.bam, AddOrReplaceReadGroups.bam, StarAlign.bam, bam]),
                input_bam_index = select_first([MarkDuplicates.bai, AddOrReplaceReadGroups.bai, StarAlign.bai, bai]),
                base_name = sample_id + ".split",
                ref_fasta = fasta,
                ref_fasta_index = fasta_index,
                ref_dict = sequence_dict,
                gatk_path = gatk_path,
                memory = split_n_cigar_reads_memory,
                docker = docker,
                preemptible = preemptible
        }

        if(apply_bqsr && defined(db_snp_vcf)) {
            call BaseRecalibrator {
                input:
                    input_bam = SplitNCigarReads.bam,
                    input_bam_index = SplitNCigarReads.bam_index,
                    recal_output_file = sample_id + ".recal_data.csv",
                    db_snp_vcf = db_snp_vcf,
                    db_snp_vcf_index = db_snp_vcf_index,
                #                known_indels_sites = known_indels_sites,
                #                known_indels_sites_indices = known_indels_sites_indices,

                    ref_fasta = fasta,
                    ref_fasta_index = fasta_index,
                    ref_dict = sequence_dict,
                    gatk_path = gatk_path,
                    docker = docker,
                    preemptible = preemptible
            }
            call ApplyBQSR {
                input:
                    input_bam = SplitNCigarReads.bam,
                    input_bam_index = SplitNCigarReads.bam_index,
                    base_name = sample_id + ".bqsr",
                    recalibration_report = BaseRecalibrator.recalibration_report,
                #                recalibration_plot = recalibration_plot,
                    ref_fasta = fasta,
                    ref_fasta_index = fasta_index,
                    ref_dict = sequence_dict,
                    gatk_path = gatk_path,
                    docker = docker,
                    preemptible = preemptible
            }
        }
    }

    if(!vcf_input && !variant_ready_bam && defined(extra_fasta)) {
        call CreateFastaIndex {
            input:
                input_fasta = extra_fasta,
                docker = docker,
                gatk_path = gatk_path,
                preemptible = preemptible
        }

        call SplitReads {
            input:
                input_bam = select_first([ApplyBQSR.bam, SplitNCigarReads.bam]),
                input_bam_index = select_first([ApplyBQSR.bam_index, SplitNCigarReads.bam_index]),
                extra_name = sample_id + '_' + basename(basename(select_first([extra_fasta]), ".fa"), ".fasta"),
                ref_name = basename(basename(ref_fasta, ".fa"), ".fasta"),
                extra_fasta_index = CreateFastaIndex.fasta_index,
                ref_fasta_index = ref_fasta_index,
                docker = docker,
                preemptible = preemptible
        }
        if(SplitReads.extra_bam_number_of_reads > 0) {
            call HaplotypeCaller as HaplotypeCallerExtra {
                input:
                    input_bam = SplitReads.extra_bam,
                    input_bam_index = SplitReads.extra_bai,
                    base_name = sample_id + '_' + basename(basename(select_first([extra_fasta]), ".fa"), ".fasta"),
                    ref_dict = select_first([CreateFastaIndex.dict]),
                    ref_fasta = select_first([CreateFastaIndex.fasta]),
                    ref_fasta_index = select_first([CreateFastaIndex.fasta_index]),
                    extra_args = haplotype_caller_args_for_extra_reads,
                    gatk_path = gatk_path,
                    docker = docker,
                    memory = haplotype_caller_memory,
                    preemptible = preemptible
            }
        }
    }


    if(!vcf_input) {
        
        File bam_for_variant_calls = select_first([SplitReads.ref_bam, ApplyBQSR.bam, SplitNCigarReads.bam, bam])
        File bai_for_variant_calls = select_first([SplitReads.ref_bai, ApplyBQSR.bam_index, SplitNCigarReads.bam_index, bai])
       
         if(variant_scatter_count > 1) {
            call SplitIntervals {
                input:
                    ref_fasta = ref_fasta,
                    ref_fasta_index = ref_fasta_index,
                    ref_dict=ref_dict,
                    intervals=intervals,
                    scatter_count = variant_scatter_count,
                    docker = docker,
                    gatk_path = gatk_path,
                    preemptible = preemptible
            }
            scatter (interval in SplitIntervals.interval_files) {
                call HaplotypeCaller as HaplotypeCallerInterval {
                    input:
                        input_bam = bam_for_variant_calls,
                        input_bam_index = bai_for_variant_calls,
                        base_name = sample_id,
                        interval_list = interval,
                        ref_dict = ref_dict,
                        ref_fasta = ref_fasta,
                        extra_args = haplotype_caller_args,
                        ref_fasta_index = ref_fasta_index,
                        gatk_path = gatk_path,
                        docker = docker,
                        memory = haplotype_caller_memory,
                        preemptible = preemptible
                }
            }
            call MergeVCFs {
                input:
                    input_vcfs = HaplotypeCallerInterval.output_vcf,
                    input_vcfs_indexes = HaplotypeCallerInterval.output_vcf_index,
                    output_vcf_name = sample_id + ".vcf.gz",
                    gatk_path = gatk_path,
                    docker = docker,
                    preemptible = preemptible
            }
        }
        if(variant_scatter_count <= 1) {
            call HaplotypeCaller {
                input:
                    input_bam = bam_for_variant_calls,
                    input_bam_index = bai_for_variant_calls,
                    base_name = sample_id,
                    ref_dict = ref_dict,
                    ref_fasta = ref_fasta,
                    extra_args = haplotype_caller_args,
                    ref_fasta_index = ref_fasta_index,
                    gatk_path = gatk_path,
                    docker = docker,
                    memory = haplotype_caller_memory,
                    preemptible = preemptible
            }
        }

        if(!vcf_input && defined(extra_fasta) && SplitReads.extra_bam_number_of_reads > 0) {
            call MergeVCFs as MergePrimaryAndExtraVCFs { # combine extra vcf with primary vcf for joint annotating and boosting
                input:
                    input_vcfs = select_all([select_first([MergeVCFs.output_vcf, HaplotypeCaller.output_vcf, vcf]), HaplotypeCallerExtra.output_vcf]),
                    input_vcfs_indexes = [],
                    output_vcf_name = sample_id + ".and.extra.vcf.gz",
                    gatk_path = gatk_path,
                    docker = docker,
                    preemptible = preemptible
            }
        }
     }

     File variant_vcf = select_first([MergePrimaryAndExtraVCFs.output_vcf, MergeVCFs.output_vcf, HaplotypeCaller.output_vcf, vcf])
     File variant_vcf_index = select_first([MergePrimaryAndExtraVCFs.output_vcf_index, MergeVCFs.output_vcf_index, HaplotypeCaller.output_vcf_index, vcf_index])

     if(annotate_variants && !filter_ready_vcf) {
        call VariantAnnotation.annotate_variants_wf as AnnotateVariants {
                input:
          input_vcf = variant_vcf,
          input_vcf_index = variant_vcf_index,
                    base_name = sample_id,
                    cravat_lib_tar_gz = cravat_lib_tar_gz,
                    cravat_lib_dir = cravat_lib_dir,
                    ref_fasta = ref_fasta,
                    ref_fasta_index = ref_fasta_index,
                    cosmic_vcf=cosmic_vcf,
                    cosmic_vcf_index=cosmic_vcf_index,
                    dbsnp_vcf=db_snp_vcf,
                    dbsnp_vcf_index=db_snp_vcf_index,
                    gnomad_vcf=gnomad_vcf,
                    gnomad_vcf_index=gnomad_vcf_index,
                    rna_editing_vcf=rna_editing_vcf,
                    rna_editing_vcf_index=rna_editing_vcf_index,
                    bam = select_first([MarkDuplicates.bam, AddOrReplaceReadGroups.bam, StarAlign.bam, bam]),
                    bam_index = select_first([MarkDuplicates.bai, AddOrReplaceReadGroups.bai, StarAlign.bai, bai]),
                    include_read_var_pos_annotations=include_read_var_pos_annotations,
                    repeat_mask_bed=repeat_mask_bed,
                    ref_splice_adj_regions_bed=ref_splice_adj_regions_bed,
                    scripts_path=scripts_path,
                    plugins_path=plugins_path,
                    genome_version=genome_version,
                    docker = docker,
                    preemptible = preemptible,
	                cpu = variant_annotation_cpu,
                    incl_snpEff = incl_snpEff,
                    incl_dbsnp = incl_dbsnp,
                    incl_gnomad = incl_gnomad,
                    incl_rna_editing = incl_rna_editing,
                    include_read_var_pos_annotations = include_read_var_pos_annotations,
                    incl_repeats = incl_repeats,
                    incl_homopolymers = incl_homopolymers,
                    incl_splice_dist = incl_splice_dist,
                    incl_blat_ED = incl_blat_ED,
                    incl_cosmic = incl_cosmic,
                    incl_cravat = incl_cravat
            }
            
      }
      
      if (filter_variants) {

            String base_name = if (boosting_method == "none") then sample_id  else "~{sample_id}.~{boosting_method}-~{boosting_alg_type}"

            call VariantFiltration {
                input:
                    input_vcf = select_first([AnnotateVariants.vcf, variant_vcf]),
                    input_vcf_index = select_first([AnnotateVariants.vcf_index, variant_vcf_index]),
                    base_name = base_name,
                    ref_dict = ref_dict,
                    ref_fasta = ref_fasta,
                    ref_fasta_index = ref_fasta_index,
                    boosting_alg_type = boosting_alg_type,
                    boosting_method =boosting_method,
                    boosting_attributes=boosting_attributes,
                    boosting_score_threshold=boosting_score_threshold,
                    gatk_path = gatk_path,
                    scripts_path=scripts_path,
                    cpu=variant_filtration_cpu,
                    memory=filter_memory,
                    docker = docker,
                    preemptible = preemptible
            }

            if(filter_cancer_variants) {
                call FilterCancerVariants {
                    input:
                        input_vcf = select_first([VariantFiltration.vcf, AnnotateVariants.vcf]),
                        base_name = base_name,
                        ref_fasta = ref_fasta,
                        ref_fasta_index = ref_fasta_index,
                        ref_dict = ref_dict,
                        scripts_path=scripts_path,
                        gatk_path = gatk_path,
                        docker = docker,
                        preemptible = preemptible
                }

                if(defined(ref_bed)) {
                    call CancerVariantReport {
                        input:
                            input_vcf = FilterCancerVariants.cancer_vcf,
                            base_name = base_name,
                            ref_fasta = ref_fasta,
                            ref_fasta_index = ref_fasta_index,
                            ref_dict = ref_dict,
                            ref_bed = select_first([ref_bed]),
                            bam=select_first([MarkDuplicates.bam, AddOrReplaceReadGroups.bam, StarAlign.bam, bam]),
                            bai=select_first([MarkDuplicates.bai, AddOrReplaceReadGroups.bai, StarAlign.bai, bai]),
                            docker = docker,
                            preemptible = preemptible
                    }
                }
            }
    }
    

    output {

        File? haplotype_caller_vcf = variant_vcf
        File? haplotype_caller_vcf_index = variant_vcf_index
        File? annotated_vcf = AnnotateVariants.vcf
        File? filtered_vcf = VariantFiltration.vcf
        File? aligned_bam = StarAlign.bam
        File? output_log_final =  StarAlign.output_log_final
        File? output_SJ =  StarAlign.output_SJ
        File? recalibrated_bam = ApplyBQSR.bam
        File? recalibrated_bam_index = ApplyBQSR.bam_index
        File? cancer_igv_report = CancerVariantReport.cancer_igv_report
        File? cancer_variants_tsv = FilterCancerVariants.cancer_variants_tsv
        File? cancer_vcf = FilterCancerVariants.cancer_vcf
    }
}

task FilterCancerVariants {
    input {
        String scripts_path
        File input_vcf
        File ref_dict

        String gatk_path
        File ref_fasta
        File ref_fasta_index

        String base_name
        String docker
        Int preemptible
    }


    command <<<
        set -e
        # monitor_script.sh &

        # Groom before table conversion
        ~{scripts_path}/groom_vcf.py \
        ~{input_vcf} ~{base_name}.cancer.groom.vcf

        ~{scripts_path}/filter_vcf_for_cancer_prediction_report.py \
        ~{base_name}.cancer.groom.vcf \
        ~{base_name}.cancer.groom.filt.vcf

        # Groom before table conversion
        ~{scripts_path}/groom_vcf.py \
        ~{base_name}.cancer.groom.filt.vcf \
        ~{base_name}.cancer.vcf

        # Convert filtered VCF file to tab file.

        ~{gatk_path} --java-options "-Xmx1500m" \
        VariantsToTable \
        -R \
        ~{ref_fasta} \
        -V \
        ~{base_name}.cancer.vcf \
        -F CHROM \
        -F POS \
        -F REF \
        -F ALT \
        -F GENE \
        -F DP \
        -F QUAL \
        -F MQ \
        -F clinvar_sig \
        -F TUMOR \
        -F TISSUE \
        -F COSMIC_ID \
        -F FATHMM \
        -F chasmplus_pval \
        -F vest_pval \
        -F mupit_link \
        --lenient \
        -O ~{base_name}.cancer.tsv

    >>>

    runtime {
        disks: "local-disk " + ceil((size(ref_fasta, "GB") * 3) + 30) + " HDD"
        docker: docker
        memory: "2G"
        preemptible: preemptible
        cpu : 1
    }
    output {
        File cancer_variants_tsv = "~{base_name}.cancer.tsv"
        File cancer_vcf = "~{base_name}.cancer.vcf"
    }
}

task CancerVariantReport {
    input {

        File input_vcf
        File bam
        File bai
        File ref_dict

        File ref_fasta
        File ref_fasta_index
        File ref_bed
        String base_name
        String docker
        Int preemptible
    }

    command <<<
        set -e
        # monitor_script.sh &

        create_report \
        ~{input_vcf} \
        ~{ref_fasta} \
        --flanking 1000 \
        --info-columns-prefixes COSMIC_ID \
        --info-columns GENE clinvar_sig FATHMM TISSUE TUMOR chasmplus_pval vest_pval mupit_link \
        --tracks ~{bam} \
        ~{ref_bed} \
        --output ~{base_name}.cancer.igvjs_viewer.html
    >>>

    runtime {
        disks: "local-disk " + ceil((size(bam, "GB") * 3) + 30) + " HDD"
        docker: docker
        memory: "2G"
        preemptible: preemptible
        cpu : 1
    }
    output {
        File cancer_igv_report = "~{base_name}.cancer.igvjs_viewer.html"
    }
}


task MarkDuplicates {
    input {
        File input_bam
        String base_name
        String gatk_path
        String docker
        Float memory
        Int preemptible
    }
    Int command_mem = ceil(memory*1000 - 500)

    command <<<
        set -e
        # monitor_script.sh &


        ~{gatk_path} --java-options "-Xmx~{command_mem}m" \
        MarkDuplicates \
        --INPUT ~{input_bam} \
        --OUTPUT ~{base_name}.bam  \
        --CREATE_INDEX true \
        --METRICS_FILE ~{base_name}.metrics
    >>>

    output {
        File bam = "${base_name}.bam"
        File bai = "${base_name}.bai"
        File metrics_file = "${base_name}.metrics"
    }

    runtime {
        disks: "local-disk " + ceil(((size(input_bam, "GB") + 2) * 10)) + " HDD"
        docker: docker
        memory: memory + "GB"
        preemptible: preemptible
    }

}

task AddOrReplaceReadGroups {
    input {
        File input_bam
        String sequencing_platform
        String base_name
        String gatk_path
        String docker
        Int preemptible
        String sample_id
    }
    String unique_id = sub(sample_id, "\\.", "_")

    command <<<
        set -e
        # monitor_script.sh &

        ~{gatk_path} --java-options "-Xmx500m" \
        AddOrReplaceReadGroups \
        --INPUT ~{input_bam} \
        --OUTPUT ~{base_name}.sorted.bam \
        --SORT_ORDER coordinate \
        --RGID id \
        --RGLB library \
        --RGPL ~{sequencing_platform} \
        --RGPU machine \
        --RGSM ~{unique_id}

        samtools index "~{base_name}.sorted.bam"
    >>>

    output {
        File bam = "~{base_name}.sorted.bam"
        File bai = "~{base_name}.sorted.bam.bai"
    }

    runtime {
        memory: "1G"
        disks: "local-disk " + ceil((size(input_bam, "GB") * 3) + 30) + " HDD"
        docker: docker
        preemptible: preemptible
    }
}

task BaseRecalibrator {
    input {
        File input_bam
        File input_bam_index
        String recal_output_file
        File? db_snp_vcf
        File? db_snp_vcf_index
        #        File known_indels_sites
        #        File known_indels_sites_indices
        File ref_dict
        File ref_fasta
        File ref_fasta_index
        String gatk_path
        String docker
        Int preemptible
    }

    output {
        File recalibration_report = recal_output_file
    }
    command <<<
        set -e
        # monitor_script.sh &


        ~{gatk_path} --java-options "-Xmx3500m" \
        BaseRecalibrator \
        -R ~{ref_fasta} \
        -I ~{input_bam} \
        --use-original-qualities \
        -O ~{recal_output_file} \
        -known-sites ~{db_snp_vcf}

    >>>
    runtime {
        memory: "4G"
        disks: "local-disk " + ceil((size(input_bam, "GB") * 3) + 30) + " HDD"
        docker: docker
        preemptible: preemptible
    }

}

task ApplyBQSR {
    input {
        File input_bam
        File input_bam_index
        String base_name
        File recalibration_report
        File ref_dict
        File ref_fasta
        File ref_fasta_index
        String gatk_path
        String docker
        Int preemptible
    }

    command <<<



        ~{gatk_path} --java-options "-Xmx3000m" \
        PrintReads \
        -I ~{input_bam} \
        -O tmp.bam

        ~{gatk_path} --java-options "-Xmx3000m" \
        ApplyBQSR \
        --add-output-sam-program-record \
        -R ~{ref_fasta} \
        -I tmp.bam \
        --use-original-qualities \
        -O ~{base_name}.bam \
        --bqsr-recal-file ~{recalibration_report}

    >>>

    output {
        File bam = "~{base_name}.bam"
        File bam_index = "~{base_name}.bai"
    }

    runtime {
        memory: "3500 MB"
        disks: "local-disk " + ceil((size(input_bam, "GB") * 4) + 30) + " HDD"
        preemptible: preemptible
        docker: docker
    }

}

task StarAlign {
    input {
        File? star_reference
        String? star_reference_dir
        File? fastq1
        File? fastq2
        Int cpu
        Float memory
        String base_name
        String docker
        Int preemptible
        Float extra_disk_space
        Float fastq_disk_space_multiplier
        Boolean use_ssd
        File? genomeFastaFiles
        Boolean output_unmapped_reads
    }
    Boolean is_gzip = sub(select_first([fastq1]), "^.+\\.(gz)$", "GZ") == "GZ"

    command <<<
        set -ex

        genomeDir="~{star_reference}"
        if [ "$genomeDir" == "" ]; then
            genomeDir="~{star_reference_dir}"
        fi

        if [ -f "${genomeDir}" ] ; then
            mkdir genome_dir
            #compress="pigz"

            #if [[ $genomeDir == *.bz2 ]] ; then
            #    compress="pbzip2"
            #fi
            #tar -I $compress -xf $genomeDir -C genome_dir --strip-components 1
            tar  -xf $genomeDir -C genome_dir --strip-components 1
            genomeDir="genome_dir"
        fi

        fastqs="~{fastq1} ~{fastq2}"
        readFilesCommand=""
        if [[ "~{fastq1}" == *.gz ]] ; then
            readFilesCommand="--readFilesCommand \"gunzip -c\""
        fi

        # special case for tar of fastq files
        if [[ "~{fastq1}" == *.tar.gz ]] ; then
            mkdir fastq
            tar -I pigz -xvf ~{fastq1} -C fastq
            fastqs=$(find fastq -type f)
            readFilesCommand=""
            if [[ "$fastqs" = *.gz ]] ; then
                readFilesCommand="--readFilesCommand \"gunzip -c\""
            fi
        fi

        STAR \
        --genomeDir $genomeDir \
        --runThreadN ~{cpu} \
        --readFilesIn $fastqs $readFilesCommand \
        --outSAMtype BAM SortedByCoordinate \
        --twopassMode Basic \
        --limitBAMsortRAM 30000000000 \
        --outSAMmapqUnique 60 \
        --outFileNamePrefix ~{base_name}. \
        ~{'--genomeFastaFiles ' + genomeFastaFiles} ~{true='--outReadsUnmapped Fastx' false='' output_unmapped_reads} \


        if [ "~{output_unmapped_reads}" == "true" ] ; then
            mv ~{base_name}.Unmapped.out.mate1 ~{base_name}.Unmapped.out.mate1.fastq
            mv ~{base_name}.Unmapped.out.mate2 ~{base_name}.Unmapped.out.mate2.fastq
        fi

        samtools index "~{base_name}.Aligned.sortedByCoord.out.bam"

    >>>

    output {
        File bam = "~{base_name}.Aligned.sortedByCoord.out.bam"
        File bai = "~{base_name}.Aligned.sortedByCoord.out.bam.bai"
        File output_log_final = "~{base_name}.Log.final.out"
        File output_log = "~{base_name}.Log.out"
        File output_SJ = "~{base_name}.SJ.out.tab"
        Array[File] unmapped_reads = glob("~{base_name}.Unmapped.out.*")
    }

    runtime {
        preemptible: preemptible
        disks: "local-disk " + ceil(size(fastq1, "GB")*fastq_disk_space_multiplier + size(fastq2, "GB") * fastq_disk_space_multiplier + size(star_reference, "GB")*8 + extra_disk_space) + " " + (if use_ssd then "SSD" else "HDD")
        docker: docker
        cpu: cpu
        memory: memory + "GB"
    }

}


task MergeVCFs {
    input {
        Array[File] input_vcfs
        Array[File] input_vcfs_indexes
        String output_vcf_name
        Int? disk_size = 5
        String gatk_path
        String docker
        Int preemptible
    }

    output {
        File output_vcf = output_vcf_name
        File output_vcf_index = "${output_vcf_name}.tbi"
    }
    command <<<
        set -e
        # monitor_script.sh &

        python <<CODE
        # make sure vcf index exists
        import subprocess
        import os
        input_vcfs = '~{sep=',' input_vcfs}'.split(',')
        for input_vcf in input_vcfs:
            if not os.path.exists(input_vcf + '.tbi') and not os.path.exists(input_vcf + '.csi') and not os.path.exists(input_vcf + '.idx'):
                subprocess.check_call(['bcftools', 'index', input_vcf])
        CODE



        ~{gatk_path} --java-options "-Xmx2000m" \
        MergeVcfs \
        -I ~{sep=" -I " input_vcfs} \
        -O ~{output_vcf_name}

    >>>
    runtime {
        memory: "2.5 GB"
        disks: "local-disk " + disk_size + " HDD"
        docker: docker
        preemptible: preemptible
    }

}

task MergeFastas {
    input {
        File ref_fasta
        File? extra_fasta
        String docker
        Int preemptible
        String name
        String gatk_path
    }


    command <<<
        cat ~{ref_fasta} ~{extra_fasta} > ~{name}.fa
        samtools faidx ~{name}.fa

        ~{gatk_path} --java-options "-Xmx1500m" \
        CreateSequenceDictionary \
        -R ~{name}.fa
    >>>

    runtime {
        docker: docker
        bootDiskSizeGb: 12
        memory: "2G"
        disks: "local-disk " + ceil(10 + 2*size(ref_fasta, "GB"))  + " HDD"
        preemptible: preemptible
        cpu: 1
    }

    output {
        File fasta = "~{name}.fa"
        File fasta_index = "~{name}.fa.fai"
        File sequence_dict = "~{name}.dict"
    }
}

task CreateFastaIndex {
    input {
        File? input_fasta
        String docker
        Int preemptible
        String gatk_path
    }
    String fasta_basename = basename(select_first([input_fasta]))
    String prefix_no_ext = basename(basename(select_first([input_fasta]), ".fa"), ".fasta")
    command <<<

        cp ~{input_fasta} ~{fasta_basename}
        samtools faidx ~{fasta_basename}


        ~{gatk_path} --java-options "-Xmx1500m" \
        CreateSequenceDictionary \
        -R ~{fasta_basename} \
        -O ~{prefix_no_ext}.dict
    >>>

    runtime {
        docker: docker
        bootDiskSizeGb: 12
        memory: "2G"
        disks: "local-disk " + ceil(size(input_fasta, "GB")*2) + " HDD"
        preemptible: preemptible
        cpu: 1
    }

    output {
        File fasta = fasta_basename
        File fasta_index = "~{fasta_basename}.fai"
        File dict = "~{prefix_no_ext}.dict"
    }
}

task SplitReads {
    input {
        File? input_bam
        File? input_bam_index
        File ref_fasta_index
        File? extra_fasta_index
        String docker
        Int preemptible
        String extra_name
        String ref_name
    }

    command <<<

        python <<CODE
        extra_chr = []
        ref_chr = []

        def parse_fai(path):
            values = set()
            with open(path, 'rt') as f:
                for line in f:
                    line = line.strip()
                    if line != '':
                        values.add(line.split('\t')[0])
            return values

        def to_txt(values, path):
            is_first = True
            with open(path, 'wt') as f:
                for val in values:
                    if not is_first:
                        f.write(' ')
                    f.write(val)
                    is_first = False

        extra_chr = parse_fai('~{extra_fasta_index}')
        ref_chr = parse_fai('~{ref_fasta_index}')
        ref_chr = ref_chr - extra_chr

        to_txt(ref_chr, 'ref.txt')
        to_txt(extra_chr, 'extra.txt')
        CODE

        samtools view -b ~{input_bam} $(cat extra.txt) > ~{extra_name}.bam
        samtools view -b ~{input_bam} $(cat ref.txt) > ~{ref_name}.bam

        samtools index ~{extra_name}.bam
        samtools index ~{ref_name}.bam

        samtools view -c -F 260 ~{extra_name}.bam > "~{extra_name}_nreads.txt"

    >>>

    runtime {
        docker: docker
        bootDiskSizeGb: 12
        memory: "2G"
        disks: "local-disk " + ceil(size(input_bam, "GB")*2 + size(extra_fasta_index, "GB")*2) + " HDD"
        preemptible: preemptible
        cpu: 1
    }

    output {
        File extra_bam = "~{extra_name}.bam"
        File extra_bai = "~{extra_name}.bam.bai"
        File ref_bam = "~{ref_name}.bam"
        File ref_bai = "~{ref_name}.bam.bai"
        Int extra_bam_number_of_reads = read_int("~{extra_name}_nreads.txt")
    }
}

task VariantFiltration {
    input {
        File input_vcf
        File input_vcf_index
        String base_name
        File ref_dict
        File ref_fasta
        File ref_fasta_index
        String boosting_alg_type
        String boosting_method
        Array[String] boosting_attributes
        Float boosting_score_threshold
        Float memory

        String scripts_path
        String gatk_path
        String docker
        Int preemptible
        Int cpu
    }


    String indel_alg_type = if (boosting_method == "LR") then "classifier" else "regressor"
    String output_name = if (boosting_method == "none") then "~{base_name}.filtered.vcf.gz" else "~{base_name}.vcf.gz"
    String boost_tmp = "~{boosting_method}_filtered.vcf"
    String ctat_boost_output_snp = "~{boosting_method}_~{boosting_alg_type}_ctat_boosting_snps.vcf.gz"
    String ctat_boost_output_indels = "~{boosting_method}_~{indel_alg_type}_ctat_boosting_indels.vcf.gz" # always regressor type for indels
    String ctat_boost_output = "~{boosting_method}_~{boosting_alg_type}_ctat_boosting.vcf"
    String median_replace_NA = if (boosting_method == "regressor") then "--replace_NA_w_median" else ""

    command <<<
        set -ex

        # monitor_script.sh &

        boosting_method="~{boosting_method}"

        if [ "$boosting_method" == "none" ]; then

            ~{gatk_path} --java-options "-Xmx2500m" \
            VariantFiltration \
            --R ~{ref_fasta} \
            --V ~{input_vcf} \
            --window 35 \
            --cluster 3 \
            --filter-name "FS" \
            --filter "FS > 30.0" \
            --filter-name "QD" \
            --filter "QD < 2.0" \
            --filter-name "SPLICEDIST" \
            --filter "DJ < 3" \
            -O tmp.vcf

            ~{gatk_path} --java-options "-Xmx2500m" \
            SelectVariants \
            --R ~{ref_fasta} \
            --V tmp.vcf \
            -select-type SNP \
            --exclude-filtered \
            -O ~{output_name}

        else

            ##############
            ## snps first:
            ~{scripts_path}/annotated_vcf_to_feature_matrix.py \
                --vcf ~{input_vcf} \
                --features ~{sep=',' boosting_attributes} \
                --snps ~{median_replace_NA} \
                --output ~{boosting_method}.snps.feature_matrix
      

            ~{scripts_path}/VariantBoosting/Apply_ML.py \
                --feature_matrix ~{boosting_method}.snps.feature_matrix \
                --snps \
                --features ~{sep=',' boosting_attributes} \
                --predictor ~{boosting_alg_type} \
                --model ~{boosting_method} \
                --output ~{boosting_method}.~{boosting_alg_type}.snps.feature_matrix.wPreds



            ##############
            ## indels next
            ~{scripts_path}/annotated_vcf_to_feature_matrix.py \
                --vcf ~{input_vcf} \
                --features ~{sep=',' boosting_attributes} \
                --indels ~{median_replace_NA} \
                --output ~{boosting_method}.indels.feature_matrix
      

            ~{scripts_path}/VariantBoosting/Apply_ML.py \
                --feature_matrix ~{boosting_method}.indels.feature_matrix \
                --indels \
                --features ~{sep=',' boosting_attributes} \
                --predictor ~{boosting_alg_type} \
                --model ~{boosting_method} \
                --output ~{boosting_method}.~{boosting_alg_type}.indels.feature_matrix.wPreds


             #########
             ## combine predictions into single output vcf
      
             ~{scripts_path}/extract_boosted_vcf.py \
                 --vcf_in ~{input_vcf} \
                 --boosted_variants_matrix ~{boosting_method}.~{boosting_alg_type}.snps.feature_matrix.wPreds ~{boosting_method}.~{boosting_alg_type}.indels.feature_matrix.wPreds\
                 --vcf_out ~{boosting_method}.~{boosting_alg_type}.vcf

             bgzip -c ~{boosting_method}.~{boosting_alg_type}.vcf > ~{output_name}
      
            
        fi
    >>>

    runtime {
        docker: docker
        cpu: cpu
        memory: memory + "GB"
        disks: "local-disk " + ceil((size(input_vcf, "GB") * 4) + (size(ref_fasta, "GB") * 2) + 30) + " HDD"
        preemptible: preemptible
    }

    output {
        File vcf = "${output_name}"
    }

}


task SplitIntervals {
    input {
        File? intervals
        File ref_fasta
        File ref_fasta_index
        File ref_dict
        Int scatter_count
        String docker
        Int preemptible
        String gatk_path
    }

    command <<<
        set -e
        # monitor_script.sh &

        mkdir interval-files
        ~{gatk_path} --java-options "-Xmx1500m" \
        SplitIntervals \
        -R ~{ref_fasta} \
        -scatter ~{scatter_count} \
        -O interval-files \
        ~{"-L " + intervals}

        cp interval-files/*.interval_list .
    >>>

    runtime {
        docker: docker
        bootDiskSizeGb: 12
        memory: "2G"
        disks: "local-disk " + ceil(size(ref_fasta, "GB")*2) + " HDD"
        preemptible: preemptible
        cpu: 1
    }

    output {
        Array[File] interval_files = glob("*.interval_list")
    }
}

task CreateBamIndex {
    input {
        File input_bam
        String docker
        Int preemptible
        String memory
    }
    String name = basename(input_bam)

    output {
        File bai = "~{name}.bai"
    }
    command <<<
        set -e

        # monitor_script.sh &

        samtools index ~{input_bam}

        mv ~{input_bam}.bai .
    >>>

    runtime {
        disks: "local-disk " + ceil(1+size(input_bam, "GB")*1.5) + " HDD"
        docker: docker
        memory: memory
        preemptible: preemptible
    }
}

task SplitNCigarReads {
    input {
        File input_bam
        File input_bam_index
        String base_name
        File ref_fasta
        File ref_fasta_index
        File ref_dict
        String gatk_path
        String docker
        Int preemptible
        Float memory
    }
    Int command_mem = ceil(memory*1000 - 500)

    output {
        File bam = "${base_name}.bam"
        File bam_index = "${base_name}.bai"
    }
    command <<<
        set -e
        # monitor_script.sh &



        ~{gatk_path} --java-options "-Xmx~{command_mem}m" \
        SplitNCigarReads \
        -R ~{ref_fasta} \
        -I ~{input_bam} \
        --read-validation-stringency LENIENT \
        -O ~{base_name}.bam
    >>>

    runtime {
        disks: "local-disk " + ceil(((size(input_bam, "GB") + 1) * 10 + size(ref_fasta, "GB") * 2)) + " HDD"
        docker: docker
        memory: memory + "GB"
        preemptible: preemptible
    }
}

task HaplotypeCaller {
    input {
        File input_bam
        File input_bam_index
        String base_name
        File? interval_list
        File ref_dict
        File ref_fasta
        File ref_fasta_index
        String gatk_path
        String docker
        Int preemptible
        Float memory
        String? extra_args
    }
    Int command_mem = ceil(memory*1000 - 500)

    output {
        File output_vcf = "${base_name}.vcf.gz"
        File output_vcf_index = "${base_name}.vcf.gz.tbi"
    }
    command <<<
        set -e
        # monitor_script.sh &


        ~{gatk_path} --java-options "-Xmx~{command_mem}m" \
        HaplotypeCaller \
        -R ~{ref_fasta} \
        -I ~{input_bam} \
        -O ~{base_name}.vcf.gz \
        ~{"" + extra_args} \
        ~{"-L " + interval_list}
    >>>
    runtime {
        docker: docker
        memory: memory + "GB"
        disks: "local-disk " + ceil((size(input_bam, "GB") * 2) + 30) + " HDD"
        preemptible: preemptible
    }
}




