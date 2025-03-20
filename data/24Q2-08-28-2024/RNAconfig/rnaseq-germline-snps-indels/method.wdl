## Copyright Broad Institute, 2018
##
## Workflows for processing RNA data for germline short variant discovery with GATK (v3+v4) and related tools 
##
## Requirements/expectations :
## - BAM 
##
## Output :
## - A BAM file and its index.
## - A VCF file and its index. 
## - A Filtered VCF file and its index. 
##
## Runtime parameters are optimized for Broad's Google Cloud Platform implementation.
## For program versions, see docker containers.
##
## LICENSING :
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs. 
 
 workflow RNAseq {


	File bam_input
    File bai_input
    Boolean make_gvcf
	
	String sampleName = basename(bam_input,".bam")

	File refFasta
	File refFastaIndex
	File refDict

	String? gatk4_docker_override
	String gatk4_docker = select_first([gatk4_docker_override, "broadinstitute/gatk:latest"])
	String? gatk_path_override
	String gatk_path = select_first([gatk_path_override, "/gatk/gatk"])
	String? gitc_docker_override
	String gitc_docker = select_first([gitc_docker_override, "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"])

	File wgsCallingIntervalList

	Array[File] knownVcfs
	Array[File] knownVcfsIndices

	File dbSnpVcf
	File dbSnpVcfIndex

	Int? minConfidenceForVariantCalling

  
	## Optional user optimizations
	Int? haplotypeScatterCount
	Int scatterCount = select_first([haplotypeScatterCount, 6])
	Boolean? use_gatk4_for_all_tools
	Boolean use_all_gatk4 = select_first([use_gatk4_for_all_tools, false])


	Int? preemptible_tries
	Int preemptible_count = select_first([preemptible_tries, 3])


	if (!use_all_gatk4){
		call SplitNCigarReads {
			input:
				input_bam = bam_input,
				input_bam_index = bai_input,
				base_name = sampleName + ".split",
				ref_fasta = refFasta,
				ref_fasta_index = refFastaIndex,
				ref_dict = refDict,
				interval_list = wgsCallingIntervalList,
				preemptible_count = preemptible_count,
				docker = gitc_docker
		}
	}

	if (use_all_gatk4){
               call SplitNCigarReads_GATK4 {
                        input:
                          input_bam = bam_input,
                          input_bam_index = bai_input,
                          base_name = sampleName + ".split",
                          ref_fasta = refFasta,
                          ref_fasta_index = refFastaIndex,
                          ref_dict = refDict,
                          interval_list = wgsCallingIntervalList,
                          preemptible_count = preemptible_count,
                          docker = gatk4_docker,
                          gatk_path = gatk_path
                }
	}

	call BaseRecalibrator {
		input:
			input_bam = select_first([SplitNCigarReads_GATK4.output_bam, SplitNCigarReads.output_bam]),
			input_bam_index = select_first([SplitNCigarReads_GATK4.output_bam_index, SplitNCigarReads.output_bam_index]),
			recal_output_file = sampleName + ".recal_data.csv",
  			dbSNP_vcf = dbSnpVcf,
  			dbSNP_vcf_index = dbSnpVcfIndex,
  			known_indels_sites_VCFs = knownVcfs,
  			known_indels_sites_indices = knownVcfsIndices,
  			ref_dict = refDict,
  			ref_fasta = refFasta,
  			ref_fasta_index = refFastaIndex,
  			preemptible_count = preemptible_count,
			docker = gatk4_docker,
			gatk_path = gatk_path
	}

	call ApplyBQSR {
		input:
			input_bam =  select_first([SplitNCigarReads_GATK4.output_bam, SplitNCigarReads.output_bam]),
			input_bam_index = select_first([SplitNCigarReads_GATK4.output_bam_index, SplitNCigarReads.output_bam_index]),
			base_name = sampleName + ".aligned.duplicates_marked.recalibrated",
			ref_fasta = refFasta,
			ref_fasta_index = refFastaIndex,
			ref_dict = refDict,
			recalibration_report = BaseRecalibrator.recalibration_report,
			preemptible_count = preemptible_count,
			docker = gatk4_docker,
			gatk_path = gatk_path
	}
        if (!use_all_gatk4){
		call ScatterIntervalList {
			input:
				interval_list = wgsCallingIntervalList,
				scatter_count = scatterCount,
				preemptible_count = preemptible_count,
				docker = gitc_docker
		}
	}

        if (use_all_gatk4){
                call ScatterIntervalList_GATK4 {
                        input:
				interval_list = wgsCallingIntervalList,
				scatter_count = scatterCount,
				preemptible_count = preemptible_count,
				docker = gatk4_docker,
				gatk_path = gatk_path
                }
        }

	scatter (i in range(select_first([ScatterIntervalList_GATK4.interval_count, ScatterIntervalList.interval_count]))) {
		if (!use_all_gatk4){
			call HaplotypeCaller {
				input:
					input_bam = ApplyBQSR.output_bam,
					input_bam_index = ApplyBQSR.output_bam_index,
					base_name = sampleName + ".hc",
					interval_list = select_first([ScatterIntervalList_GATK4.out,ScatterIntervalList.out])[i],
					ref_fasta = refFasta,
					ref_fasta_index = refFastaIndex,
					ref_dict = refDict,
					dbSNP_vcf = dbSnpVcf,
                    make_gvcf= make_gvcf,
					dbSNP_vcf_index = dbSnpVcfIndex,
					stand_call_conf = minConfidenceForVariantCalling,
					preemptible_count = preemptible_count,
					docker = gitc_docker
			}
		}
		if (use_all_gatk4){
			call HaplotypeCaller_GATK4 {
        			input:
                        input_bam = ApplyBQSR.output_bam,
                        input_bam_index = ApplyBQSR.output_bam_index,
                        base_name = sampleName + ".hc",
                        interval_list = select_first([ScatterIntervalList_GATK4.out,ScatterIntervalList.out])[i],
                        ref_fasta = refFasta,
                        ref_fasta_index = refFastaIndex,
                        ref_dict = refDict,
                        make_gvcf = make_gvcf,
                        dbSNP_vcf = dbSnpVcf,
                        dbSNP_vcf_index = dbSnpVcfIndex,
                        stand_call_conf = minConfidenceForVariantCalling,
                        preemptible_count = preemptible_count,
                        docker = gatk4_docker,
						gatk_path = gatk_path
			}
		}
		File HaplotypeCallerOutputVcf = select_first([HaplotypeCaller_GATK4.output_vcf, HaplotypeCaller.output_vcf])
		File HaplotypeCallerOutputVcfIndex = select_first([HaplotypeCaller_GATK4.output_vcf_index, HaplotypeCaller.output_vcf_index])
	}

        call MergeVCFs {
                input:
                  input_vcfs = HaplotypeCallerOutputVcf,
                  input_vcfs_indexes =  HaplotypeCallerOutputVcfIndex,
                  output_vcf_name = sampleName + ".g.vcf.gz",
                  preemptible_count = preemptible_count,
                  docker = gatk4_docker,
                  gatk_path = gatk_path
        }
	
	call VariantFiltration {
		input:
			input_vcf = MergeVCFs.output_vcf,
			input_vcf_index = MergeVCFs.output_vcf_index,
			base_name = sampleName + ".variant_filtered.vcf.gz",
			ref_fasta = refFasta,
			ref_fasta_index = refFastaIndex,
			ref_dict = refDict,
			preemptible_count = preemptible_count,
			docker = gatk4_docker,
			gatk_path = gatk_path
	}

	output {
		File recalibrated_bam = ApplyBQSR.output_bam
		File recalibrated_bam_index = ApplyBQSR.output_bam_index
		File merged_vcf = MergeVCFs.output_vcf
		File merged_vcf_index = MergeVCFs.output_vcf_index
		File variant_filtered_vcf = VariantFiltration.output_vcf
		File variant_filtered_vcf_index = VariantFiltration.output_vcf_index
	}
}


## Not validated in GATK4 
task SplitNCigarReads {

	File input_bam
	File input_bam_index
	String base_name
	File interval_list

	File ref_fasta
	File ref_fasta_index
	File ref_dict

    String docker
	Int preemptible_count
    Int? memory

    command <<<
    	java -Xmx4096m -jar /usr/gitc/GATK35.jar \
    		-T SplitNCigarReads \
    		-R ${ref_fasta} \
    		-I ${input_bam} \
    		-o ${base_name}.bam \
    		-rf ReassignOneMappingQuality \
    		-RMQF 255 \
    		-RMQT 60 \
    		-U ALLOW_N_CIGAR_READS
    >>>

 	output {
 		File output_bam = "${base_name}.bam"
 		File output_bam_index = "${base_name}.bai"
 	}

    runtime {
    	disks: "local-disk " + sub(((size(input_bam,"GB")+1)*6 + size(ref_fasta,"GB")),"\\..*","") + " HDD"
		docker: docker
		memory: select_first([memory,16])+" GB"
    	preemptible: preemptible_count
    }
}

task SplitNCigarReads_GATK4 {

  File input_bam
  File input_bam_index
  String base_name
  File interval_list

  File ref_fasta
  File ref_fasta_index
  File ref_dict
  
  String gatk_path
  String docker
  Int preemptible_count
  Int? memory

    command <<<
        ${gatk_path} \
        		--java-options "-Xmx6G" \
                SplitNCigarReads \
                -R ${ref_fasta} \
                -I ${input_bam} \
                -O ${base_name}.bam 
    >>>

        output {
                File output_bam = "${base_name}.bam"
                File output_bam_index = "${base_name}.bai"
        }

    runtime {
        disks: "local-disk " + sub(((size(input_bam,"GB")+1)*6 + size(ref_fasta,"GB")),"\\..*","") + " HDD"
        docker: docker
        memory: select_first([memory,16])+" GB"
        preemptible: preemptible_count
    }
}

task BaseRecalibrator {

    File input_bam
    File input_bam_index
    String recal_output_file

    File dbSNP_vcf
    File dbSNP_vcf_index
    Array[File] known_indels_sites_VCFs
    Array[File] known_indels_sites_indices

    File ref_dict
    File ref_fasta
    File ref_fasta_index

    String gatk_path

    String docker
    Int preemptible_count
    Int? memory
    

    command <<<
        ${gatk_path} --java-options "-XX:GCTimeLimit=100 -XX:GCHeapFreeLimit=5 -XX:+PrintFlagsFinal \
            -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails \
            -Xloggc:gc_log.log -Xms8192m" \
            BaseRecalibrator \
            -R ${ref_fasta} \
            -I ${input_bam} \
            --use-original-qualities \
            -O ${recal_output_file} \
            -known-sites ${dbSNP_vcf} \
            -known-sites ${sep=" --known-sites " known_indels_sites_VCFs}
    >>>

    output {
        File recalibration_report = recal_output_file
    }

    runtime {
        memory: select_first([memory,16])+" GB"
        disks: "local-disk " + sub((size(input_bam,"GB")*3)+100, "\\..*", "") + " HDD"
        docker: docker
        preemptible: preemptible_count
    }
}


task ApplyBQSR {

    File input_bam
    File input_bam_index
    String base_name
    File recalibration_report

    File ref_dict
    File ref_fasta
    File ref_fasta_index

    String gatk_path

    String docker
    Int preemptible_count
    Int? memory

    command <<<
        ${gatk_path} \
            --java-options "-Xms8196m" \
            ApplyBQSR \
            --add-output-sam-program-record \
            -R ${ref_fasta} \
            -I ${input_bam} \
            --use-original-qualities \
            -O ${base_name}.bam \
            --bqsr-recal-file ${recalibration_report}
    >>>

    output {
        File output_bam = "${base_name}.bam"
        File output_bam_index = "${base_name}.bai"
    }

    runtime {
        memory: select_first([memory,16])+" GB"
        disks: "local-disk 300 HDD"
        preemptible: preemptible_count
        docker: docker
    }
}

task HaplotypeCaller {

	File input_bam
	File input_bam_index
	String base_name
    Boolean make_gvcf

	File interval_list

	File ref_dict
	File ref_fasta
	File ref_fasta_index

	File dbSNP_vcf
	File dbSNP_vcf_index
	String docker
	Int preemptible_count
	Int? stand_call_conf
    Int? memory

	command <<<
		java -jar /usr/gitc/GATK35.jar \
		    -T HaplotypeCaller \
		    -R ${ref_fasta} \
		    -I ${input_bam} \
		    -L ${interval_list} \
		    -dontUseSoftClippedBases \
		    -stand_call_conf ${default=20 stand_call_conf} \
		    --dbsnp ${dbSNP_vcf} \
		    -o ${base_name}.vcf.gz \
            ${true="-ERC GVCF" false="" make_gvcf}
	>>>

    output {
        File output_vcf = "${base_name}.vcf.gz"
        File output_vcf_index = "${base_name}.vcf.gz.tbi"
    }

	runtime {
		docker: docker
		memory: select_first([memory,6])+" GB"
		disks: "local-disk " + sub((size(input_bam,"GB")*2)+100, "\\..*", "") + " HDD"
		preemptible: preemptible_count
	}
}

task HaplotypeCaller_GATK4 {

	File input_bam
	File input_bam_index
	String base_name
    Boolean make_gvcf

	File interval_list

	File ref_dict
	File ref_fasta
	File ref_fasta_index

	File dbSNP_vcf
	File dbSNP_vcf_index

	String gatk_path
	String docker
	Int preemptible_count
    Int? memory

	Int? stand_call_conf

	command <<<
		${gatk_path} --java-options "-Xms8192m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
		HaplotypeCaller \
		-R ${ref_fasta} \
		-I ${input_bam} \
		-L ${interval_list} \
		-O ${base_name}.vcf.gz \
		-dont-use-soft-clipped-bases \
		--standard-min-confidence-threshold-for-calling ${default=20 stand_call_conf} \
		--dbsnp ${dbSNP_vcf} \
        ${true="-ERC GVCF" false="" make_gvcf}
	>>>

	output {
		File output_vcf = "${base_name}.vcf.gz"
		File output_vcf_index = "${base_name}.vcf.gz.tbi"
	}

	runtime {
		docker: docker
		memory: select_first([memory,10])+" GB"
		disks: "local-disk " + sub((size(input_bam,"GB")*2)+100, "\\..*", "") + " HDD"
		preemptible: preemptible_count
	}
}

task VariantFiltration {

	File input_vcf
	File input_vcf_index
	String base_name

 	File ref_dict
 	File ref_fasta
 	File ref_fasta_index

	String gatk_path
	String docker
 	Int preemptible_count
    Int? memory

	command <<<
		 ${gatk_path} \
		    VariantFiltration \
			--R ${ref_fasta} \
			--V ${input_vcf} \
			--window 35 \
			--cluster 3 \
			--filter-name "FS" \
			--filter "FS > 30.0" \
			--filter-name "QD" \
			--filter "QD < 2.0" \
			-O ${base_name}
	>>>

	output {
    	File output_vcf = "${base_name}"
    	File output_vcf_index = "${base_name}.tbi"
	}

	runtime {
		docker: docker
        bootDiskSizeGb: "32"
		memory: select_first([memory,6])+" GB"
		disks: "local-disk " + sub((size(input_vcf,"GB")*3)+180, "\\..*", "") + " HDD"
		preemptible: preemptible_count
	}
}

task MergeVCFs {
    Array[File] input_vcfs
    Array[File] input_vcfs_indexes
    String output_vcf_name

    Int? disk_size = 10
    Int? memory

    String gatk_path

    String docker
    Int preemptible_count

    # Using MergeVcfs instead of GatherVcfs so we can create indices
    # See https://github.com/broadinstitute/picard/issues/789 for relevant GatherVcfs ticket
    command <<<
        ${gatk_path} --java-options "-Xms4096m"  \
            MergeVcfs \
            --INPUT ${sep=' --INPUT=' input_vcfs} \
            --OUTPUT ${output_vcf_name}
    >>>

    output {
        File output_vcf = output_vcf_name
        File output_vcf_index = "${output_vcf_name}.tbi"
    }

    runtime {
        memory: select_first([memory,6])+" GB"
        disks: "local-disk " + disk_size + " HDD"
        docker: docker
        preemptible: preemptible_count
    }
}

task ScatterIntervalList {

    File interval_list
    Int scatter_count

    String docker
    Int preemptible_count
    Int? memory

    command <<<
        set -e
        mkdir out
        java -Xms1g -jar /usr/gitc/picard.jar \
            IntervalListTools \
            SCATTER_COUNT=${scatter_count} \
            SUBDIVISION_MODE=BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \
            UNIQUE=true \
            SORT=true \
            INPUT=${interval_list} \
            OUTPUT=out

        python3 <<CODE
        import glob, os
        # Works around a JES limitation where multiples files with the same name overwrite each other when globbed
        intervals = sorted(glob.glob("out/*/*.interval_list"))
        for i, interval in enumerate(intervals):
          (directory, filename) = os.path.split(interval)
          newName = os.path.join(directory, str(i + 1) + filename)
          os.rename(interval, newName)
        print(len(intervals))
        CODE
    >>>

    output {
        Array[File] out = glob("out/*/*.interval_list")
        Int interval_count = read_int(stdout())
    }

    runtime {
        disks: "local-disk 1 HDD"
        memory: select_first([memory,4])+" GB"
        docker: docker
        preemptible: preemptible_count
    }
}

task ScatterIntervalList_GATK4 {

	File interval_list
	Int scatter_count
	String gatk_path
	String docker
	Int preemptible_count
    Int? memory

    command <<<
        set -e
        mkdir out
        ${gatk_path} --java-options "-Xms2g" \
            IntervalListTools \
            --SCATTER_COUNT=${scatter_count} \
            --SUBDIVISION_MODE=BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \
            --UNIQUE=true \
            --SORT=true \
            --INPUT=${interval_list} \
            --OUTPUT=out
	
        python3 <<CODE
        import glob, os
        # Works around a JES limitation where multiples files with the same name overwrite each other when globbed
        intervals = sorted(glob.glob("out/*/*.interval_list"))
        for i, interval in enumerate(intervals):
          (directory, filename) = os.path.split(interval)
          newName = os.path.join(directory, str(i + 1) + filename)
          os.rename(interval, newName)
        print(len(intervals))
        f = open("interval_count.txt", "w+")
        f.write(str(len(intervals)))
        f.close()
        CODE
    >>>

    output {
        Array[File] out = glob("out/*/*.interval_list")
        Int interval_count = read_int("interval_count.txt")
    }

    runtime {
        disks: "local-disk 1 HDD"
        memory: select_first([memory,4])+" GB"
        docker: docker
        preemptible: preemptible_count
    }
}

task RevertSam {
    File input_bam
    String base_name
    String sort_order

    String gatk_path

    String docker
    Int preemptible_count

    command <<<
        ${gatk_path} \
        	RevertSam \
        	--INPUT ${input_bam} \
        	--OUTPUT ${base_name}.bam \
            --VALIDATION_STRINGENCY SILENT \
        	--ATTRIBUTE_TO_CLEAR FT \
        	--ATTRIBUTE_TO_CLEAR CO \
        	--SORT_ORDER ${sort_order}
    >>>

    output {
        File output_bam = "${base_name}.bam"
    }

    runtime {
        docker: docker
        disks: "local-disk " + sub(((size(input_bam,"GB")+1)*5),"\\..*","") + " HDD"
        memory: "8 GB"
        preemptible: preemptible_count
    }
}
