version 1.0
# PureCN repo, where dockerfile for image used here can be found: https://github.com/lima1/PureCN

task PureCN {
	
	input {
		# File-related inputs
		String sampleID
		File segFile
		File vcf
		File intervals

		# Method configuration inputs
		String genome="hg38"
		Int maxCopyNumber=8
		Float minPurity=0.90
		Float maxPurity=0.99
		String funSegmentation="Hclust"
		Int maxSegments=1000
		String otherArguments="--post-optimize --model-homozygous --min-total-counts 20"
		String purecn_docker="markusriester/purecn:2.2.0"

		# Hardware-related inputs
		Int hardware_disk_size_GB = 250
		Int hardware_memory_GB = 32
		Int hardware_preemptible_tries = 2
		Int num_threads = 1
		Int max_retries = 0
	}

	command {
		Rscript /opt/PureCN/PureCN.R \
			--out /cromwell_root/"${sampleID}" \
			--sampleid "${sampleID}" \
			--seg-file "${segFile}" \
			--vcf "${vcf}" \
			--intervals "${intervals}" \
			--genome "${genome}" \
			${"--max-purity " + maxPurity} \
			${"--min-purity " + minPurity} \
			${"--max-copy-number " + maxCopyNumber} \
			${"--fun-segmentation " + funSegmentation} \
			${"--max-segments " + maxSegments} \
			${otherArguments}
		Rscript -e "write.table(read.csv('${sampleID}.csv'),'table.txt',sep='\n',row.names=F,col.names=F,quote=F)"

		git clone -b dev https://github.com/broadinstitute/depmap_omics.git
		Rscript depmap_omics/WGS_pipeline/PureCN_pipeline/call_wgd_and_cin.R "${sampleID}_loh.csv" "${sampleID}.csv"
	}

	output {
		File solutions_pdf = "${sampleID}.pdf"
		File chromosomes_pdf = "${sampleID}_chromosomes.pdf"
		File rds = "${sampleID}.rds"
		File dnacopy = "${sampleID}_dnacopy.seg"
		File variants = "${sampleID}_variants.csv"
		File loh = "${sampleID}_loh.csv"
		File genes = "${sampleID}_genes.csv"
		File segmentation = "${sampleID}_segmentation.pdf"
		File log = "${sampleID}.log"
		File selected_solution = "${sampleID}.csv"
		File local_optima_pdf = "${sampleID}_local_optima.pdf"
		Array[String] table = read_lines('table.txt')
		String purity = table[1]
		String ploidy = table[2]
		String contamination = table[4]
		String flagged = table[5]
		String curated = table[7]
		String comment = table[8]
		Array[String] wgd_table = read_lines("out.txt")
		String wgd = wgd_table[0]
		String loh_fraction = wgd_table[1]
		String cin = wgd_table[2]
		String cin_allele_specific = wgd_table[3]
		String cin_ploidy_robust = wgd_table[4]
		String cin_allele_specific_ploidy_robust = wgd_table[5]
	}

	runtime {
		docker: purecn_docker
		bootDiskSizeGb: 32
		disks: "local-disk ${hardware_disk_size_GB} SSD"
		memory: "${hardware_memory_GB} GB"
		cpu: "${num_threads}"
		continueOnReturnCode: true
		preemptible: "${hardware_preemptible_tries}"
		maxRetries: "${max_retries}"
	}

}

workflow run_PureCN {

	input {
		String sampleID
		File segFile
		File vcf
		File intervals

		# Method configuration inputs
		String genome
	}

	call PureCN {
		input:
			sampleID=sampleID,
			segFile=segFile,
			vcf=vcf,
			intervals=intervals,
			genome=genome
	}

	output {
		File solutions_pdf = PureCN.solutions_pdf
		File chromosomes_pdf = PureCN.chromosomes_pdf
		File rds = PureCN.rds
		File dnacopy = PureCN.dnacopy
		File variants = PureCN.variants
		File loh = PureCN.loh
		File genes = PureCN.genes
		File segmentation = PureCN.segmentation
		File log = PureCN.log
		File selected_solution = PureCN.selected_solution
		File local_optima_pdf = PureCN.local_optima_pdf
		String purity = PureCN.purity
		String ploidy = PureCN.ploidy
		String contamination = PureCN.contamination
		String flagged = PureCN.flagged
		String curated = PureCN.curated
		String comment = PureCN.comment
		String wgd = PureCN.wgd
		String loh_fraction = PureCN.loh_fraction
		String cin = PureCN.cin
		String cin_allele_specific = PureCN.cin_allele_specific
		String cin_ploidy_robust = PureCN.cin_ploidy_robust
		String cin_allele_specific_ploidy_robust = PureCN.cin_allele_specific_ploidy_robust
	}
}
