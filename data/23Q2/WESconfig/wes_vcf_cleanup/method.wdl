version 1.0

task RemoveFiltered {
	input {
		File input_vcf
		String sample_id
		String bcftools_exclude_string = 'FILTER~"weak_evidence"||FILTER~"map_qual"||FILTER~"strand_bias"||FILTER~"slippage"||FILTER~"clustered_events"||FILTER~"base_qual"'

		String docker_image="dceoy/bcftools"
		Int preemptible=3
		Int boot_disk_size=10
		Int cpu = 2
		Int mem = 2
	}

	command {
		bcftools view --exclude '~{bcftools_exclude_string}' --output-type b -o ~{sample_id}.filtered.vcf.gz --threads ~{cpu} ~{input_vcf} && \
		bcftools index ~{sample_id}.filtered.vcf.gz --threads ~{cpu}
	}

	runtime {
		disks: "local-disk 20 HDD"
		memory: "~{mem} GB"
		cpu: cpu
		preemptible: preemptible
		bootDiskSizeGb: boot_disk_size
		docker: docker_image
	}

	output {
		File output_vcf = "~{sample_id}.filtered.vcf.gz"
		File output_vcf_idx = "~{sample_id}.filtered.vcf.gz.csi"
		
	}
}

task bcftools_fix_ploidy {
	input {
		File vcf
		String sample_id
	
		Int memory = 4
		Int boot_disk_size = 10
		Int num_threads = 1
		Int num_preempt = 5
		Int disk_space = 40
		String docker = "dceoy/bcftools"
	}

	command {
		set -euo pipefail

		bcftools +setGT ${vcf} -- -t q -i'INFO/DP>8 & AF>0.9' -n c:'m|m' > ${sample_id}.vcf
		bcftools +setGT ${sample_id}.vcf -- -t q -i'AD[*:0]=0 & INFO/DP>8 & GT="0/1/2"' -n c:'1/2' > ${sample_id}.vcf.1
		bcftools +setGT ${sample_id}.vcf.1 -- -t q -i'AD[*:0]=0 & INFO/DP>8 & GT="0/1/2/3"' -n c:'1/2/3' > ${sample_id}.vcf
		bcftools +setGT ${sample_id}.vcf -- -t q -i'AD[*:0]=0 & INFO/DP>8 & GT="0/1/2/3/4"' -n c:'1/2/3/4' > ${sample_id}.vcf.1
		bcftools +setGT ${sample_id}.vcf.1 -- -t q -i'AD[*:0]=0 & INFO/DP>8 & GT="0/1/2/3/4/5"' -n c:'1/2/3/4/5' > ${sample_id}.vcf
		bcftools +setGT ${sample_id}.vcf -- -t q -i'AD[*:0]=0 & INFO/DP>8 & GT="0/1/2/3/4/5/6"' -n c:'1/2/3/4/5/6' > ${sample_id}_fixedploidy.vcf
		
		gzip ${sample_id}_fixedploidy.vcf
	}

	output {
		File vcf_fixedploid="${sample_id}_fixedploidy.vcf.gz"
	}

	runtime {
		docker: docker
		bootDiskSizeGb: "${boot_disk_size}"
		memory: "${memory}GB"
		disks: "local-disk ${disk_space} HDD"
		cpu: "${num_threads}"
		preemptible: "${num_preempt}"
	}

	meta {
		author: "Jeremie Kalfon"
	}
}

task fix_mutect2 {
	input {
		File vcf_file
		String sample_id

		Int memory = 4
		Int boot_disk_size = 10
		Int num_threads = 1
		Int num_preempt = 5
		Int disk_space = 50
		Int size_tocheck = 100
		String docker = "python"
	}

	command {
		git clone https://github.com/broadinstitute/depmap_omics.git
		pip install bgzip

		python depmap_omics/WGS_pipeline/fix_mutect2.py ${vcf_file} ${sample_id}_fixed.vcf.gz ${size_tocheck}
	}

	output {
		File vcf_fixed="${sample_id}_fixed.vcf.gz"
	}

	runtime {
		docker: docker
		bootDiskSizeGb: "${boot_disk_size}"
		memory: "${memory} GB"
		disks: "local-disk ${disk_space} HDD"
		cpu: "${num_threads}"
		preemptible: "${num_preempt}"
	}

	meta {
		author: "Jeremie Kalfon"
	}
}

workflow wes_vcf_cleanup {
	input {
		String sample_id
		File vcf
	}

	call bcftools_fix_ploidy {
		input:
		sample_id=sample_id,
		vcf=vcf,
	}

	call fix_mutect2 {
		input:
			sample_id=sample_id,
			vcf_file=bcftools_fix_ploidy.vcf_fixedploid
	}
  
	call RemoveFiltered {
		input:
			sample_id=sample_id,
			input_vcf=fix_mutect2.vcf_fixed
	}


	output {
		File out_vcf=RemoveFiltered.output_vcf
	}
}