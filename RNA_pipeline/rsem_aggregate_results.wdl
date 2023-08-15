version 1.0

workflow rsem_aggregate_results_workflow {

	input {
		Array[File] rsem_isoforms
		Array[File] rsem_genes
		String prefix

		Int memory
		Int disk_space
		Int num_threads
		Int num_preempt
	}

	call rsem_aggregate_results {
		input:
			rsem_isoforms = rsem_isoforms,
			rsem_genes = rsem_genes,
			prefix = prefix,
			memory = memory,
			disk_space = disk_space,
			num_threads = num_threads,
			num_preempt = num_preempt,
	}

	output {
		File transcripts_tpm=rsem_aggregate_results.transcripts_tpm
		File transcripts_isopct=rsem_aggregate_results.transcripts_isopct
		File transcripts_expected_count=rsem_aggregate_results.transcripts_expected_count
		File genes_tpm=rsem_aggregate_results.genes_tpm
		File genes_expected_count=rsem_aggregate_results.genes_expected_count
		File genes_effective_length=rsem_aggregate_results.genes_effective_length
	}
}

task rsem_aggregate_results {

	input {
		Array[File] rsem_isoforms
		Array[File] rsem_genes
		String prefix

		Int memory
		Int disk_space
		Int num_threads
		Int num_preempt
	}


	command {
		git clone https://github.com/broadinstitute/depmap_omics.git 
		cd depmap_omics && git checkout add-effective-length && git pull && cd ..
		
		echo $(date +"[%b %d %H:%M:%S] Combining transcript-level output")
		python3 depmap_omics/RNA_pipeline/aggregate_rsem_results.py ${write_lines(rsem_isoforms)} TPM IsoPct expected_count ${prefix}
		echo $(date +"[%b %d %H:%M:%S] Combining gene-level output")
		python3 depmap_omics/RNA_pipeline/aggregate_rsem_results.py ${write_lines(rsem_genes)} TPM effective_length expected_count ${prefix}
	}

	output {
		File transcripts_tpm="${prefix}.rsem_transcripts_tpm.txt.gz"
		File transcripts_isopct="${prefix}.rsem_transcripts_isopct.txt.gz"
		File transcripts_expected_count="${prefix}.rsem_transcripts_expected_count.txt.gz"
		File genes_tpm="${prefix}.rsem_genes_tpm.txt.gz"
		File genes_expected_count="${prefix}.rsem_genes_expected_count.txt.gz"
		File genes_effective_length="${prefix}.rsem_genes_effective_length.txt.gz"
	}

	runtime {
		docker: "docker.io/jkobject/ccle_rnaseq:latest"
		memory: "${memory}GB"
		disks: "local-disk ${disk_space} HDD"
		cpu: "${num_threads}"
		preemptible: "${num_preempt}"
	}

	meta {
		author: "Francois Aguet, Simone Zhang"
	}
}

