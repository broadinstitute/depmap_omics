version 1.0
task run_update_solution {

  input {
	# File-related inputs
	File purecn_rds

	# Method configuration inputs
	Int solution_num
	String sampleID

	# Hardware-related inputs
	Int? hardware_disk_size_GB = 50
	Int? hardware_memory_GB = 1
	Int? hardware_preemptible_tries = 2
	String purecn_docker="markusriester/purecn:2.2.0"
  }

  command {
	git clone -b dev https://github.com/broadinstitute/depmap_omics.git

	Rscript depmap_omics/WGS_pipeline/PureCN_pipeline/select_solution.R ${purecn_rds} /cromwell_root/"${sampleID}" ${solution_num}
	Rscript -e "write.table(read.csv('${sampleID}.csv'),'table.txt',sep='\n',row.names=F,col.names=F,quote=F)"
	Rscript depmap_omics/WGS_pipeline/PureCN_pipeline/call_wgd_and_cin.R "${sampleID}_loh.csv" "${sampleID}.csv"
  }

  output {
	File dnacopy = "${sampleID}_dnacopy.seg"
	File variants = "${sampleID}_variants.csv"
	File loh = "${sampleID}_loh.csv"
	File genes = "${sampleID}_genes.csv"
	File selected_solution = "${sampleID}.csv"
	Array[String] table = read_lines("table.txt")
	String purity = table[1]
	String ploidy = table[2]
	String contamination = table[4]
	String flagged = table[5]
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
	memory: "${hardware_memory_GB}GB"
	cpu: 1
	continueOnReturnCode: true
	preemptible: hardware_preemptible_tries
	maxRetries: 0
  }

}

workflow update_solution {

  input {
	File purecn_rds
	Int solution_num
	String sampleID
  }

  call run_update_solution {
	  input:
	    purecn_rds = purecn_rds,
		solution_num = solution_num,
		sampleID = sampleID
  }

  output {
	File variants = run_update_solution.variants
	File loh = run_update_solution.loh
	File genes = run_update_solution.genes
	File dnacopy = run_update_solution.dnacopy
	File selected_solution = run_update_solution.selected_solution
	String purity = run_update_solution.purity
	String ploidy = run_update_solution.ploidy
	String contamination = run_update_solution.contamination
	String flagged = run_update_solution.flagged
	String comment = run_update_solution.comment
	String wgd = run_update_solution.wgd
	String loh_fraction = run_update_solution.loh_fraction
	String cin = run_update_solution.cin
	String cin_allele_specific = run_update_solution.cin_allele_specific
	String cin_ploidy_robust = run_update_solution.cin_ploidy_robust
	String cin_allele_specific_ploidy_robust = run_update_solution.cin_allele_specific_ploidy_robust
  }
}
