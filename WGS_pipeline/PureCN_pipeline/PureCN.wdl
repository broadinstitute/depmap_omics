task run_PureCN {

  # File-related inputs
  String sampleID
  File segFile
  File vcf
  File intervals
  File call_wgd_and_cin_script

  # Method configuration inputs
  String genome
  Int? maxCopyNumber
  Float? minPurity
  Float? maxPurity
  String? funSegmentation
  Int? maxSegments
  String? otherArguments

  # Hardware-related inputs
  Int? hardware_disk_size_GB = 50
  Int? hardware_memory_GB = 4
  Int? hardware_preemptible_tries = 2

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
    Rscript ${call_wgd_and_cin_script} "${sampleID}_loh.csv" "${sampleID}.csv"
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
    Array[String] table = read_lines("table.txt")
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
    docker: "markusriester/purecn:latest"
    bootDiskSizeGb: 32
    disks: "local-disk ${hardware_disk_size_GB} HDD"
    memory: "${hardware_memory_GB}GB"
    cpu: 1
    continueOnReturnCode: true
    preemptible: hardware_preemptible_tries
    maxRetries: 0
  }

}

workflow PureCN {

  call run_PureCN

  output {
    File solutions_pdf = run_PureCN.solutions_pdf
    File chromosomes_pdf = run_PureCN.chromosomes_pdf
    File rds = run_PureCN.rds
    File dnacopy = run_PureCN.dnacopy
    File variants = run_PureCN.variants
    File loh = run_PureCN.loh
    File genes = run_PureCN.genes
    File segmentation = run_PureCN.segmentation
    File log = run_PureCN.log
    File selected_solution = run_PureCN.selected_solution
    File local_optima_pdf = run_PureCN.local_optima_pdf
    String purity = run_PureCN.purity
    String ploidy = run_PureCN.ploidy
    String contamination = run_PureCN.contamination
    String flagged = run_PureCN.flagged
    String curated = run_PureCN.curated
    String comment = run_PureCN.comment
    String wgd = run_PureCN.wgd
    String loh_fraction = run_PureCN.loh_fraction
    String cin = run_PureCN.cin
    String cin_allele_specific = run_PureCN.cin_allele_specific
    String cin_ploidy_robust = run_PureCN.cin_ploidy_robust
    String cin_allele_specific_ploidy_robust = run_PureCN.cin_allele_specific_ploidy_robust
    }
}