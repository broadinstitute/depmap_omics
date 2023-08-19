version 1.0

workflow TelomereWorkFlow {
  input {
      String sample_name
      File tumor_bam
      File tumor_bam_index
      File? normal_bam
      File? normal_bam_index
      File? cytoband
      Int preemptible=2
      Int boot_disk_size=100
      Int disk_space=200
      Int cpu = 4
      Int mem = 80
  }

  call runTelomereHunter {
      input:
          sample_name=sample_name,
          tumor_bam=tumor_bam,
          tumor_bam_index=tumor_bam_index,
          normal_bam=normal_bam,
          normal_bam_index=normal_bam_index,
          cytoband=cytoband,
          preemptible=preemptible,
          boot_disk_size=boot_disk_size,
          disk_space=disk_space,
          cpu = cpu,
          mem = mem
  }

  output {
      File TVR_top_contexts = runTelomereHunter.TVR_top_contexts
      File all_plots_merged = runTelomereHunter.all_plots_merged
      File normalized_TVR_counts = runTelomereHunter.normalized_TVR_counts
      File singletons = runTelomereHunter.singletons
      File summary = runTelomereHunter.summary
      File summary_plot = runTelomereHunter.summary_plot
      File all_plots = runTelomereHunter.all_plots
      File tumor_telomere_counts = runTelomereHunter.tumor_telomere_counts
      File? control_telomere_counts = runTelomereHunter.control_telomere_counts
  }
}

task runTelomereHunter {
  input {
      String sample_name
      File tumor_bam
      File tumor_bam_index
      File? normal_bam
      File? normal_bam_index
      File? cytoband
      Int preemptible=2
      Int boot_disk_size=60
      Int disk_space=60
      Int cpu = 4
      Int mem = 80
  }

  output {
      File TVR_top_contexts = "${sample_name}/${sample_name}_TVR_top_contexts.tsv"
      File all_plots_merged = "${sample_name}/${sample_name}_all_plots_merged.pdf"
      File normalized_TVR_counts = "${sample_name}/${sample_name}_normalized_TVR_counts.tsv"
      File singletons = "${sample_name}/${sample_name}_singletons.tsv"
      File summary = "${sample_name}/${sample_name}_summary.tsv"
      File summary_plot = "${sample_name}/${sample_name}_telomerehunter_summary_plot.pdf"
      File all_plots = "${sample_name}.telomerehunter.plots.output.tar.gz"
      File tumor_telomere_counts = "${sample_name}.telomerehunter.tumorCnt.output.tar.gz"
      File? control_telomere_counts = "${sample_name}.telomerehunter.controlCnt.output.tar.gz"
  }

  command {
    telomerehunter -ibt ${tumor_bam} \
      ~{"-ibc " + normal_bam} \
      -o . \
      -p ${sample_name} \
      ~{"-b " + cytoband}\
      --parallel

    tar -czvf ${sample_name}.telomerehunter.plots.output.tar.gz ${sample_name}/plots
    tar -czvf ${sample_name}.telomerehunter.tumorCnt.output.tar.gz ${sample_name}/tumor_TelomerCnt_${sample_name}
    if [ -d "${sample_name}/control_TelomerCnt_${sample_name}" ]; then
        tar -czvf ${sample_name}.telomerehunter.controlCnt.output.tar.gz ${sample_name}/control_TelomerCnt_${sample_name}
    fi
  }

  runtime {
      docker: "quay.io/wtsicgp/cgp-telomerehunter:1.1.0"
      disks: "local-disk ~{disk_space} HDD"
      memory: "~{mem} GB"
      cpu: cpu
      preemptible: preemptible
      bootDiskSizeGb: boot_disk_size
  }
}
