version 1.0

workflow TelomereWorkFlow {
  call runTelomereHunter
}

task runTelomereHunter {

  input {
      String sample_name
      File tumor_bam
      File tumor_bam_index
      File? normal_bam
      File? normal_bam_index
      File? cytoband
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
  }
}
