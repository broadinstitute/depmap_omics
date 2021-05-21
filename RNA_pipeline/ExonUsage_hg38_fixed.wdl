task ExonUsage {
  String sample_id
  File exon_usage_script
  File juncReadFN
  File exonsFN
  Int minCov

  command {
    set -euo pipefail

    touch ${sample_id}_outputFN.txt
    touch ${sample_id}_outputFN.RData
    Rscript ${exon_usage_script} ${juncReadFN} ${exonsFN} ${sample_id}_outputFN.txt ${sample_id}_outputFN.RData ${minCov}
  }

  output {
    File outFN = "${sample_id}_outputFN.txt"
    File outRobjFN = "${sample_id}_outputFN.RData"
  }

  runtime {
    docker: "us-docker.pkg.dev/depmap-omics/public/ccle_docker_r"
    memory: "10GB"
    disks: "local-disk 50 HDD"
    preemptible: "5"
  }
}

workflow ExonUsage_workflow {
  call ExonUsage
}
