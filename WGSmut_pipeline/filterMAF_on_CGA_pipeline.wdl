workflow filterMAF_for_CGA_workflow {
    call filterMaf
}

task filterMaf {
  File inMAFfn
  Float minAF
  Float maxExAC_AF
  Float maxExAC_AF_COSMIC
  Float maxPon_loglike
  Boolean onlyCoding
  Float minCoverage
  Float minAltReads
  File TCGAhotspotMutFN
  Float TCGAhotspotMinCnt
  Float COSMIChotspotMinCnt
  File blackListFN
  File? intersectWithMafFN
  String sample_id

  command <<<
    set -euo pipefail
    touch ${sample_id}_outMAFfn.txt
    touch ${sample_id}_outMAFannotatedFN.txt

    git clone https://github.com/broadinstitute/ccle_processing.git

    Rscript ccle_processing/mutation_pipeline/add_column.R ${inMAFfn}

    Rscript ccle_processing/mutation_pipeline/filterMAF.R CLEANED.maf ${sample_id}"_outMAFfn.txt" \
    ${sample_id}"_outMAFannotatedFN.txt" ${minAF} ${maxExAC_AF} \
    ${maxExAC_AF_COSMIC} ${maxPon_loglike} ${onlyCoding} \
    ${minCoverage} ${minAltReads} ${TCGAhotspotMutFN} ${TCGAhotspotMinCnt} \
    ${COSMIChotspotMinCnt} ${blackListFN} ${default="NULL" intersectWithMafFN}
  >>>

  output {
    File outMAFfn = "${sample_id}_outMAFfn.txt"
    File outMAFannotatedFN = "${sample_id}_outMAFannotatedFN.txt"
  }

  runtime {
    docker: "us.gcr.io/cds-docker-containers/omics-ccle_rnaseq:latest"
    memory: "10GB"
    disks: "local-disk 50 HDD"
    preemptible: "5"
  }
}