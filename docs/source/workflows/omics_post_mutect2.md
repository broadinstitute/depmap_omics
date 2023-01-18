
## omics_post_mutect2

### Inputs

#### Required

  * `annotators` (Array[String], **required**)
  * `sample_id` (String, **required**)
  * `vcf` (File, **required**)

#### Optional

  * `open_cravat.oc_modules` (File?)
  * `open_cravat.oncokb_api_key` (File?)

#### Defaults

  * `run_open_cravat` (Boolean, default=false)
  * `RemoveFiltered.bcftools_exclude_string` (String, default='FILTER~"weak_evidence"||FILTER~"map_qual"||FILTER~"strand_bias"||FILTER~"slippage"||FILTER~"clustered_events"||FILTER~"base_qual"')
  * `RemoveFiltered.boot_disk_size` (Int, default=10)
  * `RemoveFiltered.cpu` (Int, default=2)
  * `RemoveFiltered.docker_image` (String, default="dceoy/bcftools")
  * `RemoveFiltered.mem` (Int, default=2)
  * `RemoveFiltered.preemptible` (Int, default=3)
  * `my_vcf_to_depmap.boot_disk_size` (Int, default=10)
  * `my_vcf_to_depmap.cpu` (Int, default=4)
  * `my_vcf_to_depmap.disk_space` (Int, default=40)
  * `my_vcf_to_depmap.docker_image` (String, default="python")
  * `my_vcf_to_depmap.force_keep` (Array[String], default=['oc_brca1_func_assay__class', 'oc_brca1_func_assay__score', 'clinvar_vcf_ssr', 'clinvar_vcf_clndisdbincl', 'clinvar_vcf_clnsigincl', 'clinvar_vcf_clndnincl', 'oc_oncokb__all', 'oc_oncokb__highestdiagnosticimplicationlevel', 'cgc_other_germline_mut', 'cgc_other_syndrome/disease', 'clinvar_vcf_clnsigconf', 'center', 'cosmicfusion_fusion_id', 'tumor_barcode', 'source"', 'clinvar_vcf_filter', 'clinvar_vcf_dbvarid', 'normal_barcode', 'gencode_34_ncbibuild', 'oc_oncokb_dm__highestprognosticimplicationlevel', 'strandq', 'oc_oncokb_dm__highestsensitivelevel', 'ocm', 'oc_oncokb_dm__all', 'seqq', 'nlod', 'contq', 'nalod', 'oc_base__note', 'oc_cancer_hotspots__samples', 'oc_oncokb_dm__highestresistancelevel', 'oc_oncokb_dm__tumorsummary', 'oc_oncokb_dm__highestdiagnosticimplicationlevel', 'oc_hess_drivers__signature', 'oc_hess_drivers__is_driver'])
  * `my_vcf_to_depmap.mem` (Int, default=32)
  * `my_vcf_to_depmap.n_rows` (Int, default=300000)
  * `my_vcf_to_depmap.preemptible` (Int, default=3)
  * `my_vcf_to_depmap.use_multi` (Boolean, default=false)
  * `my_vcf_to_depmap.whitelist` (Boolean, default=false)
  * `open_cravat.boot_disk_size` (Int, default=100)
  * `open_cravat.disk_space` (Int, default=20)
  * `open_cravat.docker` (String, default="karchinlab/opencravat")
  * `open_cravat.format` (String, default="vcf")
  * `open_cravat.genome` (String, default="hg38")
  * `open_cravat.memory` (Int, default=16)
  * `open_cravat.modules_options` (String, default="vcfreporter.type=separate")
  * `open_cravat.num_preempt` (Int, default=2)
  * `open_cravat.num_threads` (Int, default=4)
  * `open_cravat.retries` (Int, default=1)

### Outputs

  * `main_output` (Array[File])
  * `oc_error_files` (File?)
  * `oc_log_files` (File?)
  * `somatic_maf` (File)
