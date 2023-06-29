
## run_vcf_to_depmap

### Inputs

#### Required

  * `input_vcf` (File, **required**)
  * `sample_id` (String, **required**)

#### Defaults

  * `vcf_to_depmap.boot_disk_size` (Int, default=10)
  * `vcf_to_depmap.cpu` (Int, default=4)
  * `vcf_to_depmap.disk_space` (Int, default=40)
  * `vcf_to_depmap.docker_image` (String, default="depmapomics:test")
  * `vcf_to_depmap.force_keep` (Array[String], default=['oc_brca1_func_assay__class', 'oc_brca1_func_assay__score', 'clinvar_vcf_ssr', 'clinvar_vcf_clndisdbincl', 'clinvar_vcf_clnsigincl', 'clinvar_vcf_clndnincl', 'oc_oncokb__all', 'oc_oncokb__highestdiagnosticimplicationlevel', 'cgc_other_germline_mut', 'cgc_other_syndrome/disease', 'clinvar_vcf_clnsigconf', 'center', 'cosmicfusion_fusion_id', 'tumor_barcode', 'source"', 'clinvar_vcf_filter', 'clinvar_vcf_dbvarid', 'normal_barcode', 'gencode_34_ncbibuild', 'oc_oncokb_dm__highestprognosticimplicationlevel', 'strandq', 'oc_oncokb_dm__highestsensitivelevel', 'ocm', 'oc_oncokb_dm__all', 'seqq', 'nlod', 'contq', 'nalod', 'oc_base__note', 'oc_cancer_hotspots__samples', 'oc_oncokb_dm__highestresistancelevel', 'oc_oncokb_dm__tumorsummary', 'oc_oncokb_dm__highestdiagnosticimplicationlevel', 'oc_hess_drivers__signature', 'oc_hess_drivers__is_driver'])
  * `vcf_to_depmap.mem` (Int, default=32)
  * `vcf_to_depmap.n_rows` (Int, default=100000)
  * `vcf_to_depmap.preemptible` (Int, default=3)
  * `vcf_to_depmap.use_multi` (Boolean, default=false)
  * `vcf_to_depmap.whitelist` (Boolean, default=false)

### Outputs

  * `full_file` (Array[File])
  * `depmap_maf` (File)

## vcf_to_depmap

### Inputs

#### Required

  * `input_vcf` (File, **required**)
  * `sample_id` (String, **required**)

#### Defaults

  * `boot_disk_size` (Int, default=10)
  * `cpu` (Int, default=4)
  * `disk_space` (Int, default=40)
  * `docker_image` (String, default="depmapomics:test")
  * `force_keep` (Array[String], default=['oc_brca1_func_assay__class', 'oc_brca1_func_assay__score', 'clinvar_vcf_ssr', 'clinvar_vcf_clndisdbincl', 'clinvar_vcf_clnsigincl', 'clinvar_vcf_clndnincl', 'oc_oncokb__all', 'oc_oncokb__highestdiagnosticimplicationlevel', 'cgc_other_germline_mut', 'cgc_other_syndrome/disease', 'clinvar_vcf_clnsigconf', 'center', 'cosmicfusion_fusion_id', 'tumor_barcode', 'source"', 'clinvar_vcf_filter', 'clinvar_vcf_dbvarid', 'normal_barcode', 'gencode_34_ncbibuild', 'oc_oncokb_dm__highestprognosticimplicationlevel', 'strandq', 'oc_oncokb_dm__highestsensitivelevel', 'ocm', 'oc_oncokb_dm__all', 'seqq', 'nlod', 'contq', 'nalod', 'oc_base__note', 'oc_cancer_hotspots__samples', 'oc_oncokb_dm__highestresistancelevel', 'oc_oncokb_dm__tumorsummary', 'oc_oncokb_dm__highestdiagnosticimplicationlevel', 'oc_hess_drivers__signature', 'oc_hess_drivers__is_driver'])
  * `mem` (Int, default=32)
  * `n_rows` (Int, default=100000)
  * `preemptible` (Int, default=3)
  * `use_multi` (Boolean, default=false)
  * `whitelist` (Boolean, default=false)

### Outputs

  * `full_file` (Array[File])
  * `depmap_maf` (File)
