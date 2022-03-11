version 1.0

# Given a set of samples, combine segment files into a single file
# more information available at https://open-cravat.readthedocs.io/en/latest/2.-Command-line-usage.html
import "bcftools.wdl" as setGT
import "fix_mutect2col.wdl" as fixCol
import "opencravat.wdl" as openCravat
import "fix_mutect2_clust.wdl" as fixClust
# import "filter_to_maf.wdl" as filtmaf

workflow omics_post_mutect2 {
  input {
    String sample_id
    File vcf
    String annotators="spliceai alfa cscape civic mavedb uniprot loftool fitcons dann dida funseq2 genehancer gwas_catalog pharmgkb provean revel chasmplus oncokb"
  }

  call setGT.bcftools_fix_ploidy as set_GT {
    input:
      sample_id=sample_id,
      vcf=vcf,
      disk_space=50
  }

  call fixClust.fix_mutect_clust as fix_clust {
      input:
        sample_id=sample_id,
        vcf_file=set_GT.vcf_fixedploid,
        disk_space=50
  }

  call openCravat.opencravat as open_cravat {
      input:
        sample_id=sample_id,
        vcf=fix_clust.vcf_fixed,
        annotators_to_use=annotators

  }

  call fixCol.fix_column as fix_col {
    input:
      sample_id=sample_id,
      vcf_file=open_cravat.oc_main_files,
      disk_space=20
  }

  # to test
  # call filtmaf.filter_to_maf as filter_to_maf {
  #   input:
  #     sample_id=sample_id,
  #     vcf=fix_col.vcf_fixedcol,
  #     disk_space=20
  # }

  output {
    File out_vcf=fix_col.vcf_fixedcol
    File oc_error_files=open_cravat.oc_error_files
    File oc_log_files=open_cravat.oc_log_files
    File oc_sql_files=open_cravat.oc_sql_files
    #File somatic_maf=filter_to_maf.out_maf
  }
}