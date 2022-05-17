version 1.0

# more information available at https://open-cravat.readthedocs.io/en/latest/2.-Command-line-usage.html
import "merge_vcfs.wdl" as merge_vcfs
import "opencravat.wdl" as openCravat
# import "filter_to_maf.wdl" as filtmaf

workflow omics_post_mutect2 {
    input {
        String sample_id
        Array[File] vcfs
        String annotators="spliceai alfa cscape civic mavedb uniprot loftool fitcons dann dida funseq2 genehancer gwas_catalog pharmgkb provean revel chasmplus oncokb"
        File oncokb_api_key="gs://jkobject/oncokb_key.txt"
    }

    call merge_vcfs.merge_vcfs as runm_merge_vcfs {
        input:
            sample_id=sample_id,
            vcfs=vcfs,
            merge_mode='all'
    }

    call openCravat.opencravat as open_cravat {
        input:
            sample_id=sample_id,
            vcf=runm_merge_vcfs.vcf_merged,
            annotators_to_use=annotators,
            oncokb_api_key=oncokb_api_key
    }
    #call filter_annotate.filter_annotate as filter_annotate {
    #    input:
    #        vcf=open_cravat.oc_main_files,
    #}
    
  # to test
  # call filtmaf.filter_to_maf as filter_to_maf {
  #   input:
  #     sample_id=sample_id,
  #     vcf=fix_col.vcf_fixedcol,
  #     disk_space=20
  # }

    output {
        File out_vcf=open_cravat.oc_main_files
        File out_vcf_index=runm_merge_vcfs.vcf_merged_index
        File oc_sql_files=open_cravat.oc_sql_files
        #File somatic_maf=filter_to_maf.out_maf
    }
}