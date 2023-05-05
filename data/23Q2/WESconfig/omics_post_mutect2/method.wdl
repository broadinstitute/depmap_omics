version 1.0

# more information available at https://open-cravat.readthedocs.io/en/latest/2.-Command-line-usage.html
import "opencravat_dm.wdl" as https://raw.githubusercontent.com/broadinstitute/depmap_omics/dev/WGS_pipeline/opencravat_dm.wdl
import "remove_filtered.wdl" as https://raw.githubusercontent.com/broadinstitute/depmap_omics/dev/WGS_pipeline/remove_filtered.wdl
import "vcf_to_depmap.wdl" as https://raw.githubusercontent.com/broadinstitute/depmap_omics/dev/WGS_pipeline/vcf_to_depmap.wdl
import "bcftools.wdl" as https://raw.githubusercontent.com/broadinstitute/depmap_omics/dev/WGS_pipeline/bcftools.wdl
import "fix_mutect2.wdl" as https://raw.githubusercontent.com/broadinstitute/depmap_omics/dev/WGS_pipeline/fix_mutect2.wdl

workflow omics_post_mutect2 {
    input {
        String sample_id
        File vcf
        Array[String] annotators
        #File oncokb_api_key="gs://jkobject/oncokb_key.txt"
        Boolean run_open_cravat=false
    }

    #call merge_vcfs.merge_vcfs as run_merge_vcfs {
    #    input:
    #        sample_id=sample_id,
    #        vcfs=vcfs,
    #        merge_mode='all'
    #}
    
    #call bcftools.bcftools_fix_ploidy as fix_ploid {
    #    input:
    #        sample_id=sample_id,
    #        vcf=vcf,
    #}
#
    #call fix_mutect2.fix_mutect2 as fixm2 {
    #    input:
    #        sample_id=sample_id,
    #        vcf_file=fix_ploid.vcf_fixedploid
    #}
    #
    call removeFiltered.RemoveFiltered as RemoveFiltered {
        input:
            sample_id=sample_id,
            input_vcf=vcf
    }

    if (run_open_cravat){
        call openCravat.opencravat as open_cravat {
            input:
                vcf=RemoveFiltered.output_vcf,
                annotators_to_use=annotators
        }
    }

    call vcf_to_depmap.vcf_to_depmap as my_vcf_to_depmap {
        input:
            input_vcf=select_first([open_cravat.oc_main_file, RemoveFiltered.output_vcf]),
            sample_id=sample_id,
    }

    output {
        Array[File] main_output=my_vcf_to_depmap.full_file
        File? oc_error_files=open_cravat.oc_error_file
        File? oc_log_files=open_cravat.oc_log_file
        #File? oc_sql_files=open_cravat.oc_sql_files
        File somatic_maf=my_vcf_to_depmap.depmap_maf
    }
}