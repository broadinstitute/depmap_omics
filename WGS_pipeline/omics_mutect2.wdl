version 1.0

# more information available at https://open-cravat.readthedocs.io/en/latest/2.-Command-line-usage.html
import "mutect2_v4.2.6.1.wdl" as mutect2
import "bcftools.wdl" as setGT
import "fix_mutect2.wdl" as fixmutect2
import "merge_vcfs.wdl" as mergevcfs
import "opencravat_dm.wdl" as openCravat
import "remove_filtered.wdl" as removeFiltered
import "vcf_to_depmap.wdl" as vcf_to_depmap

workflow omics_mutect2 {
    input {
        String sample_id
        String gatk_docker="broadinstitute/gatk:4.2.6.1"
        String gcs_project_for_requester_pays= "broad-firecloud-ccle"
        File ref_dict
        File ref_fai
        File ref_fasta

        File tumor_reads
        File tumor_reads_index

        Boolean run_open_cravat=true
        Array[String] annotators=["spliceai", "alfa", "cscape", "civic", "dann", "dida", "funseq2", "gwas_catalog", "pharmgkb", "provean", "revel", "brca1_func_assay", "ccre_screen", "gtex"]
        File oncokb_api_key="gs://jkobject/oncokb_key.txt"

        File? intervals

        Int M2scatter=10

        File gnomad="gs://gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz"
        File gnomad_idx="gs://gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi"
        String m2_extra_args="--genotype-germline-sites true --genotype-pon-sites true --emit-ref-confidence GVCF"
        String? m2_filter_args
        File pon="gs://gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz"
        File pon_idx="gs://gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz.tbi"
    }

    call mutect2.Mutect2 as mutect2 {
        input:
            gatk_docker=gatk_docker,
            ref_dict=ref_dict,
            ref_fai=ref_fai,
            ref_fasta=ref_fasta,
            scatter_count=M2scatter,
            tumor_reads=tumor_reads,
            tumor_reads_index=tumor_reads_index,
            tumor_name=sample_id,
            intervals=intervals,
            gcs_project_for_requester_pays=gcs_project_for_requester_pays,
            compress_vcfs=true,
            filter_funcotations=false,
            funco_compress=true,
            funco_filter_funcotations=false,
            funco_output_format="VCF",
            funco_reference_version="hg38",
            gnomad=gnomad,
            gnomad_idx=gnomad_idx,
            m2_extra_args=m2_extra_args,
            m2_extra_filtering_args=m2_filter_args,
            make_bamout=false,
            pon=pon,
            pon_idx=pon_idx,
            run_funcotator=true,
            run_orientation_bias_mixture_model_filter=true,
    }

    call setGT.bcftools_fix_ploidy as set_GT {
        input:
            sample_id=sample_id,
            vcf=select_first([mutect2.funcotated_file, mutect2.base_vcf]),
    }

    call fixmutect2.fix_mutect2 as fix_mutect2 {
        input:
            sample_id=sample_id,
            vcf_file=set_GT.vcf_fixedploid
    }

    call removeFiltered.RemoveFiltered as RemoveFiltered {
        input:
            sample_id=sample_id,
            input_vcf=fix_mutect2.vcf_fixed
    }

    if (run_open_cravat){
        call openCravat.opencravat as open_cravat {
            input:
                vcf=RemoveFiltered.output_vcf,
                annotators_to_use=annotators,
                oncokb_api_key=oncokb_api_key
        }
    }

    call vcf_to_depmap.vcf_to_depmap as my_vcf_to_depmap {
        input:
            input_vcf=select_first([open_cravat.oc_main_file, RemoveFiltered.output_vcf]),
            sample_id=sample_id,
            annotators=flatten([annotators, ["hess_et_al"]])
    }

    output {
        Array[File] main_output=my_vcf_to_depmap.full_file
        File full_vcf=fix_mutect2.vcf_fixed
        File full_vcf_idx=select_first([mutect2.funcotated_file_index, mutect2.base_vcf_idx])
        File? oc_error_files=open_cravat.oc_error_file
        File? oc_log_files=open_cravat.oc_log_file
       # File oc_sql_files=open_cravat.oc_sql_file
        File somatic_maf=my_vcf_to_depmap.depmap_maf
    }
}