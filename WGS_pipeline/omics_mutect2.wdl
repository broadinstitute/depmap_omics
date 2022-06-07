version 1.0

# Given a set of samples, combine segment files into a single file
# more information available at https://open-cravat.readthedocs.io/en/latest/2.-Command-line-usage.html
import "gatk_mutect2_v21.wdl" as mutect2
import "bcftools.wdl" as setGT
import "fix_mutect2.wdl" as fixmutect2
import "opencravat.wdl" as openCravat
# import "filter_to_maf.wdl" as filtmaf

workflow omics_mutect2 {
  input {
    String sample_id
    String gatk_docker="broadinstitute/gatk:4.2.6.0"
    String gcs_project_for_requester_pays
    File ref_dict
    File ref_fai
    File ref_fasta

    File tumor_reads
    File tumor_reads_index

    File? intervals

    Int M2scatter=10

    File gnomad="gs://gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz"
    File gnomad_idx="gs://gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi"
    String m2_extra_args="--genotype-germline-sites true --genotype-pon-sites true"
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
      run_orientation_bias_mixture_model_filter=true
  }

  call setGT.bcftools_fix_ploidy as set_GT {
    input:
      sample_id=sample_id,
      vcf=select_first([mutect2.funcotated_file, mutect2.filtered_vcf]),
  }

  call fixmutect2.fix_mutect2 as fix_mutect2 {
      input:
        sample_id=sample_id,
        vcf_file=set_GT.vcf_fixedploid
  }

  #call openCravat.opencravat as open_cravat {
  #    input:
  #      sample_id=sample_id,
  #      vcf=fix_mutect2.vcf_fixed
  #}

  # to test
  # call filtmaf.filter_to_maf as filter_to_maf {
  #   input:
  #     sample_id=sample_id,
  #     vcf=open_cravat.oc_main_files,
  #     disk_space=20
  # }

  output {
    File out_vcf=fix_mutect2.vcf_fixed
    #File out_vcf_index=select_first([mutect2.funcotated_file_index, mutect2.filtered_vcf_idx])
    #File oc_error_files=open_cravat.oc_error_files
    #File oc_log_files=open_cravat.oc_log_files
    #File oc_sql_files=open_cravat.oc_sql_files
    #File somatic_maf=filter_to_maf.out_maf
  }
}