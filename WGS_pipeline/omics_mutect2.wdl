# Given a set of samples, combine segment files into a single file
# more information available at https://open-cravat.readthedocs.io/en/latest/2.-Command-line-usage.html
import "https://raw.githubusercontent.com/broadinstitute/gatk/master/scripts/mutect2_wdl/mutect2.wdl" as mutect2
import "bcftools.wdl" as setGT
import "fix_mutect2col.wdl" as fixCol

workflow omics_mutect2 {
  String sample_id
  String gatk_docker="broadinstitute/gatk:4.2.4.0"

  File ref_dict
  File ref_fai
  File ref_fasta

  File tumor_reads
  File tumor_reads_index

  File bait_intervals

  String gcs_project_for_requester_pays
  Int M2cpu=4
  Int M2mem=32
  Int M2scatter=10

  File gnomad="gs://gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz"
  File gnomad_idx="gs://gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi"
  String m2_extra_args="--genotype-germline-sites true --genotype-pon-sites true"

  File pon="gs://gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz"
  File pon_idx="gs://gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz.tbi"

  call mutect2 {
    input:
      gatk_docker=gatk_docker,
      ref_dict=ref_dict,
      ref_fai=ref_fai,
      ref_fasta=ref_fasta,
      scatter_count=M2scatter,
      tumor_reads=tumor_reads,
      tumor_reads_index=tumor_reads_index,
      intervals=bait_intervals,
      cpu=M2cpu,
      mem=M2mem,
      compress_vcfs=true,
      filter_funcotations=false,
      funco_compress=true,
      funco_filter_funcotations=false,
      funco_output_format="VCF",
      funco_reference_version="hg38",
      gnomad=gnomad,
      gnomad_idx=gnomad_idx,
      m2_extra_args=m2_extra_args,
      make_bamout=false,
      pon=pon,
      pon_idx=pon_idx,
      run_funcotator=true,
      run_orientation_bias_mixture_model_filter=true
  }

  call setGT {
    input:
      sample_id=sample_id,
      vcf=mutect2.funcotated_file,
      disk_space=20
  }

  call fixCol {
    input:
      sample_id=sample_id,
      vcf=setGT.vcf_fixedploid,
      disk_space=20
  }

  output {
    File out_vcf=fixCol.vcf_fixedcol
    File out_vcf_index=mutect2.funcotated_file_index
  }
}