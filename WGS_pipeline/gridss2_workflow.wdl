version 1.0

import "https://raw.githubusercontent.com/biowdl/tasks/develop/gridss.wdl" as gridss_tasks

struct BwaIndex {
    File fastaFile
    Array[File] indexFiles
}

workflow run_gridss2 {
    input {
        File bam
        File bai
        String sample_id

        File ref_fasta
        File ref_fasta_index
        File ref_dict
        File? ref_alt
        File ref_sa
        File ref_ann
        File ref_bwt
        File ref_pac
        File ref_amb
    }

    BwaIndex ref_struct = BwaIndex {
            fastaFile: ref_fasta,
            indexFiles: [ref_fasta_index, ref_dict, ref_alt, ref_sa, ref_ann, ref_bwt, ref_pac, ref_amb]
        }

    call gridss_tasks.GRIDSS as GRIDSS{
        input:
            tumorBam = [bam],
            tumorBai = [bai],
            tumorLabel = [sample_id],
            reference = ref_struct
    }

    output {
        File vcf = GRIDSS.vcf
        File vcfIndex = GRIDSS.vcfIndex
        File assembly = GRIDSS.assembly
        File assemblyIndex = GRIDSS.assemblyIndex
    }
}