version 1.0

import "https://raw.githubusercontent.com/biowdl/tasks/develop/stringtie.wdl" as stringtie

workflow run_stringtie {
    input {
        File bam
        File bai
        String output_dir
    }

    call stringtie.Stringtie as Stringtie {
        input:
            bam=bam,
            bamIndex=bai,
            assembledTranscriptsFile=output_dir
    }

    output {
        File assembledTranscripts = Stringtie.assembledTranscriptsFile
        File? geneAbundance = Stringtie.geneAbundanceFile
    }
}