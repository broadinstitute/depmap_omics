task rmats_turbo {

    File bam
    File star_index
    String prefix

    Int memory
    Int disk_space
    Int num_threads
    Int num_preempt

    command {
        set -euo pipefail
        
        wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.primary_assembly.annotation.gtf.gz
        gunzip gencode.v38.primary_assembly.annotation.gtf.gz
        
        mkdir star_index
        tar -xvvf ${star_index} -C star_index --strip-components=1
        
        mkdir tmp
        samtools collate -O -@ 15 ${bam} tmp | samtools fastq -@ ${num_threads} -1 ${prefix}.R1.fastq.gz -2 ${prefix}.R2.fastq.gz -0 /dev/null -s /dev/null -n
		rm -rf tmp
        rm ${bam}
        
        mkdir tmp
        echo ${prefix}.R1.fastq.gz:${prefix}.R2.fastq.gz > s1.txt
        run_rmats --s1 s1.txt --gtf gencode.v38.primary_assembly.annotation.gtf --bi star_index -t paired --libType fr-unstranded --readLength 101 --nthread ${num_threads} --od ${prefix}_rmats --tmp tmp --statoff
    	rm -rf tmp
    }

    output {
        File A3SS_JCEC_output="${prefix}_rmats/A3SS.MATS.JCEC.txt"
        File A5SS_JCEC_output="${prefix}_rmats/A5SS.MATS.JCEC.txt"
        File MXE_JCEC_output="${prefix}_rmats/MXE.MATS.JCEC.txt"
        File RI_JCEC_output="${prefix}_rmats/RI.MATS.JCEC.txt"
        File SE_JCEC_output="${prefix}_rmats/SE.MATS.JCEC.txt"
        
        File A3SS_JC_output="${prefix}_rmats/A3SS.MATS.JC.txt"
        File A5SS_JC_output="${prefix}_rmats/A5SS.MATS.JC.txt"
        File MXE_JC_output="${prefix}_rmats/MXE.MATS.JC.txt"
        File RI_JC_output="${prefix}_rmats/RI.MATS.JC.txt"
        File SE_JC_output="${prefix}_rmats/SE.MATS.JC.txt"
    }

    runtime {
        docker: "docker.io/davidwu20/rmats_turbo:latest"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    meta {
        author: "David Wu"
    }
}


workflow rmats_turbo_workflow {
    call rmats_turbo
}