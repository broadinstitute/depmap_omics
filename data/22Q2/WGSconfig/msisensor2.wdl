# Given a set of samples, combine segment files into a single file
version 1.0

task msisensor2 {
	
    input {
        String sample_id
        File bam
        File bai
        
        Int memory # needs to be default
        Int disk_space
        Int num_threads
        Int num_preempt
    }

    command {
    	set -euo pipefail
    
    	mv ${bai} .
        mv ${bam} .
        
        echo ${bam} | rev | cut -d'/' -f1 | rev > new_bam_path.txt
     
        new_bam_path=$(cat new_bam_path.txt)
    	msisensor2 msi -M /msisensor2/models_hg38 -t $new_bam_path -o ${sample_id}.msisensor2.output # one param per line, also what is /msisensor2/models_hg38 ?
        head -2 ${sample_id}.msisensor2.output | tail -1 | cut -f3 > ${sample_id}.msisensor2.score
    }

    output {
    	Float msisensor2_score=read_float("${sample_id}.msisensor2.score")
        File msisensor2_output="${sample_id}.msisensor2.output"
        File msisensor2_output_dis="${sample_id}.msisensor2.output_dis"
        File msisensor2_output_somatic="${sample_id}.msisensor2.output_somatic"
    }

    runtime {
        docker: "davidwu20/msisensor2:1"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }
    
    meta {
        author: "David Wu"
    }
}

workflow msisensor2_workflow {
    input {
        String sample_id
        File bam
        File bai
        
        Int memory
        Int disk_space
        Int num_threads
        Int num_preempt
    }

    call msisensor2 {
        input {
            sample_id=sample_id,
            bam=bam,
            bai=bai,
            memory=memory,
            disk_space=disk_space,
            num_threads=num_threads,
            num_preempt=num_preempt
        }
    }

    output {
        Float msisensor2_score=msisensor2.msisensor2_score
        File msisensor2_output=msisensor2.msisensor2_output
        File msisensor2_output_dis=msisensor2.msisensor2_output_dis
        File msisensor2_output_somatic=msisensor2.msisensor2_output_somatic
    }
}