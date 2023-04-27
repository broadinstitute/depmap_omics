# Given a set of samples, combine segment files into a single file
workflow msisensor2_workflow {
    call run_msisensor2
}

task run_msisensor2 {
	String sample_id
    File bam
    File bai
    
    Int memory
    Int disk_space
    Int num_threads
    Int num_preempt

    command {
    	set -euo pipefail
    
    	mv ${bai} .
        mv ${bam} .
        
        echo ${bam} | rev | cut -d'/' -f1 | rev > new_bam_path.txt
     
        new_bam_path=$(cat new_bam_path.txt)
    	msisensor2 msi -M /msisensor2/models_hg38 -t $new_bam_path -o ${sample_id}.msisensor2.output
        head -2 ${sample_id}.msisensor2.output | tail -1 | cut -f3 > ${sample_id}.msisensor2.score
    }

    output {
    	Float msisensor2_score=read_float("${sample_id}.msisensor2.score")
        File msisensor2_output="${sample_id}.msisensor2.output"
        File msisensor2_output_dis="${sample_id}.msisensor2.output_dis"
        File msisensor2_output_somatic="${sample_id}.msisensor2.output_somatic"
    }

    runtime {
        docker: "davidwu20/msisensor2:latest"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }
    
    meta {
        author: "David Wu"
    }
}
