# Given a set of samples, aggregate the MSISensor2 repeats per site
version 1.0

workflow microsatellite_repeats {
    input {
    	Array[File] msisensor2_output_dis
        String sample_set_id
    
        Int memory=6
        Int num_threads=1
        Int disk_space=20
        Int num_preempt=1
    }

    call find_n_repeats {
        input:
            msisensor2_output_dis=msisensor2_output_dis,
            sample_set_id=sample_set_id,
            memory=memory,
            num_threads=num_threads,
            disk_space=disk_space,
            num_preempt=num_preempt
    }

    output {
        File ms_repeats=find_n_repeats.ms_repeats
    }
}


task find_n_repeats {
    input {
    	Array[File] msisensor2_output_dis
        String sample_set_id
    
        Int memory=6
        Int num_threads=1
        Int disk_space=20
        Int num_preempt=1
    }

    
    command{
        git clone --branch microsat https://github.com/broadinstitute/depmap_omics.git
        python3 depmap_omics/WGS_pipeline/aggregate_msisensor2_repeats.py ${write_lines(msisensor2_output_dis)} ${sample_set_id}
    }

    output {
        File ms_repeats="${sample_set_id}.microsatellite_repeats.csv"
    }

    runtime {
        docker: "iboylebroad/microsat:latest"
        cpu: "${num_threads}"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} SSD"
        preemptible: "${num_preempt}"
    }

    meta {
        author: "Isabella Boyle"
    }
}
