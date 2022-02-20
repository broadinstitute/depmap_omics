version 1.0

# Given a set of samples, combine segment files into a single file
# more information available at https://open-cravat.readthedocs.io/en/latest/2.-Command-line-usage.html
workflow run_filter_to_maf {
    input {
        File vcf
        String sample_id 
    }
    
    call filter_to_maf {
        input:
        vcf_file=vcf,
        sample_id=sample_id
    }
    output {
        File somatic_maf=filter_to_maf.somatic_maf
    }
}

task filter_to_maf {
    input {
        File vcf_file
        String sample_id

        Int memory = 2
        Int boot_disk_size = 10
        Int num_threads = 1
        Int num_preempt = 5
        Int disk_space = 10
        String docker = "python"
    }

    command{
        git clone https://github.com/broadinstitute/depmap_omics.git
       
        ## filter on gnomad mutations

        ## remove if ... and  .. and ...

        ## merge columns into ... and subset columns


        ## transform to maf file

    }

    output {
        File somatic_maf="${sample_id}_somatic.maf"
    }

    runtime {
        docker: docker
        bootDiskSizeGb: "${boot_disk_size}"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    meta {
        author: "Jeremie Kalfon"
    }
}