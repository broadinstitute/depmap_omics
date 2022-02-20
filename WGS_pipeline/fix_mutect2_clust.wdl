version 1.0

# Given a set of samples, combine segment files into a single file
# more information available at https://open-cravat.readthedocs.io/en/latest/2.-Command-line-usage.html
workflow run_fix_mutect_clust {
    input {
        File vcf
        String sample_id 
    }
    
    call fix_mutect_clust {
        input:
            vcf_file=vcf,
            sample_id=sample_id
    }
    output {
        File vcf_fixed=fix_mutect_clust.vcf_fixed
    }
}

task fix_mutect_clust {
    input {
        File vcf_file
        String sample_id

        Int memory = 2
        Int boot_disk_size = 10
        Int num_threads = 1
        Int num_preempt = 5
        Int disk_space = 10
        String docker = "jkobject/simple_r"
    }

    command{
        git clone https://github.com/broadinstitute/depmap_omics.git

        Rscript WGS_pipeline/correct_mutect2_clusteredevent.R ${vcf_file} ${sample_id}
    }
    output {
        File vcf_fixed="${sample_id}_clustercorrected.vcf.gz"
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