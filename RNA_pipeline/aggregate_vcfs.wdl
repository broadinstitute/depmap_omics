# Given a set of samples, combine vcf files into a single file
workflow aggregate_vcfs {
    call aggregate
}


task aggregate {
	Array[File] vcf_files
    Array[File]? vcf_indexes
    String sample_set_id    
    Int memory
    Int num_threads
    Int disk_space
    Int num_preempt
    String option

    
    command{
        bcftools merge -m ${option} --force-samples --threads ${num_threads} ${sep=" " vcf_files} -o "${sample_set_id}.called.vcf"
    }

    output {
        File merged_vcf="${sample_set_id}.called.vcf"
    }

    runtime {
        docker: "us-docker.pkg.dev/depmap-omics/public/bcftools:v1.9-1-deb_cv1"
        cpu: "${num_threads}"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        preemptible: "${num_preempt}"
    }

    meta {
        author: "Jeremie Kalfon"
    }
}