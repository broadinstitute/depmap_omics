version 1.0


# this correct know issues in Mutect2 using the script ./fix_mutect2.py

# the AS_filter_status field in the vcf file contains “|” and “,”. 
# But their meanings are swapped compared to other columns in the VCF file, 
# so we swap these back everywhere to keep the same meaning and be able to parse the file easily.

# There is a known issue with filtering combined somatic and germline calls from Mutect2: 
#https://gatk.broadinstitute.org/hc/en-us/community/posts/4404184803227-Mutect2-genotype-germline-sites-filtering-discrepancy-
# Germline variants should not be considered when filtering out clustered events.
# The somatic variant on the left is flagged as a clustered_event because 
# it is near two germline variants. This issue affects about 2.5% of our Mutect2 somatic calls. 
# Unfortunately, we can't just ignore the clustered events filter since it removes a large number of
# sequencing and mapping errors. So we remove the clustered_event flag if less than 2 events in 100bp region are somatic.

workflow run_fix_mutect2 { 
    input {
        File vcf
        String sample_id 
    }
    
    call fix_mutect2 {
        input:
        vcf_file=vcf,
        sample_id=sample_id
    }
    output {
        File vcf_fixed=fix_mutect2.vcf_fixed
    }
}

task fix_mutect2 {
    input {
        File vcf_file
        String sample_id

        Int memory = 4
        Int boot_disk_size = 10
        Int num_threads = 1
        Int num_preempt = 5
        Int disk_space = 50
        Int size_tocheck = 100
        String docker = "python:3.10.7-bullseye"
    }

    command {
        git clone https://github.com/broadinstitute/depmap_omics.git
        pip install wheel
        pip install bgzip

        python depmap_omics/WGS_pipeline/fix_mutect2.py ${vcf_file} ${sample_id}_fixed.vcf.gz ${size_tocheck}
    }

    output {
        File vcf_fixed="${sample_id}_fixed.vcf.gz"
    }

    runtime {
        docker: docker
        bootDiskSizeGb: "${boot_disk_size}"
        memory: "${memory} GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    meta {
        author: "Jeremie Kalfon"
    }
}