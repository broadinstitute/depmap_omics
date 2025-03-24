version 1.0
task Manta {
    input {
        String sample_name
        File tumor_bam
        File tumor_bam_index
        File? normal_bam
        File? normal_bam_index
        File ref_fasta
        File ref_fasta_index
        File? interval_bed
        File? interval_bed_index

        Boolean is_cram
        String manta_docker
        String config_manta
        Boolean is_major_contigs_only
        File major_contig_bed="gs://ccleparams/manta_major_contigs.bed.gz"
        File major_contig_bed_index="gs://ccleparams/manta_major_contigs.bed.gz.tbi"

        Int preemptible = 2
        Int max_retries = 0
        Int cpu = 8
        Float mem_per_job_gb = 0.4
        Int additional_disk_gb = 0
    }

    Float jobs_per_cpu = 1.3
    Int num_jobs = round(cpu * jobs_per_cpu)
    Int mem_gb = ceil(num_jobs * mem_per_job_gb)
    Int disk_space = ceil(size(tumor_bam, "GiB")) + 10 + additional_disk_gb

    command {
        EXTENSION="bam"
        EXTENSION_INDEX="bai"
        if [ ${is_cram} = true ]
        then
           EXTENSION="cram"
           EXTENSION_INDEX="crai"
        fi
        
        # sample links
        ln -vs ${tumor_bam} "tumor.$EXTENSION"
        ln -vs ${tumor_bam_index} "tumor.$EXTENSION_INDEX"
        
        if [[ -f "${normal_bam}" ]]; then
            ln -vs ${normal_bam} "normal.$EXTENSION"
            ln -vs ${normal_bam_index} "normal.$EXTENSION_INDEX"
            normal_command_line="--normalBam normal.$EXTENSION"
        fi

        # reference links
        ln -vs ${ref_fasta} reference.fasta
        ln -vs ${ref_fasta_index} reference.fasta.fai
        ln -vs ${major_contig_bed} major_contigs.bed.gz
        ln -vs ${major_contig_bed_index} major_contigs.bed.gz.tbi

        if [ ${is_major_contigs_only} = true ]
        then
            major_contig_line="--callRegions major_contigs.bed.gz"
        fi

        ${config_manta} --tumorBam "tumor.$EXTENSION" \
                        $normal_command_line \
                        --referenceFasta reference.fasta \
                        --runDir . \
                        $major_contig_line

        ./runWorkflow.py --mode local \
                         --jobs ${default=32 cpu} \
                         --memGb ${(${num_jobs} * 2)}

        # change the default names with sample prefix
        if [[ -f "${normal_bam}" ]]; then
           mv results/variants/diploidSV.vcf.gz ${sample_name}.diploidSV.vcf.gz
           mv results/variants/diploidSV.vcf.gz.tbi ${sample_name}.diploidSV.vcf.gz.tbi
           mv results/variants/somaticSV.vcf.gz ${sample_name}.somaticSV.vcf.gz
           mv results/variants/somaticSV.vcf.gz.tbi ${sample_name}.somaticSV.vcf.gz.tbi
        else
           touch ${sample_name}.diploidSV.vcf.gz
           touch ${sample_name}.diploidSV.vcf.gz.tbi
           mv results/variants/tumorSV.vcf.gz ${sample_name}.somaticSV.vcf.gz
           mv results/variants/tumorSV.vcf.gz.tbi ${sample_name}.somaticSV.vcf.gz.tbi
        fi
        mv results/variants/candidateSV.vcf.gz ${sample_name}.candidateSV.vcf.gz
        mv results/variants/candidateSV.vcf.gz.tbi ${sample_name}.candidateSV.vcf.gz.tbi
        mv results/variants/candidateSmallIndels.vcf.gz ${sample_name}.candidateSmallIndels.vcf.gz
        mv results/variants/candidateSmallIndels.vcf.gz.tbi ${sample_name}.candidateSmallIndels.vcf.gz.tbi
    }
    runtime {
        docker: "${manta_docker}"
        memory: "~{mem_gb} GiB"
        cpu: cpu
        disks: "local-disk ~{disk_space} SSD"
        preemptible: preemptible
        maxRetries: max_retries
    }
    output {
        File germline_sv_vcf = "${sample_name}.diploidSV.vcf.gz"
        File germline_sv_vcf_index = "${sample_name}.diploidSV.vcf.gz.tbi"
        File somatic_sv_vcf = "${sample_name}.somaticSV.vcf.gz"
        File somatic_sv_vcf_index = "${sample_name}.somaticSV.vcf.gz.tbi"
        File candidate_sv_vcf = "${sample_name}.candidateSV.vcf.gz"
        File candidate_sv_vcf_index = "${sample_name}.candidateSV.vcf.gz.tbi"
        File candidate_indel_vcf = "${sample_name}.candidateSmallIndels.vcf.gz"
        File candidate_indel_vcf_index = "${sample_name}.candidateSmallIndels.vcf.gz.tbi"
    }
}


task ConvertToBedTabix {
    input {
        File? interval_list
        String output_bed = "interval.bed"
    }

    command <<<
        if [[ "${interval_list}" == *.interval_list ]]
        then
            # Convert Picard-style intervals to BED
            grep -v "^@" ${interval_list} | awk '{$2-=1; print $1,$2,$3,$5}' OFS="\t" > ${output_bed}
        else
            cp ${interval_list} ${output_bed}
        fi

        /usr/gitc/bgzip ${output_bed}
        /usr/gitc/tabix -p bed ${output_bed}.gz
    >>>
    runtime {
        docker: "us-docker.pkg.dev/depmap-omics/public/genomes-in-the-cloud:2.3.3-1525357118"
        memory: "2 GB"
        preemptible: 3
    }
    output {
        File out_interval = "${output_bed}.gz"
        File out_interval_index = "${output_bed}.gz.tbi"
    }
}


workflow MantaSomaticSV {
    input {
        String sample_name
        File tumor_bam
        File tumor_bam_index
        File? normal_bam
        File? normal_bam_index
        File ref_fasta
        File ref_fasta_index
        File? interval_list
        
        String manta_docker
        String config_manta

        Boolean is_exome = defined(interval_list)
        Boolean is_cram
        Boolean is_major_contigs_only
    }

    if (is_exome) {
        call ConvertToBedTabix {
            input: interval_list = interval_list
        }
    }
    
    call Manta {
        input: sample_name = sample_name,
               tumor_bam = tumor_bam,
               tumor_bam_index = tumor_bam_index,
               normal_bam = normal_bam,
               normal_bam_index = normal_bam_index,
               ref_fasta = ref_fasta,
               ref_fasta_index = ref_fasta_index,
               interval_bed = ConvertToBedTabix.out_interval,
               interval_bed_index = ConvertToBedTabix.out_interval_index,
               manta_docker = manta_docker,
               config_manta = config_manta,
               is_cram = is_cram,
               is_major_contigs_only = is_major_contigs_only
    }
    output {
        File germline_sv_vcf = Manta.germline_sv_vcf
        File germline_sv_vcf_index = Manta.germline_sv_vcf_index
        File somatic_sv_vcf = Manta.somatic_sv_vcf
        File somatic_sv_vcf_index = Manta.somatic_sv_vcf_index
        File candidate_sv_vcf = Manta.candidate_sv_vcf
        File candidate_sv_vcf_index = Manta.candidate_sv_vcf_index
        File candidate_indel_vcf = Manta.candidate_indel_vcf
        File candidate_indel_vcf_index = Manta.candidate_indel_vcf_index
    }
}
