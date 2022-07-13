version 1.0
# PureCN repo: https://github.com/lima1/PureCN

task PureCN {
    
    input {
        # File-related inputs
        String sampleID
        File segFile
        File vcf
        File intervals
        File call_wgd_and_cin_script

        # Method configuration inputs
        String genome="hg38"
        Int maxCopyNumber=8
        Float minPurity=0.90
        Float maxPurity=0.99
        String funSegmentation="Hclust"
        Int maxSegments=100
        String otherArguments="--post-optimize --model-homozygous --min-total-counts 20"

        # Hardware-related inputs
        Int hardware_disk_size_GB = 50
        Int hardware_memory_GB = 4
        Int hardware_preemptible_tries = 2
        Int num_threads = 1
        Int max_retries = 0
    }

    command {
        Rscript /opt/PureCN/PureCN.R \
            --out /cromwell_root/"${sampleID}" \
            --sampleid "${sampleID}" \
            --seg-file "${segFile}" \
            --vcf "${vcf}" \
            --intervals "${intervals}" \
            --genome "${genome}" \
            ${"--max-purity " + maxPurity} \
            ${"--min-purity " + minPurity} \
            ${"--max-copy-number " + maxCopyNumber} \
            ${"--fun-segmentation " + funSegmentation} \
            ${"--max-segments " + maxSegments} \
            ${otherArguments}
        Rscript -e "write.table(read.csv('${sampleID}.csv'),'table.txt',sep='\n',row.names=F,col.names=F,quote=F)"
        Rscript ${call_wgd_and_cin_script} "${sampleID}_loh.csv" "${sampleID}.csv"
    }

    Array[String] table = read_lines('${sampleID}.csv') # maybe wdl has builtin read_csv files?
    Array[String] wgd_table = read_lines("out.txt")

    output {
        File solutions_pdf = "${sampleID}.pdf"
        File chromosomes_pdf = "${sampleID}_chromosomes.pdf"
        File rds = "${sampleID}.rds"
        File dnacopy = "${sampleID}_dnacopy.seg"
        File variants = "${sampleID}_variants.csv"
        File loh = "${sampleID}_loh.csv"
        File genes = "${sampleID}_genes.csv"
        File segmentation = "${sampleID}_segmentation.pdf"
        File log = "${sampleID}.log"
        File selected_solution = "${sampleID}.csv"
        File local_optima_pdf = "${sampleID}_local_optima.pdf"
        String purity = table[1]
        String ploidy = table[2]
        String contamination = table[4]
        String flagged = table[5]
        String curated = table[7]
        String comment = table[8]
        String wgd = wgd_table[0]
        String loh_fraction = wgd_table[1]
        String cin = wgd_table[2]
        String cin_allele_specific = wgd_table[3]
        String cin_ploidy_robust = wgd_table[4]
        String cin_allele_specific_ploidy_robust = wgd_table[5]
    }

    runtime {
        docker: "markusriester/purecn:2.2.0"
        bootDiskSizeGb: 32
        disks: "local-disk ${hardware_disk_size_GB} HDD"
        memory: "${hardware_memory_GB} GB"
        cpu: "${num_threads}"
        continueOnReturnCode: true
        preemptible: "${hardware_preemptible_tries}"
        maxRetries: "${max_retries}"
    }

}

workflow run_PureCN {

    input {
        String sampleID
        File segFile
        File vcf
        File intervals
        File call_wgd_and_cin_script

        # Method configuration inputs
        String genome
    }

    call PureCN {
        input:
            sampleID=sampleID,
            segFile=segFile,
            vcf=vcf,
            intervals=intervals,
            call_wgd_and_cin_script=call_wgd_and_cin_script,
            genome=genome
    }

    output {
        File PureCN_solutions_pdf = PureCN.solutions_pdf
        File PureCN_chromosomes_pdf = PureCN.chromosomes_pdf
        File PureCN_rds = PureCN.rds
        File PureCN_dnacopy = PureCN.dnacopy
        File PureCN_variants = PureCN.variants
        File PureCN_loh = PureCN.loh
        File PureCN_genes = PureCN.genes
        File PureCN_segmentation = PureCN.segmentation
        File PureCN_log = PureCN.log
        File PureCN_selected_solution = PureCN.selected_solution
        File PureCN_local_optima_pdf = PureCN.local_optima_pdf
        String PureCN_purity = PureCN.purity
        String PureCN_ploidy = PureCN.ploidy
        String PureCN_contamination = PureCN.contamination
        String PureCN_flagged = PureCN.flagged
        String PureCN_curated = PureCN.curated
        String PureCN_comment = PureCN.comment
        String PureCN_wgd = PureCN.wgd
        String PureCN_loh_fraction = PureCN.loh_fraction
        String PureCN_cin = PureCN.cin
        String PureCN_cin_allele_specific = PureCN.cin_allele_specific
        String PureCN_cin_ploidy_robust = PureCN.cin_ploidy_robust
        String PureCN_cin_allele_specific_ploidy_robust = PureCN.cin_allele_specific_ploidy_robust
    }
}