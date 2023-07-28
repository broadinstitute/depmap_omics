
SELECT likely_lof,dbsnp_caf,CDS_ID,civic_score FROM `depmap-omics.maf_staging.merged_maf` where civic_score >8 LIMIT 1000


SELECT COUNTIF(civic_score > 0) AS CIVIC_thresh0, COUNTIF(civic_score >8) AS CIVIC_thresh8, COUNTIF(civic_score >20) as CIVIC_thresh20 FROM `depmap-omics.maf_staging.merged_maf`

