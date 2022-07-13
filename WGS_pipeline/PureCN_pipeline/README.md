# PureCN absolute copy number calling

## Resources

[PureCN best paractices](https://bioconductor.org/packages/devel/bioc/vignettes/PureCN/inst/doc/Quick.html)

## Reference files

PureCN interval files with mappability were generated with following commands. The interval files used in DepMap data generation can also be found in `gs://ccleparams/references/PureCN_intervals`.

**WES AGILENT**

`Rscript $PURECN/IntervalFile.R --in-file agilent_hg38_lifted_chrXY.no_header.bed --fasta ~/Data/VCFs/Liftover/hg38.fa --out-file agilent_hg38_intervals.txt --genome hg38 --mappability GCA_000001405.15_GRCh38_no_alt_analysis_set_100.bw`

**WES ICE**

`Rscript $PURECN/IntervalFile.R --in-file ice_hg38_lifted_chrXY.no_header.bed --fasta ~/Data/VCFs/Liftover/hg38.fa --out-file ice_hg38_intervals.txt --genome hg38 --mappability GCA_000001405.15_GRCh38_no_alt_analysis_set_100.bw`

**WGS**

Here we uniformly sample 2% of the genome to use for absolution copy number inference.

`Rscript make_wgs_intervals.R`

We have also created new intervals from the recommended pureCN workflow, including mappability information on both the wgs intervals and wes intervals
`Rscript $PURECN/IntervalFile.R --in-file wgs_hg38_intervals.bed --fasta ~/Data/VCFs/Liftover/hg38.fa --out-file wgs_hg38_2_percent_intervals.txt --genome hg38 --mappability GCA_000001405.15_GRCh38_no_alt_analysis_set_100.bw`

## Terra workflows

[PureCN-AGILENT](https://app.terra.bio/#workspaces/broad-firecloud-ccle/DepMap_WES_CN_hg38/workflows/broad-firecloud-ccle/PureCN-AGILENT) for AGILENT

[PureCN-AGILENT](https://app.terra.bio/#workspaces/broad-firecloud-ccle/DepMap_WES_CN_hg38/workflows/broad-firecloud-ccle/PureCN-ICE)  for ICE

[PureCN](https://app.terra.bio/#workspaces/broad-firecloud-ccle/DepMap_WGS_CN/workflows/broad-firecloud-ccle/PureCN) for WGS

## Manual curation

About 10% of PureCN calls need to be manually curated. The  [PureCN_Curation](https://github.com/broadinstitute/depmap_omics/blob/master/WGS_pipeline/PureCN_Curation.ipynb) notebook should be used to select the solutions that require curation, download the solution PDFs, and then update the Terra workspace to reflect manual changes. Detailed curation guidelines can be found [here](https://docs.google.com/document/d/1Rte0xKK3ZE_UV6MWepdXRIbAehUJg8FuLaDckrWhPTQ/edit).

Once manual curation is complete the PureCN output files need to be updated to reflect the newly selected solution. To do this, run [PureCN_update_solution](https://app.terra.bio/#workspaces/broad-firecloud-ccle/DepMap_WES_CN_hg38/workflows/colganwi/PureCN_update_solution) on the curated samples.

## Whole genome doubling

WGD is determined using this formula: -2*loh_frac + 3 < ploidy. The call_wgd.R script does step and is part of the PureCN and PureCN_update_solution workflows.
