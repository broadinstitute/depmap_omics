Whole genome sequencing pipeline for DepMap 
=============================================


Mutation calling pipeline
---------------------------

Mutation annotation and filtering pipeline
-------------------------------------------



CNV calling pipeline
---------------------------

Absolute CNV calling pipeline
-------------------------------

Structural variation calling pipeline
-----------------------------------------



Masking
--------------------------
Due to limitations of short read sequencing, confidence of copy number and mutation calls can be low in regions that are highly repeated or segmentally duplicated. As a result, we use RepeatMasker and Segmental Duplication downloaded from UCSC to mask potentially low quality calls that fall into these regions.

### Segmental Duplication:
The list of segmental duplication intervals is pulled from [UCSC's table browser](http://genome.ucsc.edu/cgi-bin/hgTables), inputs shown in screenshot below:

![](../segdup_ucsc.png)

We then filtered this list by removing duplications where fracMatch < 0.98, as well as duplication pairs where neither segment is on a major contig. Processed bed file can be found at `data/segDup_majorAllele_withAltContigs_98pcFracMatch_merged.bed`.

For the gene-level copy number matrix, we calculated each gene's overlap with segmentally duplicated regions above and masked those where more than 50% of the gene body overlaps with segdup regions.