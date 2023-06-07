import pysam
from pysam import VariantFile

bcf = VariantFile("CosmicCodingMuts.hg38.normal_GRCh38_processed.bcf")
print(bcf.header)

fout = open("cosmic_somatic_gnomad.stat", "w")
tab = '\t'

for record in bcf.fetch():
    tumors = record.info["TUMOR_TYPE"].split('|')
    af = record.info["gnomad_af"]
    max_af = record.info["gnomad_popmax_af"]

    for tumor in tumors:
        tumor = tumor.split(':')
        if len(tumor) >= 3: 
            fout.write(f"{tumor[0]}{tab}{tumor[1]}{tab}{tumor[2]}{tab}{af}{tab}{max_af}\n")
        else:
            print(tumors)
            break
    if len(tumor) < 3: 
        break

fout.close()
