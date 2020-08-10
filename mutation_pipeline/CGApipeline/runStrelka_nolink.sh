#!/bin/bash -l


libdir=$1
tumorBam=$2
normalBam=$3
ref=$4
indiv=$5
config=$6

echo "start"
date

echo $@

STRELKA_INSTALL_DIR=${libdir}/strelka_workflow_1.0.11

$STRELKA_INSTALL_DIR/bin/configureStrelkaWorkflow_cga.pl --normal=$normalBam --tumor=$tumorBam --ref=$ref --config=${libdir}/${config} --output-dir=./analysis

make -j 4 -C ./analysis

cp -v analysis/results/all.somatic.indels.vcf $indiv.all.somatic.indels.vcf 
cp -v analysis/results/passed.somatic.indels.vcf $indiv.passed.somatic.indels.vcf 
cp -v analysis/results/all.somatic.snvs.vcf $indiv.all.somatic.snvs.vcf 
cp -v analysis/results/passed.somatic.snvs.vcf $indiv.passed.somatic.snvs.vcf 

echo "done"
date





