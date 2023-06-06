#!/bin/bash -ex

#docker build -t depmapomics:test .

#pytest tests

#cd sandbox/vcf2maf/
#
#docker build -t vcf2maf:test .

pytest tests/depmapomics/tasks/test_vcf2maf.py -s
