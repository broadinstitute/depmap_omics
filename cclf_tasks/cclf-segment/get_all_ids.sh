#!/bin/bash -ex

for i in $(cut -f 1 23Q2.called.seg | sed 1d | sort -k 1 -u); do
echo $i CCLF WES
done

for i in $(cut -f 1 all.called.seg | sed 1d | sort -k 1 -u); do
echo $i CCLE WES
done

for i in $(cut -f 1 all.called.wgs.seg | sed 1d | sort -k 1 -u); do
echo $i CCLE WGS
done
