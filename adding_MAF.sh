#!/bin/bash


for i in $(ls *_filtered.vcf)
do

        bcftools +fill-tags ${i} -Ov -o ${i}.gz -- -t all

done
