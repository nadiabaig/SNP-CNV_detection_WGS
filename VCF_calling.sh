reference_genome="/mnt/d/PotatoTools/Agria-Ref/Agria_21082020/Potato/dAg_v1.0/Agria_assembly_final_2020_21_08.fasta"

#variant calling per chromosome
bam_list="/mnt/d/F23A430001979_PLAheanR/soapnuke/clean/list.bam"
freebayes -f $reference_genome -r Chr01 -p 4 -C 5 -m 30 -q 20 -R 0 -S 0 -F 0.2 -j -H -P 0.95 -L $bam_list  --use-best-n-alleles 4  > Chr01.vcf &&
freebayes -f $reference_genome -r Chr02 -p 4 -C 5 -m 30 -q 20 -R 0 -S 0 -F 0.2 -j -H -P 0.95 -L $bam_list  --use-best-n-alleles 4  > Chr02.vcf &&
freebayes -f $reference_genome -r Chr03 -p 4 -C 5 -m 30 -q 20 -R 0 -S 0 -F 0.2 -j -H -P 0.95 -L $bam_list  --use-best-n-alleles 4  > Chr03.vcf 

