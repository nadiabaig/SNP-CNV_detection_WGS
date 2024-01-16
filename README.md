# SNP-CNV_detection_WGS
Bash pipeline for Single Nucleotide Polymorphism (SNP) and CNV detection in WGS data

## Installations

# Fastqc
   
    sudo apt-get update
    sudo apt-get -y install fastqc 

# Picard

    wget -O picard.jar https://github.com/broadinstitute/picard/releases/download/2.26.2/picard.jar

# GATK
wget https://github.com/broadinstitute/gatk/releases/download/4.2.0.0/gatk-4.2.0.0.zip
unzip gatk-4.2.0.0.zip

# Quality checking using FastQC

```
#!/bin/bash

#Specify the base directory where your folders are located
base_dir="/mnt/d/PLAheanR/soapnuke/clean"

#Iterate through folders starting with RC_
for folder in "$base_dir"/RC_*;
do
#Extract folder name without the base directory
folder_name=$(basename "$folder")
#Run FastQC on each folder
fastqc -o "$base_dir" "$folder"/*.fq.gz

#Create a new directory for storing results
result_dir="$base_dir/Fastqc_Quality"

mkdir -p "$result_dir"

#Move FastQC results to the new directory
mv "$folder"/*.zip "$result_dir"
mv "$folder"/*.html "$result_dir"

echo "Quality control results for $folder_name saved in $result_dir"
done
```
# BAM file generation, PCR duplicates removal and BAM files preprocessing

```
#!/bin/bash

#indexing reference genome
reference_genome="/mnt/d/PotatoTools/Agria-Ref/Agria_21082020/Potato/dAg_v1.0/Agria_assembly_final_2020_21_08.fasta"
#bwa index $reference_genome
#
base_dir="/mnt/d/F23A430001979_PLAheanR/soapnuke/clean"

# Iterate through folders starting with RC_
for folder in "$base_dir"/RC_*;
do

# Extract folder name without the base directory
folder_name=$(basename "$folder")

# Create a new directory for storing intermediate files and results
result_dir="$base_dir/${folder_name}_results"
mkdir -p "$result_dir"

# BWA MEM alignment
bwa mem -M -t 14 "$reference_genome" "$folder"/*_1.fq.gz "$folder"/*_2.fq.gz | samtools view -b - > "$result_dir"/aligned_reads.bam

wait
# Sorting the BAM file
samtools sort -o "$result_dir"/sorted_reads.bam "$result_dir"/aligned_reads.bam

wait

java -jar picard.jar CollectAlignmentSummaryMetrics R=$reference_genome I="$result_dir"/sorted_reads.bam O=$result_dir/alignment_metrics.txt

wait


# Marking PCR duplicates
java -jar picard.jar MarkDuplicates INPUT="$result_dir"/sorted_reads.bam OUTPUT="$result_dir"/deduplicated_reads.bam METRICS_FILE="$result_dir"/deduplication_metrics.txt REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT

wait

#Making bam index for further processing
java -jar picard.jar BuildBamIndex INPUT="$result_dir"/deduplicated_reads.bam

wait

java -jar /mnt/d/F23A430001979_PLAheanR/soapnuke/clean/gatk-4.2.0.0/gatk-package-4.2.0.0-local.jar \
        -T RealignerTargetCreator \
        -R $reference_genome \
        -I "$result_dir"/deduplicated_reads.bam \
        -o "$result_dir"/target_intervals.list

wait
#Create Realignment Targets. Local realignment around indels to reduce wrong variant calls
java -jar /mnt/d/F23A430001979_PLAheanR/soapnuke/clean/gatk-4.2.0.0/gatk-package-4.2.0.0-local.jar \
        -T IndelRealigner \
        -R $reference_genome \
        -I "$result_dir"/deduplicated_reads.bam \
        -targetIntervals "$result_dir"/target_intervals.list \
        -o "$result_dir"/realigned_reads.bam

wait
java -jar picard.jar BuildBamIndex INPUT="$result_dir"/realigned_reads.bam

wait

# Indexing the deduplicated BAM file
samtools index "$result_dir"/realigned_reads.bam

wait

# getting uniquely aligned reads
samtools view -h -q 20 -b "$result_dir"/realigned_reads.bam > "$result_dir"/filtered_realigned_reads.bam

wait
samtools index "$result_dir"/filtered_realigned_reads.bam

wait
samtools stats "$result_dir"/filtered_realigned_reads.bam > "$result_dir"/filtered_realigned_reads.stats

done
```
    
