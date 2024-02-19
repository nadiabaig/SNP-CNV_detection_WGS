#!/bin/bash

#indexing reference genome
reference_genome="/dAg_v1.0/Agria_assembly_final_2020_21_08.fasta"
#bwa index $reference_genome
#
base_dir="/mnt/d/PLAheanR/soapnuke/clean"

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

java -jar /mnt/d/gatk-4.2.0.0/gatk-package-4.2.0.0-local.jar \
        -T RealignerTargetCreator \
        -R $reference_genome \
        -I "$result_dir"/deduplicated_reads.bam \
        -o "$result_dir"/target_intervals.list

wait
#Create Realignment Targets. Local realignment around indels to reduce wrong variant calls
java -jar /gatk-4.2.0.0/gatk-package-4.2.0.0-local.jar \
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
