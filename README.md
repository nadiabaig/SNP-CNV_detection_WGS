# SNP-CNV_detection_WGS
Bash pipeline for Single Nucleotide Polymorphism (SNP) and CNV detection in WGS data

## Installations

# Fastqc
   
    sudo apt-get update
    sudo apt-get -y install fastqc 





# Quality checking using FastQC

#!/bin/bash

# Specify the base directory where your folders are located
base_dir="/mnt/d/F23A430001979_PLAheanR/soapnuke/clean"

# Iterate through folders starting with RC_
for folder in "$base_dir"/RC_*;
do
# Extract folder name without the base directory
folder_name=$(basename "$folder")
# Run FastQC on each folder
fastqc -o "$base_dir" "$folder"/*.fq.gz

# Create a new directory for storing results
result_dir="$base_dir/Fastqc_Quality"

mkdir -p "$result_dir"

# Move FastQC results to the new directory
mv "$folder"/*.zip "$result_dir"
mv "$folder"/*.html "$result_dir"

echo "Quality control results for $folder_name saved in $result_dir"
done

    
