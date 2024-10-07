#!/bin/bash

# Define the source directory where the .bam.bai files are located
source_dir="mtawdarous@rstudio2:/collignon/Tawdarous/atac_seq_lab/results/fastqc"

# Define the destination directory where the .bam.bai files should be copied
destination_dir="/Users/mike/Desktop/atac_seq/fastqc"

# Ensure the destination directory exists
mkdir -p "$destination_dir"

# Find and copy all .bam.bai files from the source to the destination directory
find "$source_dir" -type f -name "*.html" -exec sh -c '
    for file do
        rsync -avz "$file" "'$destination_dir'/" && echo "Copied $file" || echo "Failed to copy $file"
    done
' sh {} +

echo "All HTML FastQC report files have been copied from $source_dir to $destination_dir."
