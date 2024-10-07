#!/bin/bash

# Define directories
input_dir="/collignon/Tawdarous/atac_seq_lab/results/alignment/rmChrM__AND_duplicate_AND_filtered_quality"
sorted_dir="$input_dir/sorted"

# Ensure the directory for sorted files exists
mkdir -p "$sorted_dir"

# Loop through each filtered BAM file
for filtered_bam in "$input_dir"/*.filtered.bam; do
    base_name=$(basename "$filtered_bam" .filtered.bam)
    sorted_bam="$sorted_dir/${base_name}.sorted.bam"

    # Sort the filtered BAM file
    echo "Sorting $filtered_bam..."
    samtools sort -o  "$sorted_bam" "$filtered_bam"
    if [ $? -eq 0 ]; then
        echo "Successfully sorted $filtered_bam to $sorted_bam"
    else
        echo "Failed to sort $filtered_bam"
    fi
done

echo "All filtered BAM files have been sorted."
