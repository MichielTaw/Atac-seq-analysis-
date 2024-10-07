#!/bin/bash

# Define directories
input_dir="/collignon/Tawdarous/atac_seq_lab/results/alignment/rmChrM__AND_duplicate_AND_filtered_quality"
sorted_dir="$input_dir/sorted"

# Ensure the directory for sorted files exists
mkdir -p "$sorted_dir"

# List of specific samples to sort
declare -a samples=("t24_WNT-06h_L_rep1_sorted.rmChrM_RG.marked.filtered.bam")

# Loop through each specified sample
for sample in "${samples[@]}"; do
    filtered_bam="$input_dir/$sample"
    base_name=$(basename "$filtered_bam" .filtered.bam)
    sorted_bam="$sorted_dir/${base_name}.sorted.bam"

    # Check if the filtered BAM file exists
    if [ -f "$filtered_bam" ]; then
        # Sort the filtered BAM file
        echo "Sorting $filtered_bam..."
        samtools sort -o "$sorted_bam" "$filtered_bam"
        if [ $? -eq 0 ]; then
            echo "Successfully sorted $filtered_bam to $sorted_bam"
        else
            echo "Failed to sort $filtered_bam"
        fi
    else
        echo "Filtered BAM file does not exist: $filtered_bam"
    fi
done

echo "Specified filtered BAM files have been sorted."
