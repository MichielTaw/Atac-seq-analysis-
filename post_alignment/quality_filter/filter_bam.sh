#!/bin/bash

# Define directories
input_dir="/collignon/Tawdarous/atac_seq_lab/results/alignment/rmChrM__AND_duplicate_removed"
output_dir="/collignon/Tawdarous/atac_seq_lab/results/alignment/rmChrM__AND_duplicate_AND_filtered_quality"

# Ensure output directory exists
mkdir -p "$output_dir"

# Loop through each BAM file in the input directory
for bam_file in "$input_dir"/*.marked.bam; do
    base_name=$(basename "$bam_file" .bam)
    filtered_bam="$output_dir/${base_name}.filtered.bam"

    # Filter the BAM file
    echo "Filtering $bam_file..."
    samtools view -h -b -f 2 -F 1548 -q 30 "$bam_file" > "$filtered_bam"
    if [ $? -eq 0 ]; then
        echo "Successfully filtered $bam_file"
    else
        echo "Failed to filter $bam_file"
        rm "$filtered_bam"  # Remove potentially incomplete file
    fi
done

echo "All BAM files have been filtered."
