#!/bin/bash

# Define directories and file paths
input_dir="/collignon/Tawdarous/atac_seq_lab/results/post_alignment/encode"
output_dir="/collignon/Tawdarous/atac_seq_lab/results/post_alignment/shift_read_coordinates"
blacklist_file="/collignon/Tawdarous/data/Tosic19/atac_seq_PRJNA369204/raw_data/Blacklist/lists/mm10-blacklist.v2.bed.gz"
threads=30  # Define how many threads each task should use

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Loop through each BAM file in the input directory
for bam_file in "$input_dir"/*.bam; do
    base_name=$(basename "$bam_file" .bam)
    shifted_bam="$output_dir/${base_name}.shifted.bam"

    # Run alignmentSieve with specified threading and other options
    echo "Processing $bam_file..."
    alignmentSieve --numberOfProcessors max  --ATACshift --blackListFileName "$blacklist_file" --bam "$bam_file" -o "$shifted_bam"   

    echo "Processed $shifted_bam"
done

echo "All BAM files have been processed using alignmentSieve."
