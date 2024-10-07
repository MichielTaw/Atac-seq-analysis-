#!/bin/bash

# Directory containing SAM files
sam_dir="/collignon/Tawdarous/atac_seq_lab/results/alignment/trimmed_alignment"

# Output directory for sorted BAM files
bam_dir="/collignon/Tawdarous/atac_seq_lab/results/alignment/bam_files"

# Number of threads to use for sorting
threads=30

# Ensure the output directory exists
mkdir -p "$bam_dir"

# Loop through each SAM file in the directory
for sam_file in "$sam_dir"/*.sam; do
    # Extract the base name of the file without the directory or .sam extension
    base_name=$(basename "$sam_file" .sam)

    # Define the output file name for the unsorted BAM
    unsorted_bam_file="$bam_dir/${base_name}.bam"

    # Convert SAM to BAM
    echo "Converting $sam_file to $unsorted_bam_file"
    samtools view -bS "$sam_file" > "$unsorted_bam_file"
    
    # Check if samtools view was successful
    if [ $? -ne 0 ]; then
        echo "Failed to convert $sam_file to BAM" >&2
        continue # Skip this file and go to the next iteration
    fi

    # Define the output file name for the sorted BAM
    sorted_bam_file="$bam_dir/${base_name}_sorted.bam"

    # Sort BAM using multiple threads
    echo "Sorting $unsorted_bam_file to $sorted_bam_file"
    samtools sort -@ $threads -o "$sorted_bam_file" "$unsorted_bam_file"
    
    # Check if samtools sort was successful
    if [ $? -eq 0 ]; then
        echo "Successfully converted and sorted $sam_file to $sorted_bam_file"
        # Optionally remove the unsorted BAM file
        rm "$unsorted_bam_file"
    else
        echo "Failed to sort $unsorted_bam_file" >&2
    fi
done

echo "All SAM files have been converted to BAM and sorted."
