#!/bin/bash

# Directory containing sorted BAM files
bam_dir="/collignon/Tawdarous/atac_seq_lab/results/alignment/bam_files"

# Loop through each sorted BAM file in the directory
for bam_file in "$bam_dir"/*_sorted.bam; do
    echo "Indexing $bam_file..."
    # Index the BAM file
    samtools index "$bam_file"
    
    # Optional: Check if the indexing was successful
    if [ $? -eq 0 ]; then
        echo "Successfully indexed $bam_file"
    else
        echo "Failed to index $bam_file" >&2
    fi
done

echo "All sorted BAM files have been indexed."
