#!/bin/bash

# Define the directory containing sorted BAM files
sorted_dir="/collignon/Tawdarous/atac_seq_lab/results/alignment/rmChrM__AND_duplicate_AND_filtered_quality/sorted"

# Loop through each sorted BAM file
for sorted_bam in "$sorted_dir"/*.sorted.bam; do
    echo "Indexing $sorted_bam..."
    samtools index "$sorted_bam"
    if [ $? -eq 0 ]; then
        echo "Successfully indexed $sorted_bam"
    else
        echo "Failed to index $sorted_bam"
    fi
done

echo "All sorted BAM files have been indexed."
