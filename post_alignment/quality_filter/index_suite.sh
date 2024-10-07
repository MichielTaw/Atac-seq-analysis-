#!/bin/bash

# Define the directory containing sorted BAM files
sorted_dir="/collignon/Tawdarous/atac_seq_lab/results/alignment/rmChrM__AND_duplicate_AND_filtered_quality/sorted"

# Specify the exact file to index
sorted_bam="$sorted_dir/t24_WNT-06h_L_rep1_sorted.rmChrM_RG.marked.sorted.bam"

echo "Indexing $sorted_bam..."
samtools index "$sorted_bam"
if [ $? -eq 0 ]; then
    echo "Successfully indexed $sorted_bam"
else
    echo "Failed to index $sorted_bam"
fi

echo "Indexing operation completed."

