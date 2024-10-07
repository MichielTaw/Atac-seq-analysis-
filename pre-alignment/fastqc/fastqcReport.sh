#!/bin/bash

# Base directory for your FASTQ files
fastq_base_dir="/collignon/Tawdarous/atac_seq_lab/data/fastq"

# Output directory for FastQC reports
output_dir="/collignon/Tawdarous/atac_seq_lab/results/fastqc"
mkdir -p "$output_dir"  # Ensure the output directory exists

# Number of threads for FastQC to use
threads=30

# Loop through each sample directory within the FASTQ base directory
for sample_dir in "$fastq_base_dir"/t*; do
    # Navigate into each sample directory
    if [ -d "$sample_dir" ]; then
        echo "Processing directory: $sample_dir"
        
        # Run FastQC for each FASTQ file in the sample directory using multiple threads
        fastqc -t "$threads" -o "$output_dir" "$sample_dir"/*.fq.gz

        echo "FastQC done for $(basename "$sample_dir")"
    else
        echo "Skipping $sample_dir - not a directory"
    fi
done

echo "FastQC analysis is complete for all samples."
