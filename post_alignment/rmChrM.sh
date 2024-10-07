#!/bin/bash

# Directory containing the BAM files
bam_dir="/collignon/Tawdarous/atac_seq_lab/results/alignment/bam_files"

# Directory to store the BAM files without mitochondrial reads

output_dir="/collignon/Tawdarous/atac_seq_lab/results/alignment/rmChrM_bam_files"

mkdir -p "$output_dir"  # Ensure the output directory exists

# Number of threads to use for samtools operations
threads=30

# Loop through each BAM file in the directory
for bam_file in "$bam_dir"/*_sorted.bam; do
    # Extract the base name of the file without the directory or .bam extension
    base_name=$(basename "$bam_file" .bam)
    
    # Define the output file name
    output_bam="$output_dir/${base_name}.rmChrM.bam"
    
    # Filter out mitochondrial reads and sort the output using multiple threads
    echo "Processing $bam_file..."
    samtools view -h "$bam_file" | grep -v NC_005089.1 | samtools sort -O bam -@ "$threads" -o "$output_bam" -T "$output_dir/tmp-${base_name}"
    
    echo "Processed $bam_file"
done

echo "All BAM files have been processed."
