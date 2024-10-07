#!/bin/bash

# Directory containing sorted BAM files
bam_dir="/collignon/Tawdarous/atac_seq_lab/results/alignment/rmChrM_bam_files"

# Output directory for idxstats files
output_dir="/collignon/Tawdarous/atac_seq_lab/results/alignment/rmChrM_bam_files/idx_rmChrM"

mkdir -p "$output_dir"  # Create the output directory if it doesn't exist

# Loop through each sorted BAM file in the directory
for bam_file in "$bam_dir"/*rmChrM.bam; do
    # Extract the base name of the file without the directory or _sorted.bam extension
    base_name=$(basename "$bam_file" "_sorted.bam")
    
    # Define the output file for idxstats
    idxstats_file="$output_dir/${base_name}.idxstats"
    
    echo "Generating index stats for $bam_file..."
    # Generate index statistics and redirect output to a file
    samtools idxstats "$bam_file" > "$idxstats_file"
    
    # Optional: Verify that the idxstats were successfully generated
    if [ $? -eq 0 ]; then
        echo "Successfully generated index stats for $bam_file"
    else
        echo "Failed to generate index stats for $bam_file" >&2
    fi
done

echo "All sorted BAM files have been processed for index stats."
