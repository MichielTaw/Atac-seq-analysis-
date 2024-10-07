#!/bin/bash

# Define the directory containing the FASTQ files
directory_path="/collignon/Tawdarous/atac_seq_lab/data/fastq"
# Define where the output checksum file should be stored
output_checksum_file="/collignon/Tawdarous/atac_seq_lab/results/md5checksum.txt"

# Clear the output file to ensure it's empty before starting
> "$output_checksum_file"

# Use GNU Parallel to calculate md5sums in parallel
find "$directory_path" -mindepth 2 -maxdepth 2 -name "*.fq.gz" -print0 |
  parallel -0 --lb md5sum {} >> "$output_checksum_file"

echo "All directories have been processed. Checksums are stored in $output_checksum_file."
