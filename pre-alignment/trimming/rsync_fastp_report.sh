#!/bin/bash
# Define the source directory on the remote server where the .html files are located
source_dir="mtawdarous@rstudio2:/collignon/Tawdarous/atac_seq_lab/results/trimmed_fastq/report"

# Define the destination directory on the local machine where the files should be copied
destination_dir="/Users/mike/Desktop/atac_seq/trimming_report"

# Use rsync to securely copy all .html files from the source to the destination directory
rsync -avz -e ssh "${source_dir}/*.html" "${destination_dir}/"
rsync -avz -e ssh "${source_dir}/*.json" "${destination_dir}/"


echo "All html and json  fastp report files have been copied from $source_dir to $destination_dir."
