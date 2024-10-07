#!/bin/bash

# Define directories
input_dir="/collignon/Tawdarous/atac_seq_lab/data/fastq"  # Directory containing the raw FASTQ files
output_dir="/collignon/Tawdarous/atac_seq_lab/results/trimmed_fastq"  # Directory to save the trimmed FASTQ files
report_dir="/collignon/Tawdarous/atac_seq_lab/results/trimmed_fastq/report"  # Directory to save the reports

# Create necessary directories
mkdir -p "$output_dir" "$report_dir"

# Export necessary variables and functions
export input_dir output_dir report_dir

# Function to process each sample
process_fastq_pair() {
    dir_path="$1"
    sample_name=$(basename "$dir_path")
    r1="$dir_path/${sample_name}_L1_1.fq.gz"
    r2="$dir_path/${sample_name}_L1_2.fq.gz"
    trimmed_r1="$output_dir/${sample_name}_L1_1.trimmed.fq.gz"
    trimmed_r2="$output_dir/${sample_name}_L1_2.trimmed.fq.gz"
    json_report="$report_dir/${sample_name}.fastp.json"
    html_report="$report_dir/${sample_name}.fastp.html"

    fastp -i "$r1" -I "$r2" -o "$trimmed_r1" -O "$trimmed_r2" --detect_adapter_for_pe -j "$json_report" -h "$html_report" -w 60


    echo "Processed $sample_name: $r1, $r2"
}

export -f process_fastq_pair

# Find all directories containing FASTQ files and process them
find "$input_dir" -mindepth 1 -maxdepth 1 -type d | parallel process_fastq_pair

echo "All samples have been processed."

