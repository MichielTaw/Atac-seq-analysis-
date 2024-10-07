#!/bin/bash

# Directory containing BAM files
bam_dir="/collignon/Tawdarous/atac_seq_lab/results/alignment/rmChrM__AND_duplicate_removed"

# Output CSV file
output_csv="/collignon/Tawdarous/atac_seq_lab/results/alignment/rmChrM__AND_duplicate_removed/read_counts.csv"



threads=30

# Initialize the CSV with headers
echo "Sample,Total Reads,High Quality Reads (>Q30),Percentage" > "$output_csv"

# Loop through each BAM file in the directory
for bam_file in "$bam_dir"/*.bam; do
    # Extract the base name of the file without the directory or '.bam' extension
    base_name=$(basename "$bam_file" .bam)

    # Count total reads
    total_reads=$(samtools view -@ $threads -c "$bam_file")

    # Count high-quality reads (quality > 30)
    high_quality_reads=$(samtools view -@ $threads -c -q 30 "$bam_file")

    # Calculate percentage of high-quality reads
    if [ $total_reads -gt 0 ]; then
        percentage_hq=$(echo "scale=2; $high_quality_reads * 100 / $total_reads" | bc)
    else
        percentage_hq=0
    fi

    # Append results to CSV
    echo "$base_name,$total_reads,$high_quality_reads,$percentage_hq" >> "$output_csv"

    echo "Processed $bam_file"
done

echo "All BAM files have been processed and results are saved to $output_csv."
