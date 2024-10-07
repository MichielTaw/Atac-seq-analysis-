

#!/bin/bash

# Define directories
input_dir="/collignon/Tawdarous/atac_seq_lab/results/alignment/rmChrM__AND_duplicate_AND_filtered_quality/sorted"
output_dir="/collignon/Tawdarous/atac_seq_lab/results/post_alignment/encode"
metrics_dir="/collignon/Tawdarous/atac_seq_lab/results/post_alignment/encode/metrics"
blacklist_bed="/collignon/Tawdarous/data/Tosic19/atac_seq_PRJNA369204/raw_data/Blacklist/lists/mm10-blacklist.v2.bed.gz"  # Path to the blacklist BED file

# Create necessary directories
mkdir -p "$output_dir" "$metrics_dir" "$output_dir/tmp"

# Number of cores to use for parallel operations
threads=30

export TMPDIR=$output_dir/tmp

process_bam() {
    bam_file=$1
    base_name=$(basename "$bam_file" .bam)
    
    # Define output file names
    output_bam="$output_dir/${base_name}.blacklist-filtered.bam"
    temp_bam="$output_dir/${base_name}.tmp.bam"
    metrics_file="$metrics_dir/${base_name}.dup.metrics"
    
    # Filtering with bedtools
    bedtools intersect -nonamecheck -v -abam "$bam_file" -b "$blacklist_bed" > "$temp_bam"
    
    # Sorting and indexing the bam file
    samtools sort -@ $threads -o "$output_bam" "$temp_bam"
    samtools index -@ $threads "$output_bam"
    
    # Remove temporary file
    rm "$temp_bam"
    
    echo "Processed $bam_file"
}

export -f process_bam
export threads
export output_dir
export blacklist_bed

# Find all bam files and process them in parallel
find "$input_dir" -name "*.bam" | parallel -j $threads process_bam

echo "All BAM files have been processed for duplicate marking."
