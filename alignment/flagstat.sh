
#!/bin/bash

# Path to the directory containing sorted BAM files
bam_path="/collignon/Tawdarous/atac_seq_lab/results/alignment/bam_files"

# Path to the directory where the flagstat results will be saved
output_dir="/collignon/Tawdarous/atac_seq_lab/results/alignment/flagstats"

# Set the number of threads to use for samtools
threads=30

# Check if the BAM directory exists
if [ ! -d "$bam_path" ]; then
    echo "Directory does not exist: $bam_path"
    exit 1
fi

# Check if there are any sorted BAM files in the directory
if [ $(ls "$bam_path"/*_sorted.bam 2> /dev/null | wc -l) -eq 0 ]; then
    echo "No sorted BAM files found in $bam_path"
    exit 1
fi

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Loop through each sorted BAM file in the directory
for sorted_bam in "$bam_path"/*_sorted.bam; do
    # Extract the base name of the file without the directory or '.bam' extension
    base_name=$(basename "$sorted_bam" .bam)

    # Output file for flagstat in the designated output directory
    output_file="${output_dir}/${base_name}.flagstat"

    echo "Processing $sorted_bam..."
    # Run samtools flagstat with multi-threading and direct the output to a file
    samtools flagstat -@ $threads "$sorted_bam" > "$output_file"

    # Check if samtools flagstat was successful
    if [ $? -eq 0 ]; then
        echo "Successfully generated flagstat for $sorted_bam"
    else
        echo "Failed to generate flagstat for $sorted_bam" >&2
    fi
done

echo "All sorted BAM files have been processed for flagstats."
