#!/bin/bash

# Define directories and file paths
input_dir="/collignon/Tawdarous/atac_seq_lab/results/alignment/reading_group"
output_dir="/collignon/Tawdarous/atac_seq_lab/results/alignment/rmChrM__AND_duplicate_removed"
metrics_dir="/collignon/Tawdarous/atac_seq_lab/results/alignment/rmChrM__AND_duplicate_removed/metrics"
mkdir -p "$output_dir" "$metrics_dir" "$output_dir/tmp"

# Set the path to the Picard jar file
picard_jar_path="/home/mtawdarous/softwares/picard/build/libs/picard-3.1.1-20-gc8b2c06b2-SNAPSHOT-all.jar"

# Increase Java heap size and ensure use of JDK's deflater/inflater
export _JAVA_OPTIONS="-Xmx64g -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.compression_level=1 -Dsamjdk.use_jdk_deflater=true -Dsamjdk.use_jdk_inflater=true"


# Loop through each BAM file that has mitochondrial reads removed
for bam_file in "$input_dir"/*rmChrM_RG.bam; do
    # Extract the base name of the file without the directory or file extension
    base_name=$(basename "$bam_file" .bam)
    
    # Define the output file names
    output_bam="$output_dir/${base_name}.marked.bam"
    metrics_file="$metrics_dir/${base_name}.dup.metrics"
    
    # Run Picard MarkDuplicates
    echo "Processing $bam_file..."
    java -jar "$picard_jar_path" MarkDuplicates \
        QUIET=true \
        INPUT="$bam_file" \
        OUTPUT="$output_bam" \
        METRICS_FILE="$metrics_file" \
        REMOVE_DUPLICATES=false \
        CREATE_INDEX=true \
        VALIDATION_STRINGENCY=LENIENT \
        TMP_DIR="$output_dir/tmp" \
        USE_JDK_DEFLATER=true \
        USE_JDK_INFLATER=true
    
    echo "Processed $bam_file"
done

echo "All BAM files have been processed for duplicate marking."
