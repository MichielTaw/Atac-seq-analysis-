#!/bin/bash

# Directory containing original BAM files
bam_dir="/collignon/Tawdarous/atac_seq_lab/results/alignment/rmChrM_bam_files"
# Output directory for BAM files with added read groups
output_dir="/collignon/Tawdarous/atac_seq_lab/results/alignment/reading_group"
mkdir -p "$output_dir"

# Set the path to the Picard jar file
picard_jar_path="/home/mtawdarous/softwares/picard/build/libs/picard-3.1.1-20-gc8b2c06b2-SNAPSHOT-all.jar"

# Number of concurrent jobs
max_jobs=30
running_jobs=0

# Loop over each BAM file in the input directory
for bam_file in "$bam_dir"/*rmChrM.bam; do
    # Wait for other jobs if the limit is reached
    if (( running_jobs >= max_jobs )); then
        wait -n  # Wait for any job to finish
        ((running_jobs--))
    fi

    # Extract the base name of the file without the directory or file extension
    base_name=$(basename "$bam_file" .bam)
    output_bam="$output_dir/${base_name}_RG.bam"

    # Execute Picard's AddOrReplaceReadGroups in background
    java -jar "$picard_jar_path" AddOrReplaceReadGroups \
        I="$bam_file" \
        O="$output_bam" \
        RGID="$base_name" \
        RGLB="${base_name}_library" \
        RGPL="ILLUMINA" \
        RGPU="${base_name}_unit" \
        RGSM="${base_name}_sample" &
    
    ((running_jobs++))
done

# Wait for all background jobs to finish
wait
