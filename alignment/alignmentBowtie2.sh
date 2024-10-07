#!/bin/bash

# Directory containing FASTQ files
fastq_dir="/collignon/Tawdarous/atac_seq_lab/results/trimmed_fastq"

# Output directory for SAM files and logs
output_dir="/collignon/Tawdarous/atac_seq_lab/results/alignment/trimmed_alignment"
log_dir="${output_dir}/logs"
mkdir -p "$output_dir" "$log_dir"

# Bowtie2 Index Base
index_base="/collignon/Tawdarous/ref_genome/GRCm39/GRCm39"

# Process each pair of FASTQ files
find "$fastq_dir" -name "*_L1_1.trimmed.fq.gz" | sed 's/_L1_1.trimmed.fq.gz$//' | while read -r sample; do
    sample_name=$(basename "$sample")
    echo "Processing ${sample_name} at $(date)..."

    # Set up output paths
    sam_output="${output_dir}/${sample_name}.sam"
    unaligned_output="${output_dir}/${sample_name}_unaligned.fq"
    aligned_output="${output_dir}/${sample_name}_aligned.fq"
    log_file="${log_dir}/${sample_name}.log"

    # Check for the presence of input files
    if [[ -f "${sample}_L1_1.trimmed.fq.gz" && -f "${sample}_L1_2.trimmed.fq.gz" ]]; then
        # Running Bowtie2 with additional options
        bowtie2 --local --very-sensitive --no-mixed --no-discordant -X 700 -p 30 -x "$index_base" \
                -1 "${sample}_L1_1.trimmed.fq.gz" -2 "${sample}_L1_2.trimmed.fq.gz" \
                -S "$sam_output" \
                --un-conc "${unaligned_output}" \
                --al-conc "${aligned_output}" \
                2> "$log_file"

        # Check if Bowtie2 ran successfully
        if [ $? -eq 0 ]; then
            echo "Bowtie2 alignment for ${sample_name} completed successfully at $(date)."
        else
            echo "Error: Bowtie2 alignment for ${sample_name} failed at $(date). See log: $log_file"
        fi
    else
        echo "Missing files for ${sample_name}. Check availability of trimmed FASTQ files."
    fi
done

echo "All samples processed."
