#!/bin/bash

input_folder="/collignon/Tawdarous/atac_seq_lab/results/post_alignment/shift_read_coordinates"
output_folder="/collignon/Tawdarous/atac_seq_lab/results/post_alignment/shift_chr_renamed"

# Create output folder if it doesn't exist
mkdir -p "$output_folder"

# Loop through all .bam files in the input folder
for bam_file in "$input_folder"/*.bam; do
    # Get the base name of the file
    base_name=$(basename "$bam_file" .bam)
    
    # Create the output file name
    output_file="$output_folder/${base_name}_chr_renamed.bam"
    
    # Process the BAM file
    samtools view -H "$bam_file" | \
    sed -e 's/SN:NC_000067.7/SN:chr1/' \
        -e 's/SN:NC_000068.8/SN:chr2/' \
        -e 's/SN:NC_000069.7/SN:chr3/' \
        -e 's/SN:NC_000070.7/SN:chr4/' \
        -e 's/SN:NC_000071.7/SN:chr5/' \
        -e 's/SN:NC_000072.7/SN:chr6/' \
        -e 's/SN:NC_000073.7/SN:chr7/' \
        -e 's/SN:NC_000074.7/SN:chr8/' \
        -e 's/SN:NC_000075.7/SN:chr9/' \
        -e 's/SN:NC_000076.7/SN:chr10/' \
        -e 's/SN:NC_000077.7/SN:chr11/' \
        -e 's/SN:NC_000078.7/SN:chr12/' \
        -e 's/SN:NC_000079.7/SN:chr13/' \
        -e 's/SN:NC_000080.7/SN:chr14/' \
        -e 's/SN:NC_000081.7/SN:chr15/' \
        -e 's/SN:NC_000082.7/SN:chr16/' \
        -e 's/SN:NC_000083.7/SN:chr17/' \
        -e 's/SN:NC_000084.7/SN:chr18/' \
        -e 's/SN:NC_000085.7/SN:chr19/' \
        -e 's/SN:NC_000086.8/SN:chrX/' \
        -e 's/SN:NC_000087.8/SN:chrY/' | \
    samtools reheader - "$bam_file" > "$output_file"
    
    echo "Processed $bam_file -> $output_file"
done

echo "All BAM files have been processed."


