#!/bin/bash

# Specify the folder containing your .bam files
BAM_FOLDER="/collignon/Tawdarous/atac_seq_lab/results/post_alignment/encode_chr_renamed"


# Loop through all .bam files in the specified folder
for bam_file in "$BAM_FOLDER"/*.bam
do
    if [ -f "$bam_file" ]; then
        echo "Indexing $bam_file"
        samtools index "$bam_file"
    fi
done

echo "Indexing complete."
