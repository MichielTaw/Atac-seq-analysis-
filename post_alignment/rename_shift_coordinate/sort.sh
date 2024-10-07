BAM_FOLDER="/collignon/Tawdarous/atac_seq_lab/results/post_alignment/shift_chr_renamed"

# Loop through all .bam files in the specified folder
for bam_file in "$BAM_FOLDER"/*.bam
do
    if [ -f "$bam_file" ]; then
        echo "Processing $bam_file"
        
        # Sort the BAM file (regardless of the current sorting status)
        echo "Sorting $bam_file"
        sorted_bam="${bam_file%.bam}.sorted.bam"
        samtools sort -o "$sorted_bam" "$bam_file"
        mv "$sorted_bam" "$bam_file"
        
        # Index the BAM file
        echo "Indexing $bam_file"
        samtools index "$bam_file"
    fi
done
