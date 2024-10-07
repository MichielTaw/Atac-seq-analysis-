#!/bin/bash

# Define full paths to your original checksum files
original_md5="/collignon/Tawdarous/atac_seq_lab/data/fastq/md5.txt"  # Path to the unsorted original MD5 file
calculated_md5="/collignon/Tawdarous/atac_seq_lab/results/md5checksum/md5checksum.txt"  # Path to the unsorted calculated MD5 file

# Paths where the sorted checksum files will be stored (overwriting existing files)
calculated_md5_sorted="/collignon/Tawdarous/atac_seq_lab/results/md5checksum/sortedMD5/calculated_md5_sorted.txt"
original_md5_sorted="/collignon/Tawdarous/atac_seq_lab/results/md5checksum/sortedMD5/original_md5_sorted.txt"

# Path for the diff output
diff_output="/collignon/Tawdarous/atac_seq_lab/results/md5checksum/sortedMD5/diff_output.txt"

# Ensure the original checksum files exist
if [ ! -f "$original_md5" ] || [ ! -f "$calculated_md5" ]; then
  echo "One or both of the original MD5 files do not exist. Please check the paths."
  exit 1
fi

# Sort and store the checksum files, overwriting any existing sorted files
sort -k 2 "$calculated_md5" | awk '{print $1}' > "$calculated_md5_sorted"
sort -k 2 "$original_md5" | awk '{print $1}' > "$original_md5_sorted"

# Use diff to compare the sorted files and output the differences
diff "$original_md5_sorted" "$calculated_md5_sorted" > "$diff_output"

# Output results
if [ -s "$diff_output" ]; then
    echo "Differences found. See '$diff_output' for details."
else
    echo "No differences found."
fi
