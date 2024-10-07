#!/bin/bash

idx_path="/collignon/Tawdarous/atac_seq_lab/results/alignment/bam_files"
mkdir -p /collignon/Tawdarous/atac_seq_lab/results/alignment/bam_files/idx_report

for idx_file in "$idx_path"/*idxstats; do
    mv "$idx_file"  /collignon/Tawdarous/data/Tosic19/atac_seq_PRJNA369204/results/alignment/post_alignment_quality_control/idx_report
done
