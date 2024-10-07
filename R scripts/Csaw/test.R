library(GenomicAlignments)
library(GenomicRanges)
library(Rsamtools)

# Define the function to apply sequence level mapping, UCSC style conversion, and pruning
process_bam_files <- function(bam.files, refseq_to_ucsc) {
  # Function to process a single BAM file
  process_bam_file <- function(bamfile, refseq_to_ucsc) {
    # Scan BAM file header to get sequence levels
    bam_header <- scanBamHeader(bamfile)
    current_seqlevels <- names(bam_header[[1]]$targets)
    
    # Read BAM file using readGAlignments
    gal <- readGAlignments(bamfile)
    
    # Define the NT_ and NW_ mapping
    nt_nw_seqlevels <- current_seqlevels[grepl("^NT_|^NW_", current_seqlevels)]
    nt_nw_to_ucsc <- setNames(paste0("chrUn_", nt_nw_seqlevels), nt_nw_seqlevels)
    
    # Combine mappings
    full_mapping <- c(refseq_to_ucsc, nt_nw_to_ucsc)
    
    # For any remaining seqlevels, keep the original name
    for (seqlevel in current_seqlevels) {
      if (!(seqlevel %in% names(full_mapping))) {
        full_mapping[seqlevel] <- seqlevel
      }
    }
    
    # Apply the mapping
    new_seqlevels <- full_mapping[current_seqlevels]
    seqlevels(gal) <- new_seqlevels
    
    # Set seqlevels style
    seqlevelsStyle(gal) <- "UCSC"
    
    # Define the standard chromosomes
    standard_chromosomes <- paste0("chr", c(1:19, "X", "Y"))
    
    # Apply keepStandardChromosomes with pruning.mode
    gal <- keepStandardChromosomes(gal, pruning.mode = "coarse")
    
    # Prune the non-standard chromosomes with the specified pruning mode
    gal <- keepSeqlevels(gal, standard_chromosomes, pruning.mode = "coarse")
    
    return(gal)
  }
  
  # Process each BAM file
  gals <- lapply(bam.files, process_bam_file, refseq_to_ucsc)
  
  return(gals)
}

# Define the reference sequence to UCSC mapping
refseq_to_ucsc <- c(
  "NC_000067.7" = "chr1", "NC_000068.8" = "chr2", "NC_000069.7" = "chr3",
  "NC_000070.7" = "chr4", "NC_000071.7" = "chr5", "NC_000072.7" = "chr6",
  "NC_000073.7" = "chr7", "NC_000074.7" = "chr8", "NC_000075.7" = "chr9",
  "NC_000076.7" = "chr10", "NC_000077.7" = "chr11", "NC_000078.7" = "chr12",
  "NC_000079.7" = "chr13", "NC_000080.7" = "chr14", "NC_000081.7" = "chr15",
  "NC_000082.7" = "chr16", "NC_000083.7" = "chr17", "NC_000084.7" = "chr18",
  "NC_000085.7" = "chr19", "NC_000086.8" = "chrX", "NC_000087.8" = "chrY"
)
dir.bam= "/collignon/Tawdarous/atac_seq_lab/results/post_alignment/encode"

# Define BAM files
bam.files <- c(
  file.path(dir.bam, "t06_BMP+WNT_rep1_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.bam"),
  file.path(dir.bam, "t06_BMP+WNT_rep2_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.bam"),
  file.path(dir.bam, "t06_BMP+WNT_rep3_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.bam"),
  file.path(dir.bam, "t06_WNT_rep1_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.bam"),
  file.path(dir.bam, "t06_WNT_rep2_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.bam"),
  file.path(dir.bam, "t06_WNT_rep3_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.bam"),
  file.path(dir.bam, "t24_BMP+WNT-06h_LI_rep1_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.bam"),
  file.path(dir.bam, "t24_BMP+WNT-06h_LI_rep2_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.bam"),
  file.path(dir.bam, "t24_BMP+WNT-06h_LI_rep3_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.bam"),
  file.path(dir.bam, "t24_BMP+WNT-06h_L_rep1_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.bam"),
  file.path(dir.bam, "t24_BMP+WNT-06h_L_rep2_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.bam"),
  file.path(dir.bam, "t24_BMP+WNT-06h_L_rep3_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.bam"),
  file.path(dir.bam, "t24_BMP+WNT-06h_LX1_rep1_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.bam"),
  file.path(dir.bam, "t24_BMP+WNT-06h_LX1_rep2_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.bam"),
  file.path(dir.bam, "t24_BMP+WNT-06h_LX1_rep3_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.bam"),
  file.path(dir.bam, "t24_WNT-06h_L_rep1_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.bam"),
  file.path(dir.bam, "t24_WNT-06h_L_rep2_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.bam"),
  file.path(dir.bam, "t24_WNT-06h_L_rep3_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.bam")
)

# Process the BAM files
processed_bams <- process_bam_files(bam.files, refseq_to_ucsc)

# Print the processed BAMs to check the results
print(processed_bams)

output_dir <- "/collignon/Tawdarous/atac_seq_lab/results/post_alignment/encode_chrom_name"

# Save each processed BAM file
for (i in seq_along(bam.files)) {
  output_file <- file.path(output_dir, paste0(basename(bam.files[i]), ".processed.bam"))
  export(processed_bams[[i]], con = output_file, format = "BAM")
}

