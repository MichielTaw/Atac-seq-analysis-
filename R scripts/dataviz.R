# Install the necessary packages (if not already installed)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Gviz")
BiocManager::install("GenomicRanges")
BiocManager::install("Rsamtools")
# Load the required libraries
library(Gviz)
library(GenomicRanges)
library(Rsamtools)

# Set the path to the BAM file
bamfile <- "/collignon/Tawdarous/atac_seq_lab/results/post_alignment/encode/t24_WNT-06h_L_rep3_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.bam"


# Load the required libraries
library(Gviz)
library(GenomicRanges)
library(Rsamtools)

# Set the path to the BAM file
bamfile <- "/collignon/Tawdarous/atac_seq_lab/results/post_alignment/encode/t24_WNT-06h_L_rep3_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.bam"

# Check if BAM index exists, if not create it
if (!file.exists(paste0(bamfile, ".bai"))) {
  indexBam(bamfile)
}

# Get the actual chromosome names from the BAM file
bam_header <- scanBamHeader(bamfile)
actual_chromosomes <- names(bam_header[[1]]$targets)

# Print the actual chromosome names for verification
print(actual_chromosomes)

# Function to plot each chromosome
plot_chromosome <- function(chr, bamfile, genome) {
  tryCatch({
    # Get the range of data for this chromosome
    param <- ScanBamParam(which=GRanges(chr, IRanges(1, 536870912)))  # Use a large range
    bam_data <- scanBam(bamfile, param=param)[[1]]
    if (length(bam_data$pos) == 0) {
      message(paste("No data found for", chr))
      return(NULL)
    }
    range_start <- min(bam_data$pos)
    range_end <- max(bam_data$pos + bam_data$qwidth)
    
    dTrack <- DataTrack(range = bamfile, genome = genome, 
                        type = "h", chromosome = chr, name = paste("BAM -", chr))
    
    plotTracks(dTrack, from = range_start, to = range_end, chromosome = chr)
  }, error = function(e) {
    message(paste("Error plotting", chr, ":", e$message))
  })
}

# Plotting each chromosome
for (chr in actual_chromosomes) {
  plot_chromosome(chr, bamfile, "mm9")
}


# Load the required libraries
library(Gviz)
library(GenomicRanges)
library(Rsamtools)

# Set the path to the BAM file
bamfile <- "/collignon/Tawdarous/atac_seq_lab/results/post_alignment/encode/t24_WNT-06h_L_rep3_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.bam"

# Check if BAM index exists, if not create it
if (!file.exists(paste0(bamfile, ".bai"))) {
  indexBam(bamfile)
}

# Get the actual chromosome names from the BAM file
bam_header <- scanBamHeader(bamfile)
actual_chromosomes <- names(bam_header[[1]]$targets)

# Print the actual chromosome names and their lengths for verification
print("Chromosome names and lengths in the BAM file:")
print(bam_header[[1]]$targets)

# Function to plot each chromosome
plot_chromosome <- function(chr, bamfile, genome) {
  tryCatch({
    # Get the range of data for this chromosome
    param <- ScanBamParam(which=GRanges(chr, IRanges(1, bam_header[[1]]$targets[chr])))
    bam_data <- scanBam(bamfile, param=param)[[1]]
    if (length(bam_data$pos) == 0) {
      message(paste("No data found for", chr))
      return(NULL)
    }
    range_start <- min(bam_data$pos)
    range_end <- max(bam_data$pos + bam_data$qwidth)
    
    dTrack <- DataTrack(range = bamfile, genome = genome, 
                        type = "h", chromosome = chr, name = paste("BAM -", chr))
    
    plotTracks(dTrack, from = range_start, to = range_end, chromosome = chr)
  }, error = function(e) {
    message(paste("Error plotting", chr, ":", e$message))
  })
}

# Plotting each chromosome
for (chr in actual_chromosomes) {
  plot_chromosome(chr, bamfile, "mm10")  # Changed to mm10, adjust if necessary
}

# Print some statistics about the BAM file
print("BAM file statistics:")
print(idxstatsBam(bamfile))

