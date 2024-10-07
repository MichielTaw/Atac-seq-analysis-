#______________________CHiPseeker______________________


#________________IMPORT LIBS________________
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
library(clusterProfiler)
library(TxDb.Mmusculus.UCSC.mm10.knownGene) # Ensure you have this library installed
library(GenomicRanges)
#________________PATHS________________
dir.bam <- "/collignon/Tawdarous/atac_seq_lab/results/peak_calling/macs3"

t06_BMP_WNT_rep1 <- file.path(dir.bam, "t06_BMP+WNT_rep1_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted.bam/t06_BMP+WNT_rep1_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted.bam_peaks.broadPeak")
t06_BMP_WNT_rep2 <- file.path(dir.bam, "t06_BMP+WNT_rep2_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted.bam/t06_BMP+WNT_rep2_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted.bam_peaks.broadPeak")
t06_BMP_WNT_rep3 <- file.path(dir.bam, "t06_BMP+WNT_rep3_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted.bam/t06_BMP+WNT_rep3_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted.bam_peaks.broadPeak")
t06_WNT_rep1 <- file.path(dir.bam, "t06_WNT_rep1_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted.bam/t06_WNT_rep1_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted.bam_peaks.broadPeak")
t06_WNT_rep2 <- file.path(dir.bam, "t06_WNT_rep2_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted.bam/t06_WNT_rep2_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted.bam_peaks.broadPeak")
t06_WNT_rep3 <- file.path(dir.bam, "t06_WNT_rep3_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted.bam/t06_WNT_rep3_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted.bam_peaks.broadPeak")

t24_BMP_WNT_06h_LI_rep1 <- file.path(dir.bam, "t24_BMP+WNT-06h_LI_rep1_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted.bam/t24_BMP+WNT-06h_LI_rep1_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted.bam_peaks.broadPeak")
t24_BMP_WNT_06h_LI_rep2 <- file.path(dir.bam, "t24_BMP+WNT-06h_LI_rep2_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted.bam/t24_BMP+WNT-06h_LI_rep2_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted.bam_peaks.broadPeak")
t24_BMP_WNT_06h_LI_rep3 <- file.path(dir.bam, "t24_BMP+WNT-06h_LI_rep3_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted.bam/t24_BMP+WNT-06h_LI_rep3_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted.bam_peaks.broadPeak")

t24_BMP_WNT_06h_L_rep1 <- file.path(dir.bam, "t24_BMP+WNT-06h_L_rep1_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted.bam/t24_BMP+WNT-06h_L_rep1_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted.bam_peaks.broadPeak")
t24_BMP_WNT_06h_L_rep2 <- file.path(dir.bam, "t24_BMP+WNT-06h_L_rep2_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted.bam/t24_BMP+WNT-06h_L_rep2_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted.bam_peaks.broadPeak")
t24_BMP_WNT_06h_L_rep3 <- file.path(dir.bam, "t24_BMP+WNT-06h_L_rep3_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted.bam/t24_BMP+WNT-06h_L_rep3_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted.bam_peaks.broadPeak")

t24_BMP_WNT_06h_LX1_rep1 <- file.path(dir.bam, "t24_BMP+WNT-06h_LX1_rep1_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted.bam/t24_BMP+WNT-06h_LX1_rep1_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted.bam_peaks.broadPeak")
t24_BMP_WNT_06h_LX1_rep2 <- file.path(dir.bam, "t24_BMP+WNT-06h_LX1_rep2_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted.bam/t24_BMP+WNT-06h_LX1_rep2_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted.bam_peaks.broadPeak")
t24_BMP_WNT_06h_LX1_rep3 <- file.path(dir.bam, "t24_BMP+WNT-06h_LX1_rep3_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted.bam/t24_BMP+WNT-06h_LX1_rep3_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted.bam_peaks.broadPeak")

t24_WNT_06h_L_rep1 <- file.path(dir.bam, "t24_WNT-06h_L_rep1_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted.bam/t24_WNT-06h_L_rep1_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted.bam_peaks.broadPeak")
t24_WNT_06h_L_rep2 <- file.path(dir.bam, "t24_WNT-06h_L_rep2_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted.bam/t24_WNT-06h_L_rep2_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted.bam_peaks.broadPeak")
t24_WNT_06h_L_rep3 <- file.path(dir.bam, "t24_WNT-06h_L_rep3_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted.bam/t24_WNT-06h_L_rep3_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted.bam_peaks.broadPeak")



t06_BMP_WNT_rep1_peak <- readPeakFile(t06_BMP_WNT_rep1)
t06_BMP_WNT_rep2_peak <- readPeakFile(t06_BMP_WNT_rep2)
t06_BMP_WNT_rep3_peak <- readPeakFile(t06_BMP_WNT_rep3)

t06_WNT_rep1_peak <- readPeakFile(t06_WNT_rep1)
t06_WNT_rep2_peak <- readPeakFile(t06_WNT_rep2)
t06_WNT_rep3_peak <- readPeakFile(t06_WNT_rep3)

t24_BMP_WNT_06h_LI_rep1_peak <- readPeakFile(t24_BMP_WNT_06h_LI_rep1)
t24_BMP_WNT_06h_LI_rep2_peak <- readPeakFile(t24_BMP_WNT_06h_LI_rep2)
t24_BMP_WNT_06h_LI_rep_peak <- readPeakFile(t24_BMP_WNT_06h_LI_rep3)


t24_BMP_WNT_06h_L_rep1_peak <- readPeakFile(t24_BMP_WNT_06h_L_rep1)
t24_BMP_WNT_06h_L_rep2_peak <- readPeakFile(t24_BMP_WNT_06h_L_rep2)
t24_BMP_WNT_06h_L_rep3_peak <- readPeakFile(t24_BMP_WNT_06h_L_rep3)


t24_BMP_WNT_06h_LX1_rep1_peak <- readPeakFile(t24_BMP_WNT_06h_LX1_rep1)
t24_BMP_WNT_06h_LX1_rep2_peak <- readPeakFile(t24_BMP_WNT_06h_LX1_rep2)
t24_BMP_WNT_06h_LX1_rep3_peak <- readPeakFile(t24_BMP_WNT_06h_LX1_rep3)


t24_WNT_06h_L_rep1_peak <- readPeakFile(t24_WNT_06h_L_rep1)
t24_WNT_06h_L_rep2_peak <- readPeakFile(t24_WNT_06h_L_rep2)
t24_WNT_06h_L_rep3_peak <- readPeakFile(t24_WNT_06h_L_rep3)





sample_names <- c(
  "t06_BMP_WNT_rep1", "t06_BMP_WNT_rep2", "t06_BMP_WNT_rep3",
  "t06_WNT_rep1", "t06_WNT_rep2", "t06_WNT_rep3",
  "t24_BMP_WNT_06h_LI_rep1", "t24_BMP_WNT_06h_LI_rep2", "t24_BMP_WNT_06h_LI_rep3",
  "t24_BMP_WNT_06h_L_rep1", "t24_BMP_WNT_06h_L_rep2", "t24_BMP_WNT_06h_L_rep3",
  "t24_BMP_WNT_06h_LX1_rep1", "t24_BMP_WNT_06h_LX1_rep2", "t24_BMP_WNT_06h_LX1_rep3",
  "t24_WNT_06h_L_rep1", "t24_WNT_06h_L_rep2", "t24_WNT_06h_L_rep3"
)



#--------------
# Define the function to apply sequence level mapping, UCSC style conversion, and pruning
process_peak_file <- function(peak_file, refseq_to_ucsc, nt_nw_to_ucsc, standard_chromosomes) {
  # Apply sequence level mapping and UCSC style conversion
  current_seqlevels <- seqlevels(peak_file)
  
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
  seqlevels(peak_file) <- new_seqlevels
  
  # Set seqlevels style
  seqlevelsStyle(peak_file) <- "UCSC"
  
  # Apply keepStandardChromosomes with pruning.mode
  peak_file <- keepStandardChromosomes(peak_file, pruning.mode = "coarse")
  
  # Prune the non-standard chromosomes with the specified pruning mode
  peak_file <- keepSeqlevels(peak_file, standard_chromosomes, pruning.mode = "coarse")
  
  return(peak_file)
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

# Define the NT_ and NW_ mapping
current_seqlevels <- seqlevels(readPeakFile(t24_WNT_06h_L_rep1)) # Use one of the peak files to get current seqlevels
nt_nw_seqlevels <- current_seqlevels[grepl("^NT_|^NW_", current_seqlevels)]
nt_nw_to_ucsc <- setNames(paste0("chrUn_", nt_nw_seqlevels), nt_nw_seqlevels)

# Define the standard chromosomes
standard_chromosomes <- paste0("chr", c(1:19, "X", "Y"))

# List of peak files
peak_files <- list(
  t06_BMP_WNT_rep1, t06_BMP_WNT_rep2, t06_BMP_WNT_rep3,
  t06_WNT_rep1, t06_WNT_rep2, t06_WNT_rep3,
  t24_BMP_WNT_06h_LI_rep1, t24_BMP_WNT_06h_LI_rep2, t24_BMP_WNT_06h_LI_rep3,
  t24_BMP_WNT_06h_L_rep1, t24_BMP_WNT_06h_L_rep2, t24_BMP_WNT_06h_L_rep3,
  t24_BMP_WNT_06h_LX1_rep1, t24_BMP_WNT_06h_LX1_rep2, t24_BMP_WNT_06h_LX1_rep3,
  t24_WNT_06h_L_rep1, t24_WNT_06h_L_rep2, t24_WNT_06h_L_rep3
)

# Load and process each peak file
processed_peak_files <- lapply(peak_files, function(file_path) {
  peak_file <- readPeakFile(file_path)
  process_peak_file(peak_file, refseq_to_ucsc, nt_nw_to_ucsc, standard_chromosomes)
})

# Optional: Assign back to individual variables if needed
t06_BMP_WNT_rep1_peak <- processed_peak_files[[1]]
t06_BMP_WNT_rep2_peak <- processed_peak_files[[2]]
t06_BMP_WNT_rep3_peak <- processed_peak_files[[3]]
t06_WNT_rep1_peak <- processed_peak_files[[4]]
t06_WNT_rep2_peak <- processed_peak_files[[5]]
t06_WNT_rep3_peak <- processed_peak_files[[6]]
t24_BMP_WNT_06h_LI_rep1_peak <- processed_peak_files[[7]]
t24_BMP_WNT_06h_LI_rep2_peak <- processed_peak_files[[8]]
t24_BMP_WNT_06h_LI_rep_peak <- processed_peak_files[[9]]
t24_BMP_WNT_06h_L_rep1_peak <- processed_peak_files[[10]]
t24_BMP_WNT_06h_L_rep2_peak <- processed_peak_files[[11]]
t24_BMP_WNT_06h_L_rep3_peak <- processed_peak_files[[12]]
t24_BMP_WNT_06h_LX1_rep1_peak <- processed_peak_files[[13]]
t24_BMP_WNT_06h_LX1_rep2_peak <- processed_peak_files[[14]]
t24_BMP_WNT_06h_LX1_rep3_peak <- processed_peak_files[[15]]
t24_WNT_06h_L_rep1_peak <- processed_peak_files[[16]]
t24_WNT_06h_L_rep2_peak <- processed_peak_files[[17]]
t24_WNT_06h_L_rep3_peak <- processed_peak_files[[18]]


#----------tf---------


promoter <- getPromoters(TxDb=txdb, upstream=2000, downstream=2000)

# Calculate the tag matrix

t06_BMP_WNT_rep1_peak_tag_matrix <- getTagMatrix(t06_BMP_WNT_rep1_peak, weightCol = NULL, windows=promoter)
t06_BMP_WNT_rep2_peak_tag_matrix <- getTagMatrix(t06_BMP_WNT_rep2_peak, weightCol = NULL, windows=promoter)
t06_BMP_WNT_rep3_peak_tag_matrix <- getTagMatrix(t06_BMP_WNT_rep3_peak, weightCol = NULL, windows=promoter)

## Profile plots
plotAvgProf(t06_BMP_WNT_rep1_peak_tag_matrix, xlim=c(-2000, 2000), conf=0.95,resample=500, facet="row")
plotAvgProf(t06_BMP_WNT_rep2_peak_tag_matrix, xlim=c(-2000, 2000), conf=0.95,resample=500, facet="row")
plotAvgProf(t06_BMP_WNT_rep3_peak_tag_matrix, xlim=c(-2000, 2000), conf=0.95,resample=500, facet="row")



#______t06_WNT_rep1_peak___
t06_WNT_rep1_peak_tag_matrix <- getTagMatrix(t06_WNT_rep1_peak, weightCol = NULL, windows=promoter)
t06_WNT_rep2_peak_tag_matrix <- getTagMatrix(t06_WNT_rep2_peak, weightCol = NULL, windows=promoter)
t06_WNT_rep3_peak_tag_matrix <- getTagMatrix(t06_WNT_rep3_peak, weightCol = NULL, windows=promoter)

## Profile plots
plotAvgProf(t06_WNT_rep1_peak_tag_matrix, xlim=c(-2000, 2000), conf=0.95,resample=500, facet="row")
plotAvgProf(t06_WNT_rep2_peak_tag_matrix, xlim=c(-2000, 2000), conf=0.95,resample=500, facet="row")
plotAvgProf(t06_WNT_rep3_peak_tag_matrix, xlim=c(-2000, 2000), conf=0.95,resample=500, facet="row")

#__________t24_BMP_WNT_06h_L_rep1_peak___
t24_BMP_WNT_06h_L_rep1_peak_tag_matrix <- getTagMatrix(t24_BMP_WNT_06h_L_rep1_peak, weightCol = NULL, windows=promoter)
t24_BMP_WNT_06h_L_rep2_peak_tag_matrix <- getTagMatrix(t24_BMP_WNT_06h_L_rep2_peak, weightCol = NULL, windows=promoter)
t24_BMP_WNT_06h_L_rep3_peak_tag_matrix <- getTagMatrix(t24_BMP_WNT_06h_L_rep3_peak, weightCol = NULL, windows=promoter)

plotAvgProf(t24_BMP_WNT_06h_L_rep1_peak_tag_matrix, xlim=c(-2000, 2000), conf=0.95,resample=500, facet="row")
plotAvgProf(t24_BMP_WNT_06h_L_rep2_peak_tag_matrix, xlim=c(-2000, 2000), conf=0.95,resample=500, facet="row")
plotAvgProf(t24_BMP_WNT_06h_L_rep3_peak_tag_matrix, xlim=c(-2000, 2000), conf=0.95,resample=500, facet="row")


#_________t24_WNT_06h_L_rep1_peak________
t24_WNT_06h_L_rep1_peak_tag_matrix <- getTagMatrix(t24_WNT_06h_L_rep1_peak, weightCol = NULL, windows=promoter)
t24_WNT_06h_L_rep2_peak_tag_matrix <- getTagMatrix(t24_WNT_06h_L_rep2_peak, weightCol = NULL, windows=promoter)
t24_WNT_06h_L_rep3_peak_tag_matrix <- getTagMatrix(t24_WNT_06h_L_rep3_peak, weightCol = NULL, windows=promoter)


plotAvgProf(t24_WNT_06h_L_rep1_peak_tag_matrix, xlim=c(-2000, 2000), conf=0.95,resample=500, facet="row")
plotAvgProf(t24_WNT_06h_L_rep2_peak_tag_matrix, xlim=c(-2000, 2000), conf=0.95,resample=500, facet="row")
plotAvgProf(t24_WNT_06h_L_rep3_peak_tag_matrix, xlim=c(-2000, 2000), conf=0.95,resample=500, facet="row")

#--------------------------------------COVPLOT--------------------------------------

covplot(t06_BMP_WNT_rep1_peak, weightCol="V5", title = paste("ATAC Peaks over Chromosomes -", t06_BMP_WNT_rep1_peak))
covplot(t06_BMP_WNT_rep2_peak, weightCol="V5", title = paste("ATAC Peaks over Chromosomes -", t06_BMP_WNT_rep2_peak))
covplot(t06_BMP_WNT_rep3_peak, weightCol="V5", title = paste("ATAC Peaks over Chromosomes -", t06_BMP_WNT_rep3_peak))

covplot(t06_WNT_rep1_peak, weightCol="V5", title = paste("ATAC Peaks over Chromosomes -", t06_WNT_rep1_peak))
covplot(t06_WNT_rep2_peak, weightCol="V5", title = paste("ATAC Peaks over Chromosomes -", t06_WNT_rep2_peak))
covplot(t06_WNT_rep3_peak, weightCol="V5", title = paste("ATAC Peaks over Chromosomes -", t06_WNT_rep3_peak))

covplot(t24_BMP_WNT_06h_LI_rep1_peak, weightCol="V5", title = paste("ATAC Peaks over Chromosomes -", t24_BMP_WNT_06h_LI_rep1_peak))
covplot(t24_BMP_WNT_06h_LI_rep2_peak, weightCol="V5", title = paste("ATAC Peaks over Chromosomes -", t24_BMP_WNT_06h_LI_rep2_peak))
covplot(t24_BMP_WNT_06h_LI_rep3_peak, weightCol="V5", title = paste("ATAC Peaks over Chromosomes -", t24_BMP_WNT_06h_LI_rep3_peak))

covplot(t24_BMP_WNT_06h_L_rep1_peak, weightCol="V5", title = paste("ATAC Peaks over Chromosomes -", t24_BMP_WNT_06h_L_rep1_peak))
covplot(t24_BMP_WNT_06h_L_rep2_peak, weightCol="V5", title = paste("ATAC Peaks over Chromosomes -", t24_BMP_WNT_06h_L_rep2_peak))
covplot(t24_BMP_WNT_06h_L_rep3_peak, weightCol="V5", title = paste("ATAC Peaks over Chromosomes -", t24_BMP_WNT_06h_L_rep3_peak))

covplot(t24_BMP_WNT_06h_LX1_rep1_peak, weightCol="V5", title = paste("ATAC Peaks over Chromosomes -", t24_BMP_WNT_06h_LX1_rep1_peak))
covplot(t24_BMP_WNT_06h_LX1_rep2_peak, weightCol="V5", title = paste("ATAC Peaks over Chromosomes -", t24_BMP_WNT_06h_LX1_rep2_peak))
covplot(t24_BMP_WNT_06h_LX1_rep3_peak, weightCol="V5", title = paste("ATAC Peaks over Chromosomes -", t24_BMP_WNT_06h_LX1_rep3_peak))

covplot(t24_WNT_06h_L_rep1_peak, weightCol="V5", title = paste("ATAC Peaks over Chromosomes -", t24_WNT_06h_L_rep1_peak))
covplot(t24_WNT_06h_L_rep2_peak, weightCol="V5", title = paste("ATAC Peaks over Chromosomes -", t24_WNT_06h_L_rep2_peak))
covplot(t24_WNT_06h_L_rep3_peak, weightCol="V5", title = paste("ATAC Peaks over Chromosomes -", t24_WNT_06h_L_rep3_peak))



#--------------------------------------TagMatrix--------------------------------------



# Load the TxDb object for the mouse genome
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
promoter <- getPromoters(TxDb=txdb, upstream=1000, downstream=1000)
  
  # Generate the tagMatrix
tagMatrix_t06_BMP_WNT_rep1_peak <- getTagMatrix(t06_BMP_WNT_rep1_peak, windows=promoter)
tagHeatmap(tagMatrix_t06_BMP_WNT_rep1_peak, title = "tagMatrix_t06_BMP_WNT_rep1_peak")


plotAvgProf(tagMatrix_t06_BMP_WNT_rep1_peak, xlim=c(-1000, 1000),
            xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")





replicated_peaks <- "/collignon/Tawdarous/atac_seq_lab/results/peak_calling/macs3/replicated_peaks"
t06_BMP_WNT_rep1_pooled <- file.path(replicated_peaks, "t06_BMP+WNT.replicated_broadPeak.bed")
t06_BMP_WNT_rep1_pooled_peak <- readPeakFile(t06_BMP_WNT_rep1)
t06_BMP_WNT_rep1_pooled_peak <- process_peak_file(t06_BMP_WNT_rep1_pooled_peak, refseq_to_ucsc, nt_nw_to_ucsc, standard_chromosomes)
covplot(t06_BMP_WNT_rep1_pooled_peak, weightCol="V5", title = paste("ATAC Peaks over Chromosomes -", t06_BMP_WNT_rep1_peak))

tagMatrix_t06_BMP_WNT_rep1_pooled_peak <- getTagMatrix(t06_BMP_WNT_rep1_pooled_peak, windows=promoter)
tagHeatmap(tagMatrix_t06_BMP_WNT_rep1_pooled_peak, title = "tagMatrix_t06_BMP_WNT_rep1_pooled_peak")
plotAvgProf(tagMatrix_t06_BMP_WNT_rep1_pooled_peak, xlim=c(-1000, 1000),conf = 0.95, resample = 1000,
            xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")

