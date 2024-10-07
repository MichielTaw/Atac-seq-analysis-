library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm39.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm39.knownGene
library(clusterProfiler)
library(TxDb.Mmusculus.UCSC.mm39.knownGene)
library(GenomicRanges)
library(ReactomePA)
library(org.Mm.eg.db)
library(ArchR)
# Load gene annotations
load_gene_annotations <- function() {
  txdb <- TxDb.Mmusculus.UCSC.mm39.knownGene
  genes <- genes(txdb)
  gene_ids <- names(genes)
  gene_symbols <- mapIds(org.Mm.eg.db, 
                         keys = gene_ids,
                         column = "SYMBOL",
                         keytype = "ENTREZID",
                         multiVals = "first")
  mcols(genes)$symbol <- gene_symbols
  return(genes)
}

# Load your data
genes <- load_gene_annotations()


#-----------------------Helper function------------
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
current_seqlevels <- seqlevels(readPeakFile(t06_BMP_WNT_consesus)) # Use one of the peak files to get current seqlevels
nt_nw_seqlevels <- current_seqlevels[grepl("^NT_|^NW_", current_seqlevels)]
nt_nw_to_ucsc <- setNames(paste0("chrUn_", nt_nw_seqlevels), nt_nw_seqlevels)

# Define the standard chromosomes
standard_chromosomes <- paste0("chr", c(1:19, "X", "Y"))




#________________PATHS AND PRE-PROCESS________________
consensus_peaks <- "/collignon/Tawdarous/atac_seq_lab/results/peak_calling/macs3/replicated_peaks"

t06_BMP_WNT_consesus <- file.path(consensus_peaks, "t06_BMP+WNT.replicated_broadPeak.bed")
t06_WNT_consensus <- file.path(consensus_peaks, "t06_WNT.replicated_broadPeak.bed")
t24_BMP_WNT_06h_LI_consensus <- file.path(consensus_peaks, "t24_BMP+WNT-06h_LI.replicated_broadPeak.bed")
t24_BMP_WNT_06h_L_consensus <- file.path(consensus_peaks, "t24_BMP+WNT-06h_L.replicated_broadPeak.bed")
t24_BMP_WNT_06h_LX1_consensus <- file.path(consensus_peaks, "t24_BMP+WNT-06h_LX1.replicated_broadPeak.bed")
t24_WNT_06h_L_consensus <- file.path(consensus_peaks, "t24_WNT-06h_L.replicated_broadPeak.bed")


t06_BMP_WNT_consesus_peak <- readPeakFile(t06_BMP_WNT_consesus)
t06_WNT_consensus_peak <- readPeakFile(t06_WNT_consensus)
t24_BMP_WNT_06h_LI_consensus_peak <- readPeakFile(t24_BMP_WNT_06h_LI_consensus)
t24_BMP_WNT_06h_L_consensus_peak <- readPeakFile(t24_BMP_WNT_06h_L_consensus)
t24_BMP_WNT_06h_LX1_consensus_peak <- readPeakFile(t24_BMP_WNT_06h_LX1_consensus)
t24_WNT_06h_L_consensus_peak <- readPeakFile(t24_WNT_06h_L_consensus)


t06_BMP_WNT_consesus_peak <- process_peak_file(t06_BMP_WNT_consesus_peak, refseq_to_ucsc, nt_nw_to_ucsc, standard_chromosomes)
t06_WNT_consensus_peak <- process_peak_file(t06_WNT_consensus_peak, refseq_to_ucsc, nt_nw_to_ucsc, standard_chromosomes)
t24_BMP_WNT_06h_LI_consensus_peak <- process_peak_file(t24_BMP_WNT_06h_LI_consensus_peak, refseq_to_ucsc, nt_nw_to_ucsc, standard_chromosomes)
t24_BMP_WNT_06h_L_consensus_peak <- process_peak_file(t24_BMP_WNT_06h_L_consensus_peak, refseq_to_ucsc, nt_nw_to_ucsc, standard_chromosomes)
t24_BMP_WNT_06h_LX1_consensus_peak <- process_peak_file(t24_BMP_WNT_06h_LX1_consensus_peak, refseq_to_ucsc, nt_nw_to_ucsc, standard_chromosomes)
t24_WNT_06h_L_consensus_peak <- process_peak_file(t24_WNT_06h_L_consensus_peak, refseq_to_ucsc, nt_nw_to_ucsc, standard_chromosomes)
#-------------------------------------




#---import rna-seq data---------

library(readr)
rna_seq <- read_csv("/collignon/Tawdarous/atac_seq_lab/scripts/R scripts/rna_seq_integration/PijuanSala2019-tpm.csv")
names(rna_seq)[1] <- "genes"

preprocess_rna_data <- function(rna_data) {
  # Set genes as row names
  rna_matrix <- as.matrix(rna_data[,-1])
  rownames(rna_matrix) <- rna_data$genes
  
  # Log2 transform TPM values (adding a small constant to avoid log(0))
  log2_tpm <- log2(rna_matrix + 0.01)
  
  log2_tpm
}
preprocess_rna_data <- function(rna_data) {
  # Set genes as row names
  rna_matrix <- as.matrix(rna_data[,-1])
  rownames(rna_matrix) <- rna_data$genes
  
  # Log2 transform TPM values (adding a small constant to avoid log(0))
  log2_tpm <- log2(rna_matrix + 0.1)
  
  log2_tpm
}
rna_processed <- preprocess_rna_data(rna_seq)


#-------------Find nearby peaks for each gene------------
library(GenomicRanges)

find_nearby_peaks <- function(genes, atac_gr, window = 20000) {
  # Define standard chromosomes
  standard_chromosomes <- paste0("chr", c(1:19, "X", "Y"))
  
  # Filter and trim genes
  genes <- genes[seqnames(genes) %in% standard_chromosomes]
  genes <- trim(genes)
  
  # Filter and trim ATAC-seq peaks
  atac_gr <- atac_gr[seqnames(atac_gr) %in% standard_chromosomes]
  atac_gr <- trim(atac_gr)
  
  # Remove any remaining out-of-bound ranges
  seqlengths_atac <- seqlengths(atac_gr)
  in_bounds <- start(atac_gr) >= 1 & end(atac_gr) <= seqlengths_atac[as.character(seqnames(atac_gr))]
  atac_gr <- atac_gr[in_bounds]
  
  # Ensure atac_gr has unique identifiers
  if (is.null(names(atac_gr))) {
    names(atac_gr) <- paste0("peak_", seq_along(atac_gr))
  }
  
  # Create promoter regions
  promoters <- promoters(genes, upstream = window, downstream = window)
  
  # Find overlaps
  nearby <- findOverlaps(promoters, atac_gr, ignore.strand = TRUE)
  
  if (length(nearby) == 0) {
    warning("No nearby peaks found for any genes.")
    return(data.frame(gene = character(), peak = character()))
  }
  
  # Create data frame of results
  result <- data.frame(
    gene = names(genes)[queryHits(nearby)],
    peak = names(atac_gr)[subjectHits(nearby)]
  )
  
  return(result)
}

standard_chromosomes <- paste0("chr", c(1:19, "X", "Y"))
seqlengths(t06_BMP_WNT_consesus_peak) <- seqlengths(genes)[standard_chromosomes]

# Print diagnostic information
print("Before trimming:")
print(summary(t06_BMP_WNT_consesus_peak))

# Find nearby peaks
peak_gene_pairs <- find_nearby_peaks(genes, t06_BMP_WNT_consesus_peak)

print(head(peak_gene_pairs))
print(nrow(peak_gene_pairs))

# Calculate correlations
calculate_correlations <- function(rna_data, atac_data, peak_gene_pairs) {
  correlations <- mclapply(1:nrow(peak_gene_pairs), function(i) {
    gene <- peak_gene_pairs$gene[i]
    peak <- peak_gene_pairs$peak[i]
    
    # Check if the gene exists in rna_data
    if (!(gene %in% rownames(rna_data))) {
      return(NA)
    }
    
    # Get the gene expression data
    gene_expr <- as.numeric(rna_data[gene,])
    
    # Get the peak index
    peak_index <- which(mcols(atac_data)$V4 == peak)
    
    if (length(peak_index) == 0) {
      return(NA)
    }
    
    # Get the peak accessibility data
    peak_access <- as.numeric(mcols(atac_data)$V5[peak_index])
    
    # If peak_access is a single value, repeat it to match the length of gene_expr
    if (length(peak_access) == 1) {
      peak_access <- rep(peak_access, length(gene_expr))
    }
    
    # Check if lengths match
    if (length(gene_expr) != length(peak_access)) {
      warning(paste("Length mismatch for gene", gene, "and peak", peak))
      return(NA)
    }
    
    # Calculate correlation
    tryCatch({
      cor(gene_expr, peak_access, method = "pearson")
    }, error = function(e) {
      warning(paste("Error calculating correlation for gene", gene, "and peak", peak, ":", e$message))
      NA
    })
  }, mc.cores = detectCores() - 1)
  
  unlist(correlations)
}

correlations <- calculate_correlations(rna_processed, t06_BMP_WNT_consesus_peak, peak_gene_pairs)

valid_correlations <- correlations[!is.na(correlations)]
valid_correlations