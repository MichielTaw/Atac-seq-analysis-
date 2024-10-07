#______________________CHiPseeker______________________


#________________IMPORT LIBS________________
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm39.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm39.knownGene
library(clusterProfiler)
library(TxDb.Mmusculus.UCSC.mm39.knownGene)
library(GenomicRanges)
library(ReactomePA)
library(org.Mm.eg.db)
library(ArchR)

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

#_________tf analysis---------
#________t06_BMP_WNT_consesus_peak
motifPositions <- getPositions(t06_BMP_WNT_consesus_peak)

motifs <- c("STAT3", "SOX2", "TCF4", "FoxA2")
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
#markerMotifs <- markerMotifs[markerMotifs %ni% "SREBF1_22"]
markerMotifs
seFoot <- getFootprints(
  ArchRProj = projHeme5, 
  positions = motifPositions[markerMotifs], 
  groupBy = "Clusters2"
)


#---------------------------------------------------TAGMATRIX---------------------------------------
txdb <- TxDb.Mmusculus.UCSC.mm39.knownGene
promoter <- getPromoters(TxDb=txdb, upstream=1000, downstream=1000)


#-------------------t24_BMP_WNT_06h_L_consensus_peak VS t06_WNT---------------------
files <- list(
  t06_BMP_WNT_consesus_peak = t06_BMP_WNT_consesus_peak,
  t06_WNT_consensus_peak = t06_WNT_consensus_peak
)

# Generate the tag matrix list
tagMatrixList <- lapply(files, getTagMatrix, windows=promoter)

# Plot the average profile with appropriate names and x-coordinates
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000), xlab="Genomic Region (5'->3')", ylab="Read Count Frequency")




#---------------t06_BMP_WNT_consesus_peak vs t24_BMP_WNT_06h_L_consensus_peak----------

# Ensure the files are wrapped in a list before passing to lapply and set names
files <- list(
  t06_BMP_WNT_consesus_peak = t06_BMP_WNT_consesus_peak,
  t24_BMP_WNT_06h_L_consensus_peak = t24_BMP_WNT_06h_L_consensus_peak
)

# Generate the tag matrix list
tagMatrixList <- lapply(files, getTagMatrix, windows=promoter)

# Plot the average profile with appropriate names and x-coordinates
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000), xlab="Genomic Region (5'->3')", ylab="Read Count Frequency")



#---------------t24_WNT_06h_L VS t24_BMP.WNT_06h_LX1----------------

# Ensure the files are wrapped in a list before passing to lapply and set names
files <- list(
  t24_WNT_06h_L_consensus_peak = t24_WNT_06h_L_consensus_peak,
  t24_BMP_WNT_06h_LX1_consensus_peak = t24_BMP_WNT_06h_LX1_consensus_peak
)

# Generate the tag matrix list
tagMatrixList <- lapply(files, getTagMatrix, windows=promoter)

# Plot the average profile with appropriate names and x-coordinates
plotAvgProf(tagMatrixList, xlim=c(-300, 300), xlab="Genomic Region (5'->3')", ylab="Read Count Frequency")



#_________________________t06_BMP_WNT_consesus_peak_______________________
# Generate the tagMatrix
tagMatrix_t06_BMP_WNT_consesus_peak <- getTagMatrix(t06_BMP_WNT_consesus_peak, windows=promoter)
tagHeatmap(tagMatrix_t06_BMP_WNT_consesus_peak, title = "tagMatrix_t06_BMP_WNT_consesus_peak")
plotAvgProf(tagMatrix_t06_BMP_WNT_consesus_peak, xlim=c(-1000, 1000),
            xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")


#_________________________t06_WNT_consensus_peak_______________________
# Generate the tagMatrix
tagMatrix_t06_WNT_consensus_peak <- getTagMatrix(t06_WNT_consensus_peak, windows=promoter)
tagHeatmap(tagMatrix_t06_WNT_consensus_peak, title = "tagMatrix_t06_WNT_consensus_peak")
plotAvgProf(tagMatrix_t06_WNT_consensus_peak, xlim=c(-1000, 1000),
            xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")


#_________________________t24_BMP_WNT_06h_LI_consensus_peak_______________________
# Generate the tagMatrix
tagMatrix_t24_BMP_WNT_06h_LI_consensus_peak <- getTagMatrix(t24_BMP_WNT_06h_LI_consensus_peak, windows=promoter)
tagHeatmap(tagMatrix_t24_BMP_WNT_06h_LI_consensus_peak, title = "tagMatrix_t24_BMP_WNT_06h_LI_consensus_peak")
plotAvgProf(tagMatrix_t24_BMP_WNT_06h_LI_consensus_peak, xlim=c(-1000, 1000),
            xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")


#_________________________t24_BMP_WNT_06h_L_consensus_peak_______________________
tagMatrix_t24_BMP_WNT_06h_L_consensus_peak <- getTagMatrix(t24_BMP_WNT_06h_L_consensus_peak , windows=promoter)
tagHeatmap(tagMatrix_t24_BMP_WNT_06h_L_consensus_peak, title = "tagMatrix_t24_BMP_WNT_06h_L_consensus_peak")
plotAvgProf(tagMatrix_t24_BMP_WNT_06h_L_consensus_peak, xlim=c(-1000, 1000),
            xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")



#_________________________t24_BMP_WNT_06h_LX1_consensus_peak_______________________
tagMatrix_t24_BMP_WNT_06h_LX1_consensus_peak <- getTagMatrix(t24_BMP_WNT_06h_LX1_consensus_peak , windows=promoter)
tagHeatmap(tagMatrix_t24_BMP_WNT_06h_LX1_consensus_peak, title = "tagMatrix_t24_BMP_WNT_06h_LX1_consensus_peak")
plotAvgProf(tagMatrix_t24_BMP_WNT_06h_LX1_consensus_peak, xlim=c(-1000, 1000),
            xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")


#_________________________t24_WNT_06h_L_consensus_peak_______________________

tagMatrix_t24_WNT_06h_L_consensus_peak <- getTagMatrix(t24_WNT_06h_L_consensus_peak , windows=promoter)
tagHeatmap(tagMatrix_t24_WNT_06h_L_consensus_peak, title = "tagMatrix_t24_WNT_06h_L_consensus_peak")
plotAvgProf(tagMatrix_t24_WNT_06h_L_consensus_peak, xlim=c(-1000, 1000),
            xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")


#________________PEAK ANNOTATION________________


#______________________t06_BMP_WNT_consesus_peak_annotate_________________

t06_BMP_WNT_consesus_peak_annotate <- annotatePeak(t06_BMP_WNT_consesus_peak, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Mm.eg.db")
plotAnnoPie(t06_BMP_WNT_consesus_peak_annotate)
plotAnnoBar(t06_BMP_WNT_consesus_peak_annotate)
plotDistToTSS(t06_BMP_WNT_consesus_peak_annotate,
              title="Distribution of transcription factor-binding loci\nrelative to TSS t06_BMP_WNT_consesus_peak_annotate")




# Let's assume 'peakAnno' has a column 'geneId' with Entrez IDs
gene_ids_t06_BMP_WNT_consesus_peak_annotate <- as.data.frame(t06_BMP_WNT_consesus_peak_annotate)$geneId
# Perform pathway enrichment analysis specifying mouse as the organism
pathway1_t06_BMP_WNT_consesus_peak_annotate <- enrichPathway(gene = gene_ids_t06_BMP_WNT_consesus_peak_annotate, organism = "mouse")

# View the top results
head(pathway1_t06_BMP_WNT_consesus_peak_annotate, 2)
dotplot(pathway1_t06_BMP_WNT_consesus_peak_annotate)




#_________________________t06_WNT_consensus_peak_______________________

t06_WNT_consensus_peak_annotate <- annotatePeak(t06_WNT_consensus_peak, tssRegion=c(-3000, 3000),
                                                   TxDb=txdb, annoDb="org.Mm.eg.db")
plotAnnoPie(t06_WNT_consensus_peak_annotate)
plotAnnoBar(t06_WNT_consensus_peak_annotate)
plotDistToTSS(t06_WNT_consensus_peak_annotate,
              title="Distribution of transcription factor-binding loci\nrelative to TSS t06_WNT_consensus_peak_annotate")


gene_ids_t06_WNT_consensus_peak_annotate <- as.data.frame(t06_WNT_consensus_peak_annotate)$geneId
pathway1_t06_WNT_consensus_peak_annotate <- enrichPathway(gene = gene_ids_t06_WNT_consensus_peak_annotate, organism = "mouse")
head(pathway1_t06_WNT_consensus_peak_annotate, 2)
dotplot(pathway1_t06_WNT_consensus_peak_annotate)


#_________________________t24_BMP_WNT_06h_LI_consensus_peak_______________________


t24_BMP_WNT_06h_LI_consensus_peak_annotate <- annotatePeak(t24_BMP_WNT_06h_LI_consensus_peak, tssRegion=c(-3000, 3000),
                                                TxDb=txdb, annoDb="org.Mm.eg.db")
plotAnnoPie(t24_BMP_WNT_06h_LI_consensus_peak_annotate)
plotAnnoBar(t24_BMP_WNT_06h_LI_consensus_peak_annotate)
plotDistToTSS(t24_BMP_WNT_06h_LI_consensus_peak_annotate,
              title="Distribution of transcription factor-binding loci\nrelative to TSS t06_WNT_consensus_peak_annotate")


gene_ids_t24_BMP_WNT_06h_LI_consensus_peak_annotate <- as.data.frame(t24_BMP_WNT_06h_LI_consensus_peak_annotate)$geneId
pathway1_gene_ids_t24_BMP_WNT_06h_LI_consensus_peak_annotate <- enrichPathway(gene = gene_ids_t24_BMP_WNT_06h_LI_consensus_peak_annotate, organism = "mouse")
head(pathway1_gene_ids_t24_BMP_WNT_06h_LI_consensus_peak_annotate, 2)
dotplot(pathway1_gene_ids_t24_BMP_WNT_06h_LI_consensus_peak_annotate)


#_________________________t24_BMP_WNT_06h_L_consensus_peak_______________________

t24_BMP_WNT_06h_L_consensus_peak_annotate <- annotatePeak(t24_BMP_WNT_06h_L_consensus_peak, tssRegion=c(-3000, 3000),
                                                           TxDb=txdb, annoDb="org.Mm.eg.db")
plotAnnoPie(t24_BMP_WNT_06h_L_consensus_peak_annotate)
plotAnnoBar(t24_BMP_WNT_06h_L_consensus_peak_annotate)
plotDistToTSS(t24_BMP_WNT_06h_L_consensus_peak_annotate,
              title="Distribution of transcription factor-binding loci\nrelative to TSS t24_BMP_WNT_06h_L_consensus_peak_annotate")


gene_ids_t24_BMP_WNT_06h_L_consensus_peak_annotate <- as.data.frame(t24_BMP_WNT_06h_L_consensus_peak_annotate)$geneId
pathway1_gene_ids_t24_BMP_WNT_06h_L_consensus_peak_annotate <- enrichPathway(gene = gene_ids_t24_BMP_WNT_06h_L_consensus_peak_annotate, organism = "mouse")
head(pathway1_gene_ids_t24_BMP_WNT_06h_L_consensus_peak_annotate, 2)
dotplot(pathway1_gene_ids_t24_BMP_WNT_06h_L_consensus_peak_annotate)


#_________________________t24_BMP_WNT_06h_LX1_consensus_peak_______________________

t24_BMP_WNT_06h_LX1_consensus_peak_annotate <- annotatePeak(t24_BMP_WNT_06h_LX1_consensus_peak, tssRegion=c(-3000, 3000),
                                                          TxDb=txdb, annoDb="org.Mm.eg.db")
plotAnnoPie(t24_BMP_WNT_06h_LX1_consensus_peak_annotate)
plotAnnoBar(t24_BMP_WNT_06h_LX1_consensus_peak_annotate)
plotDistToTSS(t24_BMP_WNT_06h_LX1_consensus_peak_annotate,
              title="Distribution of transcription factor-binding loci\nrelative to TSS t24_BMP_WNT_06h_L_consensus_peak_annotate")


gene_ids_t24_BMP_WNT_06h_LX1_consensus_peak_annotate<- as.data.frame(t24_BMP_WNT_06h_LX1_consensus_peak_annotate)$geneId
pathway1_gene_ids_t24_BMP_WNT_06h_LX1_consensus_peak_annotate <- enrichPathway(gene = gene_ids_t24_BMP_WNT_06h_LX1_consensus_peak_annotate, organism = "mouse")
head(pathway1_gene_ids_t24_BMP_WNT_06h_LX1_consensus_peak_annotate, 2)
dotplot(pathway1_gene_ids_t24_BMP_WNT_06h_LX1_consensus_peak_annotate)


#_________________________t24_WNT_06h_L_consensus_peak_______________________


t24_WNT_06h_L_consensus_peak_annotate <- annotatePeak(t24_WNT_06h_L_consensus_peak, tssRegion=c(-3000, 3000),
                                                            TxDb=txdb, annoDb="org.Mm.eg.db")
plotAnnoPie(t24_WNT_06h_L_consensus_peak_annotate)
plotAnnoBar(t24_WNT_06h_L_consensus_peak_annotate)
plotDistToTSS(t24_WNT_06h_L_consensus_peak_annotate,
              title="Distribution of transcription factor-binding loci\nrelative to TSS t24_WNT_06h_L_consensus_peak_annotate")


gene_ids_t24_WNT_06h_L_consensus_peak_annotate<- as.data.frame(t24_WNT_06h_L_consensus_peak_annotate)$geneId
pathway1_gene_ids_t24_WNT_06h_L_consensus_peak_annotate <- enrichPathway(gene = gene_ids_t24_WNT_06h_L_consensus_peak_annotate, organism = "mouse")
head(pathway1_gene_ids_t24_WNT_06h_L_consensus_peak_annotate, 2)
dotplot(pathway1_gene_ids_t24_WNT_06h_L_consensus_peak_annotate)



txdb1 <- peaks[peaks$rep.logFC>1]
txdb2 <- peaks[peaks$rep.logFC<1]
region_list <- list(geneX = txdb1, geneY = txdb2)



peakHeatmap_multiple_Sets(peak = t24_BMP_WNT_06h_L_consensus_peak,
                          upstream = 1000,downstream = 1000,
                          by = c("geneX","geneY"),
                          type = "start_site",
                          TxDb = region_list,nbin = 800)



gene <- seq2gene(t24_WNT_06h_L_consensus_peak, tssRegion = c(-3000, 3000), flankDistance = 3000, TxDb=txdb)
pathway2 <- enrichPathway(gene, organism = "mouse")
head(pathway2, 2)
dotplot(pathway2)


t06_BMP_WNT_consesus_peak_annotate <- annotatePeak(t06_BMP_WNT_consesus_peak, tssRegion=c(-3000, 3000),
                                                      TxDb=txdb, annoDb="org.Mm.eg.db")
plotAnnoBar(t06_BMP_WNT_consesus_peak_annotate)

t24_BMP_WNT_06h_L_consensus_peak_annotate <- annotatePeak(t24_BMP_WNT_06h_L_consensus_peak, tssRegion=c(-3000, 3000),
                                                      TxDb=txdb, annoDb="org.Mm.eg.db")
plotAnnoBar(t24_BMP_WNT_06h_L_consensus_peak_annotate)


t06_WNT_consensus_peak_annotate <- annotatePeak(t06_WNT_consensus_peak, tssRegion=c(-3000, 3000),
                                                   TxDb=txdb, annoDb="org.Mm.eg.db")

plotAnnoBar(t06_WNT_consensus_peak_annotate)
