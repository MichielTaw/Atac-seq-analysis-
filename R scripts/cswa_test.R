#Tutorial link: https://nbisweden.github.io/workshop-archive/workshop-ChIP-seq/2018-11-07/labs/lab-csaw.html


#-------------LIBS
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm39.knownGene)
library(org.Mm.eg.db)
library(GenomicRanges)
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(edgeR)
library(csaw)


#________HELPER FUNCTION: 
# Define the function to apply sequence level mapping, UCSC style conversion, and pruning
process_peaks <- function(peaks, refseq_to_ucsc, nt_nw_to_ucsc, standard_chromosomes) {
  # Apply sequence level mapping and UCSC style conversion
  current_seqlevels <- seqlevels(peaks)
  
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
  seqlevels(peaks) <- new_seqlevels
  
  # Set seqlevels style
  seqlevelsStyle(peaks) <- "UCSC"
  
  # Apply keepStandardChromosomes with pruning.mode
  peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
  
  # Prune the non-standard chromosomes with the specified pruning mode
  peaks <- keepSeqlevels(peaks, standard_chromosomes, pruning.mode = "coarse")
  
  return(peaks)
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
current_seqlevels <- seqlevels(peaks)
nt_nw_seqlevels <- current_seqlevels[grepl("^NT_|^NW_", current_seqlevels)]
nt_nw_to_ucsc <- setNames(paste0("chrUn_", nt_nw_seqlevels), nt_nw_seqlevels)

# Define the standard chromosomes
standard_chromosomes <- paste0("chr", c(1:19, "X", "Y"))

# Process the peaks
processed_peaks <- process_peaks(peaks, refseq_to_ucsc, nt_nw_to_ucsc, standard_chromosomes)


#-------------------------------------PATHS-------------------------------------


dir.bam= "/collignon/Tawdarous/atac_seq_lab/results/post_alignment/encode"


t06_BMP_WNT.rep1=file.path(dir.bam,"t06_BMP+WNT_rep1_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.bam")
t06_BMP_WNT.rep2=file.path(dir.bam,"t06_BMP+WNT_rep2_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.bam")
t06_BMP_WNT.rep3=file.path(dir.bam,"t06_BMP+WNT_rep3_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.bam")

t06_WNT.rep1=file.path(dir.bam,"t06_WNT_rep1_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.bam")
t06_WNT.rep2=file.path(dir.bam,"t06_WNT_rep2_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.bam")
t06_WNT.rep3=file.path(dir.bam,"t06_WNT_rep3_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.bam")

t24_BMP_WNT_06h_LI.rep1=file.path(dir.bam,"t24_BMP+WNT-06h_LI_rep1_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.bam")
t24_BMP_WNT_06h_LI.rep2=file.path(dir.bam,"t24_BMP+WNT-06h_LI_rep2_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.bam")
t24_BMP_WNT_06h_LI.rep3=file.path(dir.bam,"t24_BMP+WNT-06h_LI_rep3_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.bam")

t24_BMP_WNT_06h_L.rep1=file.path(dir.bam,"t24_BMP+WNT-06h_L_rep1_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.bam")
t24_BMP_WNT_06h_L.rep2=file.path(dir.bam,"t24_BMP+WNT-06h_L_rep2_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.bam")
t24_BMP_WNT_06h_L.rep3=file.path(dir.bam,"t24_BMP+WNT-06h_L_rep3_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.bam")

t24_BMP_WNT_06h_LX1.rep1=file.path(dir.bam,"t24_BMP+WNT-06h_LX1_rep1_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.bam")
t24_BMP_WNT_06h_LX1.rep2=file.path(dir.bam,"t24_BMP+WNT-06h_LX1_rep2_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.bam")
t24_BMP_WNT_06h_LX1.rep3=file.path(dir.bam,"t24_BMP+WNT-06h_LX1_rep3_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.bam")

t24_WNT_06h_L_rep1=file.path(dir.bam,"t24_WNT-06h_L_rep1_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.bam")
t24_WNT_06h_L_rep2=file.path(dir.bam,"t24_WNT-06h_L_rep2_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.bam")
t24_WNT_06h_L_rep3=file.path(dir.bam,"t24_WNT-06h_L_rep3_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.bam")


#1)-------------------------------------Compare t24_BMP_WNT_06h_LI & t24_BMP_WNT_06h_L------------------------------------- 

bam.files_t24_BMP_WNT_LI_and_L <- c(t24_BMP_WNT_06h_LI.rep1,t24_BMP_WNT_06h_LI.rep2,t24_BMP_WNT_06h_LI.rep3,t24_BMP_WNT_06h_L.rep1,t24_BMP_WNT_06h_L.rep2, t24_BMP_WNT_06h_L.rep3) # nolint # nolint # nolint

grouping_t24_BMP_WNT_LI_and_L <- factor(c('t24_BMP_WNT_06h_LI', 't24_BMP_WNT_06h_LI', 't24_BMP_WNT_06h_LI','t24_BMP_WNT_06h_L', 't24_BMP_WNT_06h_L', 't24_BMP_WNT_06h_L'))
design_t24_BMP_WNT_LI_and_L <- model.matrix(~0 + grouping_t24_BMP_WNT_LI_and_L)
colnames(design_t24_BMP_WNT_LI_and_L) <- levels(grouping_t24_BMP_WNT_LI_and_L)

contrast_t24_BMP_WNT_LI_and_L <- makeContrasts(t24_BMP_WNT_06h_LI - t24_BMP_WNT_06h_L, levels=design_t24_BMP_WNT_LI_and_L)

data_t24_BMP_WNT_LI_and_L <- windowCounts(bam.files_t24_BMP_WNT_LI_and_L, ext=0, width=150,bin=TRUE) 
#############
# Modify the process_peaks function to handle RangedSummarizedExperiment objects
process_windows <- function(windows, refseq_to_ucsc, nt_nw_to_ucsc, standard_chromosomes) {
  # Extract the ranges from the RangedSummarizedExperiment object
  ranges <- rowRanges(windows)
  
  # Apply sequence level mapping and UCSC style conversion
  current_seqlevels <- seqlevels(ranges)
  
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
  seqlevels(ranges) <- new_seqlevels
  
  # Set seqlevels style
  seqlevelsStyle(ranges) <- "UCSC"
  
  # Apply keepStandardChromosomes with pruning.mode
  ranges <- keepStandardChromosomes(ranges, pruning.mode = "coarse")
  
  # Prune the non-standard chromosomes with the specified pruning mode
  ranges <- keepSeqlevels(ranges, standard_chromosomes, pruning.mode = "coarse")
  
  # Update the rowRanges of the windows object
  rowRanges(windows) <- ranges
  
  return(windows)
}

# Define the reference sequence to UCSC mapping (same as before)
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
current_seqlevels <- seqlevels(rowRanges(data_t24_BMP_WNT_LI_and_L))
nt_nw_seqlevels <- current_seqlevels[grepl("^NT_|^NW_", current_seqlevels)]
nt_nw_to_ucsc <- setNames(paste0("chrUn_", nt_nw_seqlevels), nt_nw_seqlevels)

# Define the standard chromosomes
standard_chromosomes <- paste0("chr", c(1:19, "X", "Y"))

# Process the windows
processed_windows <- process_windows(data_t24_BMP_WNT_LI_and_L, refseq_to_ucsc, nt_nw_to_ucsc, standard_chromosomes)
############
normFactorst24_BMP_WNT_LI_and_L <- calcNormFactors(data_t24_BMP_WNT_LI_and_L)
data_t24_BMP_WNT_LI_and_L$samples$norm.factors <- normFactorst24_BMP_WNT_LI_and_L
keep_t24_BMP_WNT_LI_and_ <- filterByExpr(data_t24_BMP_WNT_LI_and_L)
data_t24_BMP_WNT_LI_and_L <- data[keep_t24_BMP_WNT_LI_and_,]
saveRDS(data_t24_BMP_WNT_LI_and_L, file='/collignon/Tawdarous/atac_seq_lab/results/csaw_result/data_t24_BMP_WNT_LI_and_L.rds')
counts_t24_BMP_WNT_LI_and_L <- estimateDisp(data_t24_BMP_WNT_LI_and_L, design_t24_BMP_WNT_LI_and_L)
saveRDS(counts_t24_BMP_WNT_LI_and_L, file='/collignon/Tawdarous/atac_seq_lab/results/csaw_result/counts_t24_BMP_WNT_LI_and_L.rds')

fit_t24_BMP_WNT_LI_and_L <- glmQLFit(counts_t24_BMP_WNT_LI_and_L, design_t24_BMP_WNT_LI_and_L)
result_t24_BMP_WNT_LI_and_L <- glmQLFTest(fit_t24_BMP_WNT_LI_and_L, contrast=contrast_t24_BMP_WNT_LI_and_)
topgenes_t24_BMP_WNT_LI_and_L <- topTags(result_t24_BMP_WNT_LI_and_L, n=100)
topgenes_t24_BMP_WNT_LI_and_L <- topTags(result_t24_BMP_WNT_LI_and_L, n=Inf, p.value=0.05)
sum(topgenes_t24_BMP_WNT_LI_and_L$table$FDR < 0.05)

peaks_t24_BMP_WNT_LI_and_L <- GRanges(seqnames = topgenes_t24_BMP_WNT_LI_and_L$table$seqnames,
                 ranges = IRanges(start = topgenes_t24_BMP_WNT_LI_and_L$table$start, end = topgenes_t24_BMP_WNT_LI_and_L$table$end),
                 strand = "*",
                 logFC = topgenes_t24_BMP_WNT_LI_and_L$table$logFC,
                 FDR = topgenes_t24_BMP_WNT_LI_and_L$table$FDR)
txdb <- TxDb.Mmusculus.UCSC.mm39.knownGene
# Define the NT_ and NW_ mapping
current_seqlevels <- seqlevels(peaks_t24_BMP_WNT_LI_and_L)
nt_nw_seqlevels <- current_seqlevels[grepl("^NT_|^NW_", current_seqlevels)]
nt_nw_to_ucsc <- setNames(paste0("chrUn_", nt_nw_seqlevels), nt_nw_seqlevels)

# Define the standard chromosomes
standard_chromosomes <- paste0("chr", c(1:19, "X", "Y"))

# Process the peaks
processed_peaks_t24_BMP_WNT_LI_and_L <- process_peaks(peaks_t24_BMP_WNT_LI_and_L, refseq_to_ucsc, nt_nw_to_ucsc, standard_chromosomes)

peakAnno_t24_BMP_WNT_LI_and_L <- annotatePeak(processed_peaks_t24_BMP_WNT_LI_and_L, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Mm.eg.db")

annotated_peaks_df_t24_BMP_WNT_LI_and_L <- as.data.frame(peakAnno_t24_BMP_WNT_LI_and_L)
head(annotated_peaks_df_t24_BMP_WNT_LI_and_L)
colnames(annotated_peaks_df_t24_BMP_WNT_LI_and_L)
write.csv(annotated_peaks_df_t24_BMP_WNT_LI_and_L, file="/collignon/Tawdarous/atac_seq_lab/results/csaw_result/annotated_peaks_df_t24_BMP_WNT_LI_and_L.csv", row.names=FALSE)



###diff peaks go
peaks_increased_t24_BMP_WNT_LI_and_L <- processed_peaks_t24_BMP_WNT_LI_and_L[processed_peaks_t24_BMP_WNT_LI_and_L$logFC > 0]
peaks_decreased_t24_BMP_WNT_LI_and_L <- processed_peaks_t24_BMP_WNT_LI_and_L[processed_peaks_t24_BMP_WNT_LI_and_L$logFC < 0]


peakAnno_increased_t24_BMP_WNT_LI_and_L <- annotatePeak(peaks_increased_t24_BMP_WNT_LI_and_L, tssRegion=c(-3000, 3000),
                                   TxDb=txdb, annoDb="org.Mm.eg.db")
peakAnno_decreased_t24_BMP_WNT_LI_and_L <- annotatePeak(peaks_decreased_t24_BMP_WNT_LI_and_L, tssRegion=c(-3000, 3000),
                                   TxDb=txdb, annoDb="org.Mm.eg.db")

genes_increased_t24_BMP_WNT_LI_and_L <- unique(as.data.frame(peakAnno_increased_t24_BMP_WNT_LI_and_L)$geneId)
genes_decreased_t24_BMP_WNT_LI_and_L <- unique(as.data.frame(peakAnno_decreased_t24_BMP_WNT_LI_and_L)$geneId)

# Remove NA values
genes_increased_t24_BMP_WNT_LI_and_L <- genes_increased_t24_BMP_WNT_LI_and_L[!is.na(genes_increased_t24_BMP_WNT_LI_and_L)]
genes_decreased_t24_BMP_WNT_LI_and_L <- genes_decreased_t24_BMP_WNT_LI_and_L[!is.na(genes_decreased_t24_BMP_WNT_LI_and_L)]

##INCREASES
ego_BP_increased <- enrichGO(gene = genes_increased_t24_BMP_WNT_LI_and_L,
                             OrgDb = org.Mm.eg.db,
                             ont = "BP",
                             pAdjustMethod = "BH",
                             pvalueCutoff = 0.05,
                             qvalueCutoff = 0.05)
dotplot(ego_BP_increased, title = "ontologie: +++ dans t24_BMP_WNT_LI")

ego_CC <- enrichGO(gene = genes_increased_t24_BMP_WNT_LI_and_L,
                   OrgDb = org.Mm.eg.db,
                   ont = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05)

dotplot(ego_CC, showCategory=20)

##DECREASED
ego_BP_decreased <- enrichGO(gene = genes_decreased_t24_BMP_WNT_LI_and_L,
                             OrgDb = org.Mm.eg.db,
                             ont = "BP",
                             pAdjustMethod = "BH",
                             pvalueCutoff = 0.05,
                             qvalueCutoff = 0.05)
dotplot(ego_BP_decreased, title = "ontologie +++ dans t24_BMP_WNT_06h_L")

ego_CC <- enrichGO(gene = genes_decreased_t24_BMP_WNT_LI_and_L,
                   OrgDb = org.Mm.eg.db,
                   ont = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05)

dotplot(ego_CC, showCategory=20)

#-------------------------------------Compare t06_BMP_WNT & t06_WNT-------------------------------------

bam.files <- c(t06_BMP_WNT.rep1,t06_BMP_WNT.rep2,t06_BMP_WNT.rep3,t06_WNT.rep1,t06_WNT.rep2, t06_WNT.rep3)

grouping <- factor(c('t06_BMP_WNT', 't06_BMP_WNT', 't06_BMP_WNT','t06_WNT', 't06_WNT', 't06_WNT'))
design <- model.matrix(~0 + grouping)
colnames(design) <- levels(grouping)



library(edgeR)
contrast <- makeContrasts(t06_BMP_WNT - t06_WNT, levels=design)


#####################CSAW####################
data <- windowCounts(bam.files, ext=100, width=100) 
normFactors <- calcNormFactors(data)
data$samples$norm.factors <- normFactors
keep <- filterByExpr(data)
data <- data[keep,]
counts <- estimateDisp(data, design)
fit <- glmQLFit(counts, design)
result <- glmQLFTest(fit, contrast=contrast)
topgenes <- topTags(result, n=100)
topgenes <- topTags(result, n=Inf, p.value=0.05)
sum(topgenes$table$FDR < 0.05)

#_______annotate peaks:
peaks <- GRanges(seqnames = topgenes$table$seqnames,
                 ranges = IRanges(start = topgenes$table$start, end = topgenes$table$end),
                 strand = "*",
                 logFC = topgenes$table$logFC,
                 FDR = topgenes$table$FDR)

txdb <- TxDb.Mmusculus.UCSC.mm39.knownGene

peakAnno <- annotatePeak(processed_peaks, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Mm.eg.db")

annotated_peaks_df <- as.data.frame(peakAnno)
head(annotated_peaks_df)
colnames(annotated_peaks_df)
write.csv(annotated_peaks_df, file="/collignon/Tawdarous/atac_seq_lab/scripts/R scripts/Csaw/test/annotated_peaks.csv", row.names=FALSE)
#-------------gene enrichment

genes <- unique(as.data.frame(peakAnno)$geneId)
genes <- genes[!is.na(genes)]  # Remove NA values
#Biological Process
ego_BP <- enrichGO(gene = genes,
                   OrgDb = org.Mm.eg.db,
                   ont = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05)
dotplot(ego_BP, showCategory=20)



# Molecular Function
ego_MF <- enrichGO(gene = genes,
                   OrgDb = org.Mm.eg.db,
                   ont = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05)
dotplot(ego_MF, showCategory=20)


# Cellular Component
ego_CC <- enrichGO(gene = genes,
                   OrgDb = org.Mm.eg.db,
                   ont = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05)

dotplot(ego_CC, showCategory=20)



#____KEGG

kk <- enrichKEGG(gene = genes,
                 organism = 'mmu',
                 pvalueCutoff = 0.05)

# Visualize KEGG results
dotplot(kk, showCategory=20)


###diff peaks go
peaks_increased <- processed_peaks[processed_peaks$logFC > 0]
peaks_decreased <- processed_peaks[processed_peaks$logFC < 0]


peakAnno_increased <- annotatePeak(peaks_increased, tssRegion=c(-3000, 3000),
                                   TxDb=txdb, annoDb="org.Mm.eg.db")
peakAnno_decreased <- annotatePeak(peaks_decreased, tssRegion=c(-3000, 3000),
                                   TxDb=txdb, annoDb="org.Mm.eg.db")

genes_increased <- unique(as.data.frame(peakAnno_increased)$geneId)
genes_decreased <- unique(as.data.frame(peakAnno_decreased)$geneId)

# Remove NA values
genes_increased <- genes_increased[!is.na(genes_increased)]
genes_decreased <- genes_decreased[!is.na(genes_decreased)]

##INCREASES
ego_BP_increased <- enrichGO(gene = genes_increased,
                             OrgDb = org.Mm.eg.db,
                             ont = "BP",
                             pAdjustMethod = "BH",
                             pvalueCutoff = 0.05,
                             qvalueCutoff = 0.05)
dotplot(ego_BP_increased, title = "ontologie: +++ dans T06_BMP_WNT")

ego_CC <- enrichGO(gene = genes_increased,
                   OrgDb = org.Mm.eg.db,
                   ont = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05)

dotplot(ego_CC, showCategory=20)

##DECREASED
ego_BP_decreased <- enrichGO(gene = genes_decreased,
                             OrgDb = org.Mm.eg.db,
                             ont = "BP",
                             pAdjustMethod = "BH",
                             pvalueCutoff = 0.05,
                             qvalueCutoff = 0.05)
dotplot(ego_BP_decreased, title = "ontologie: +++ dans T06_WNT")

ego_CC <- enrichGO(gene = genes_decreased,
                   OrgDb = org.Mm.eg.db,
                   ont = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05)

dotplot(ego_CC, showCategory=20)

#_____________


