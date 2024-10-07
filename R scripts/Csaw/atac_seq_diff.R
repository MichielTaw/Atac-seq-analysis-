#https://ivanek.github.io/analysisOfGenomicsDataWithR/12_ATACSeq_html.html#differential-accessibility-analysis
#---------------------------------LIBRARY---------------------------------:
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm39.knownGene)
library(BSgenome.Mmusculus.UCSC.mm39)
library(org.Mm.eg.db)
library(GenomicRanges)
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(edgeR)
library(csaw)
library(RColorBrewer)
library(tidyverse)
library(ggrepel)
library(ggrastr)
options(ggrastr.default.dpi=150)
library(patchwork)
library(SummarizedExperiment)
library(rtracklayer)  
library(Gviz)   
library(GenomicAlignments)
library(GenomicRanges)



#---------------------------------HELPER FUNCTION---------------------------------: 
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
##############



#############
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


bam.files<- c(t06_BMP_WNT.rep1,t06_BMP_WNT.rep2,t06_BMP_WNT.rep3,
              t06_WNT.rep1,t06_WNT.rep2,t06_WNT.rep3,
              t24_BMP_WNT_06h_LI.rep1,t24_BMP_WNT_06h_LI.rep2,t24_BMP_WNT_06h_LI.rep3,
              t24_BMP_WNT_06h_L.rep1,t24_BMP_WNT_06h_L.rep2,t24_BMP_WNT_06h_L.rep3,
              t24_BMP_WNT_06h_LX1.rep1,t24_BMP_WNT_06h_LX1.rep2,t24_BMP_WNT_06h_LX1.rep3,
              t24_WNT_06h_L_rep1,t24_WNT_06h_L_rep2,t24_WNT_06h_L_rep3)

sample_names <- basename(bam.files)
cell_types <- sub("_rep\\d+.*", "", sample_names)

#-----------------------------------Windows count, Width=150-------------------------

param <- readParam(pe="none", dedup=TRUE)
windows <- windowCounts(bam.files, ext=0, width=150,
                        bin=TRUE, filter=0, param=param)
saveRDS(windows, file='/collignon/Tawdarous/atac_seq_lab/results/csaw_result/windows.rds')

#-----------------------------------Process Windows-------------------------

current_seqlevels <- seqlevels(rowRanges(windows))
nt_nw_seqlevels <- current_seqlevels[grepl("^NT_|^NW_", current_seqlevels)]
nt_nw_to_ucsc <- setNames(paste0("chrUn_", nt_nw_seqlevels), nt_nw_seqlevels)
standard_chromosomes <- paste0("chr", c(1:19, "X", "Y"))
processed_peak <- process_peaks(windows, refseq_to_ucsc, nt_nw_to_ucsc, standard_chromosomes)
colData(processed_peak)$CellType <- cell_types
##

background <- windowCounts(bam.files, bin=TRUE, width=10000L, param=param)
saveRDS(background, file='/collignon/Tawdarous/atac_seq_lab/results/csaw_result/background.rds')
processed_backgroud <- process_peaks(background, refseq_to_ucsc, nt_nw_to_ucsc, standard_chromosomes)
saveRDS(processed_backgroud, file='/collignon/Tawdarous/atac_seq_lab/results/csaw_result/processed_backgroud.rds')





#-----------------------------------MA plot-------------------------


# Define transcription start sites (TSS)
tss <- genes(TxDb.Mmusculus.UCSC.mm39.knownGene)
tss <- resize(tss, width = 1, fix = "start")
overlapTSS <- rowRanges(processed_peak) %over% tss
table(overlapTSS)

colData <- DataFrame(
  bam.files = bam.files,
  sample_names = sample_names,
  cell_types = cell_types
)
colData(processed_peak) <- colData
cell_types <- sub("_rep\\d+.*", "", sample_names)
colData(processed_peak)$CellType <- cell_types

colData(processed_peak)$totals <- colData(filtered.processed_peak)$totals
logCPM <-  csaw::calculateCPM(processed_peak, log=TRUE)

colnames(logCPM) <- sample_names
colnames(filtered.processed_peak) <- sample_names
colnames(logCPM) <- colnames(filtered.processed_peak)

ct <- split(colnames(logCPM), colData(processed_peak)$CellType)
#ct

df <- vapply(ct, function(x) {
  rowMeans(logCPM[,x,drop=FALSE])
}, rep(0, nrow(processed_peak))) |>
  as_tibble() |>
  mutate(overlapTSS = rowData(processed_peak)$overlapTSS,
         GC = rowData(processed_peak)[,'G|C'])
plotMA <- function(df, s1, s2) {
  ggplot(df,
         mapping = aes(x = .data[[s1]] + .data[[s2]],
                       y = .data[[s1]] - .data[[s2]])) +
    geom_hex(bins = 100, aes(fill=after_stat(density)^(1/16)), show.legend = FALSE) +
    rasterize( geom_point(data = subset(df, overlapTSS), col='red', pch='.') ) +
    rasterize( geom_smooth(data = subset(df, overlapTSS), col='plum', se = FALSE) ) +
    geom_hline(yintercept = 0, lty = 2) +
    labs(x = "A", y="M", title=sprintf("%s vs. %s", s1, s2)) +
    theme_bw(base_size = 15)
}

combi <- list(c('t06_BMP+WNT','t06_WNT'), c('t24_BMP+WNT-06h_L','t24_BMP+WNT-06h_LI'), c('t24_BMP+WNT-06h_LX1','t24_WNT-06h_L'))
gl <- lapply(combi, function(x) plotMA(df, x[1], x[2]))
patchwork::wrap_plots(gl, ncol=2)


combi2 <- list(c('t06_BMP+WNT','t06_WNT'), c('t06_BMP+WNT','t24_BMP+WNT-06h_L'), c('t06_WNT','t24_WNT-06h_L'))
gl2 <- lapply(combi2, function(x) plotMA(df, x[1], x[2]))
patchwork::wrap_plots(gl2, ncol=2)


combi2 <- list(c('t24_BMP+WNT-06h_L','t24_BMP+WNT-06h_LI'), c('t24_BMP+WNT-06h_L','t24_BMP+WNT-06h_LX1'), c('t24_BMP+WNT-06h_LI','t24_BMP+WNT-06h_LX1'), c('t24_BMP+WNT-06h_L','t24_WNT-06h_L'))
gl2 <- lapply(combi2, function(x) plotMA(df, x[1], x[2]))
patchwork::wrap_plots(gl2, ncol=2)


combi2 <- list(c('t06_BMP+WNT','t24_WNT-06h_L'), c('t06_BMP+WNT','t24_BMP+WNT-06h_LX1'), c('t06_BMP+WNT','t24_BMP+WNT-06h_LI'), c('t06_BMP+WNT','t24_BMP+WNT-06h_L'))
gl2 <- lapply(combi2, function(x) plotMA(df, x[1], x[2]))
patchwork::wrap_plots(gl2, ncol=2)


combi2 <- list(
  c('t06_WNT','t24_WNT-06h_L'),
  c('t06_WNT','t24_BMP+WNT-06h_LX1'), 
  c('t06_WNT','t24_BMP+WNT-06h_LI'), 
  c('t06_WNT','t24_BMP+WNT-06h_L'))
gl2 <- lapply(combi2, function(x) plotMA(df, x[1], x[2]))
patchwork::wrap_plots(gl2, ncol=2)

#-----------------------------------Peak calling 1-------------------------

#import again processed_peak and processed_backgroud
filterStat <- filterWindowsGlobal(processed_peak, processed_backgroud)
saveRDS(filterStat, file='/collignon/Tawdarous/atac_seq_lab/results/csaw_result/filterStat.rds')

str(filterStat)
globalBG <- median( filterStat$back.abundances )
globalBG
head(filterStat$abundances - globalBG)

keep <- filterStat$filter > log2(3)
sum(keep)
mean(keep)
min(filterStat$abundances[keep])
globalBG + log2(3)


ggplot(data = tibble(x=filterStat$abundances),
       mapping = aes(x)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = globalBG, color='red') +
  geom_vline(xintercept = globalBG + log2(3), color='blue') +
  labs(x="Abundance", y="Frequency", title="Histogram of processed Abundances") +
  theme_bw(base_size=15))

rowData(processed_peak)$enrichment <- filterStat$filter
filtered.processed_peak <- processed_peak[keep,]

rtracklayer::export(rowRanges(filtered.processed_peak),
                    '/collignon/Tawdarous/atac_seq_lab/results/csaw_result/filtered.processed_peak.bed')


filtered.processed_peak <- normOffsets(filtered.processed_peak)
filtered.processed_peak$totals <-processed_peak$totals
#filtered.processed_peak$totals
saveRDS(filtered.processed_peak, file='/collignon/Tawdarous/atac_seq_lab/results/csaw_result/filtered.processed_peak.rds')

logCPM <- csaw::calculateCPM(filtered.processed_peak, log=TRUE, use.offsets = TRUE)
colData(filtered.processed_peak)$CellType <- cell_types
colnames(logCPM) <- sample_names
ct <- split(colnames(logCPM), colData(filtered.processed_peak)$CellType)


combi2 <- list(
  c('t06_BMP+WNT','t06_WNT'), 
  c('t06_BMP+WNT','t24_BMP+WNT-06h_L'), 
  c('t06_WNT','t24_WNT-06h_L'))
#filtered.processed_peak$totals <- processed_peak$totals

df <- vapply(ct, function(x) {
  rowMeans(logCPM[,x,drop=FALSE])
}, rep(0, nrow(filtered.processed_peak))) |>
  as_tibble() |>
  mutate(overlapTSS = rowData(filtered.processed_peak)$overlapTSS,
         GC = rowData(filtered.processed_peak)[,'G|C'])

patchwork::wrap_plots(lapply(combi2, function(x) plotMA(df, x[1], x[2])), ncol=2)



#-----------------------------------Differential accessibility test-------------------------
dgel <- csaw::asDGEList(filtered.processed_peak,
                        genes = as.data.frame(rowData(filtered.processed_peak)),
                        group = colData(filtered.processed_peak)$CellType)
colnames(dgel) <- colnames(filtered.processed_peak)
rownames(dgel) <- as.character(rowRanges(filtered.processed_peak))
dgel

moma <- model.matrix(~ 0 + group, dgel$samples)
colnames(moma) <- sub('group','',colnames(moma))
moma

dgel <- estimateDisp(dgel, design=moma)
saveRDS(dgel, file='/collignon/Tawdarous/atac_seq_lab/results/csaw_result/dgel.rds')

fit <- glmQLFit(dgel, moma)
saveRDS(fit, file='/collignon/Tawdarous/atac_seq_lab/results/csaw_result/fit.rds')
saveRDS(moma, file='/collignon/Tawdarous/atac_seq_lab/results/csaw_result/moma.rds')


colnames(moma) <- gsub("-", "_", colnames(moma))
colnames(moma) <- make.names(colnames(moma), unique = TRUE)


#__________________________t24_BMP.WNT_06h_LX1 - t24_BMP.WNT_06h_L _______________________
contr <- makeContrasts(t24_BMP.WNT_06h_LX1 - t24_BMP.WNT_06h_L, levels = colnames(moma))
res <- glmQLFTest(fit, contrast = contr)
table(decideTests(res, p.value = 0.05, lfc = 1.5))
tt <- topTags(res, n=Inf, sort.by="none")$table[,c('logFC','logCPM','F','PValue','FDR')]

tt |>
  as_tibble() |>
  mutate(Windows = rownames(dgel), .before=1) |>
  arrange(FDR) |>
  head()


rowData(filtered.processed_peak) <- cbind(rowData(filtered.processed_peak), tt)
merged <- mergeWindows(rowRanges(filtered.processed_peak), tol=300L)
str(merged,1)
merged$region
summary(width(merged$region))



anno <- detailRanges(merged$region,
                     txdb=TxDb.Mmusculus.UCSC.mm39.knownGene,
                     orgdb=org.Mm.eg.db, promoter=c(1500, 500), dist=10*1e3L)

merged$region$overlap <- anno$overlap
merged$region$left <- anno$left
merged$region$right <- anno$right
merged$region
head(merged$region$overlap)

region.stats <- combineTests(merged$id, tt) |>
  as_tibble() |>
  mutate(rep.logCPM = tt[rep.test,'logCPM'],
         Contrast = colnames(contr))
region.stats |>
  head()
peaks <- merged$region 

mcols(peaks) <- c(mcols(peaks), region.stats)

peaks$name <- ifelse(peaks$rep.logFC > 1.5 & peaks$FDR < 0.05, 't24_BMP.WNT_06h_LX1',
                     ifelse(peaks$rep.logFC < (-1.5) & peaks$FDR < 0.05, 't24_BMP.WNT_06h_L', 'Common'))

table(peaks$name, peaks$direction)






data <- mcols(peaks) |>
  as_tibble() |>
  mutate(logFC = rep.logFC,
         Significant = FDR < 0.05,
         Annotated = overlap != "")


ggplot(data = data,
       mapping = aes(x = logFC,
                     y = -log10(PValue),
                     color = Significant,
                     shape = Annotated,
                     label = overlap)) +
  rasterize( geom_point() ) +
  geom_text_repel(data = data |> filter(Annotated & Significant & abs(logFC) > 1.5)) +
  geom_vline(xintercept = 0, lty = 2) +
  theme_bw()


write.table(as.data.frame(peaks), file=gzfile('/collignon/Tawdarous/atac_seq_lab/results/csaw_result/peaks_t06_BMP.WNT_VS_t24_BMP.WNT_06h_L.tsv.gz'),
            row.names=FALSE, sep='\t')
rtracklayer::export(peaks, '/collignon/Tawdarous/atac_seq_lab/results/csaw_result/peaks_t06_BMP.WNT_VS_t24_BMP.WNT_06h_L.bed')

significant_genes <- data %>%
  filter(Significant & Annotated & abs(logFC) > 4) %>%
  arrange(desc(abs(logFC)))
significant_genes %>%
  select(overlap, logFC, PValue) %>%
  head(20) %>%
  print(n = 20)

# Upregulated genes
cat("Top upregulated genes:\n")
significant_genes %>%
  filter(logFC > 0) %>%
  select(overlap, logFC, PValue) %>%
  head(10) %>%
  print(n = 10)

# Downregulated genes
cat("\nTop downregulated genes:\n")
significant_genes %>%
  filter(logFC < 0) %>%
  select(overlap, logFC, PValue) %>%
  head(10) %>%
  print(n = 10)


#__________________________t06_BMP.WNT VS t24_BMP.WNT_06h_L _______________________
contr <- makeContrasts(t06_BMP.WNT - t24_BMP.WNT_06h_L, levels = colnames(moma))
res <- glmQLFTest(fit, contrast = contr)
table(decideTests(res, p.value = 0.05, lfc = 1.5))
tt <- topTags(res, n=Inf, sort.by="none")$table[,c('logFC','logCPM','F','PValue','FDR')]

tt |>
  as_tibble() |>
  mutate(Windows = rownames(dgel), .before=1) |>
  arrange(FDR) |>
  head()


rowData(filtered.processed_peak) <- cbind(rowData(filtered.processed_peak), tt)
merged <- mergeWindows(rowRanges(filtered.processed_peak), tol=5000L)
str(merged,1)
merged$region
summary(width(merged$region))



anno <- detailRanges(merged$region,
                     txdb=TxDb.Mmusculus.UCSC.mm39.knownGene,
                     orgdb=org.Mm.eg.db, promoter=c(1500, 500), dist=10*1e3L)

merged$region$overlap <- anno$overlap
merged$region$left <- anno$left
merged$region$right <- anno$right
merged$region
head(merged$region$overlap)

region.stats <- combineTests(merged$id, tt) |>
  as_tibble() |>
  mutate(rep.logCPM = tt[rep.test,'logCPM'],
         Contrast = colnames(contr))
region.stats |>
  head()
peaks <- merged$region 

mcols(peaks) <- c(mcols(peaks), region.stats)

peaks$name <- ifelse(peaks$rep.logFC > 1.5 & peaks$FDR < 0.05, 't06_BMP.WNT',
                     ifelse(peaks$rep.logFC < (-1.5) & peaks$FDR < 0.05, 't24_BMP.WNT_06h_L', 'Common'))

table(peaks$name, peaks$direction)






data <- mcols(peaks) |>
  as_tibble() |>
  mutate(logFC = rep.logFC,
         Significant = FDR < 0.05,
         Annotated = overlap != "")


ggplot(data = data,
       mapping = aes(x = logFC,
                     y = -log10(PValue),
                     color = Significant,
                     shape = Annotated,
                     label = overlap)) +
  rasterize( geom_point() ) +
  geom_text_repel(data = data |> filter(Annotated & Significant & abs(logFC) > 1.5)) +
  geom_vline(xintercept = 0, lty = 2) +
  theme_bw()


write.table(as.data.frame(peaks), file=gzfile('/collignon/Tawdarous/atac_seq_lab/results/csaw_result/peaks_t06_BMP.WNT_VS_t24_BMP.WNT_06h_L.tsv.gz'),
            row.names=FALSE, sep='\t')
rtracklayer::export(peaks, '/collignon/Tawdarous/atac_seq_lab/results/csaw_result/peaks_t06_BMP.WNT_VS_t24_BMP.WNT_06h_L.bed')

significant_genes <- data %>%
  filter(Significant & Annotated & abs(logFC) > 4) %>%
  arrange(desc(abs(logFC)))
significant_genes %>%
  select(overlap, logFC, PValue) %>%
  head(20) %>%
  print(n = 20)

# Upregulated genes
cat("Top upregulated genes:\n")
significant_genes %>%
  filter(logFC > 0) %>%
  select(overlap, logFC, PValue) %>%
  head(10) %>%
  print(n = 10)

# Downregulated genes
cat("\nTop downregulated genes:\n")
significant_genes %>%
  filter(logFC < 0) %>%
  select(overlap, logFC, PValue) %>%
  head(10) %>%
  print(n = 10)
#__________________________t06_BMP.WNT VS t24_BMP.WNT_06h_L _______________________
contr <- makeContrasts(t06_BMP.WNT - t24_BMP.WNT_06h_L, levels = colnames(moma))
res <- glmQLFTest(fit, contrast = contr)
table(decideTests(res, p.value = 0.05, lfc = 1.5))
tt <- topTags(res, n=Inf, sort.by="none")$table[,c('logFC','logCPM','F','PValue','FDR')]

tt |>
  as_tibble() |>
  mutate(Windows = rownames(dgel), .before=1) |>
  arrange(FDR) |>
  head()


rowData(filtered.processed_peak) <- cbind(rowData(filtered.processed_peak), tt)
merged <- mergeWindows(rowRanges(filtered.processed_peak), tol=300L)
str(merged,1)
merged$region
summary(width(merged$region))



anno <- detailRanges(merged$region,
                     txdb=TxDb.Mmusculus.UCSC.mm39.knownGene,
                     orgdb=org.Mm.eg.db, promoter=c(1500, 500), dist=10*1e3L)

merged$region$overlap <- anno$overlap
merged$region$left <- anno$left
merged$region$right <- anno$right
merged$region
head(merged$region$overlap)

region.stats <- combineTests(merged$id, tt) |>
  as_tibble() |>
  mutate(rep.logCPM = tt[rep.test,'logCPM'],
         Contrast = colnames(contr))
region.stats |>
  head()
peaks <- merged$region 

mcols(peaks) <- c(mcols(peaks), region.stats)

peaks$name <- ifelse(peaks$rep.logFC > 1.5 & peaks$FDR < 0.05, 't06_BMP.WNT',
                     ifelse(peaks$rep.logFC < (-1.5) & peaks$FDR < 0.05, 't24_BMP.WNT_06h_L', 'Common'))

table(peaks$name, peaks$direction)

data <- mcols(peaks) |>
  as_tibble() |>
  mutate(logFC = rep.logFC,
         Significant = FDR < 0.05,
         Annotated = overlap != "")


ggplot(data = data,
       mapping = aes(x = logFC,
                     y = -log10(PValue),
                     color = Significant,
                     shape = Annotated,
                     label = overlap)) +
  rasterize( geom_point() ) +
  geom_text_repel(data = data |> filter(Annotated & Significant & abs(logFC) > 1.5)) +
  geom_vline(xintercept = 0, lty = 2) +
  theme_bw()


write.table(as.data.frame(peaks), file=gzfile('/collignon/Tawdarous/atac_seq_lab/results/csaw_result/peaks_t06_BMP.WNT_VS_t24_BMP.WNT_06h_L.tsv.gz'),
            row.names=FALSE, sep='\t')
rtracklayer::export(peaks, '/collignon/Tawdarous/atac_seq_lab/results/csaw_result/peaks_t06_BMP.WNT_VS_t24_BMP.WNT_06h_L.bed')

significant_genes <- data %>%
  filter(Significant & Annotated & abs(logFC) > 4) %>%
  arrange(desc(abs(logFC)))
significant_genes %>%
  select(overlap, logFC, PValue) %>%
  head(20) %>%
  print(n = 20)

# Upregulated genes
cat("Top upregulated genes:\n")
significant_genes %>%
  filter(logFC > 0) %>%
  select(overlap, logFC, PValue) %>%
  head(10) %>%
  print(n = 10)

# Downregulated genes
cat("\nTop downregulated genes:\n")
significant_genes %>%
  filter(logFC < 0) %>%
  select(overlap, logFC, PValue) %>%
  head(10) %>%
  print(n = 10)
#__________________________t24_WNT-06h_L VS t24_BMP+WNT-06h_L VS t24_BMP+WNT-06h_LI _______________________
contr <- makeContrasts(
  WNT_vs_BMP.WNT_L = t24_WNT_06h_L - t24_BMP.WNT_06h_L,
  WNT_vs_BMP.WNT_LI = t24_WNT_06h_L - t24_BMP.WNT_06h_LI,
  BMP.WNT_L_vs_LI = t24_BMP.WNT_06h_L - t24_BMP.WNT_06h_LI,
  levels = colnames(moma)
)

# Perform the test
res <- glmQLFTest(fit, contrast = contr)

# Extract results
tt <- topTags(res, n=Inf, sort.by="none")$table[,c('logFC.WNT_vs_BMP.WNT_L','logFC.WNT_vs_BMP.WNT_LI','logFC.BMP.WNT_L_vs_LI','logCPM','F','PValue','FDR')]
rowData(filtered.processed_peak) <- cbind(rowData(filtered.processed_peak), tt)

merged <- mergeWindows(rowRanges(filtered.processed_peak), tol=300L)
anno <- detailRanges(merged$region,
                     txdb=TxDb.Mmusculus.UCSC.mm39.knownGene,
                     orgdb=org.Mm.eg.db, promoter=c(1500, 500), dist=10*1e3L)

merged$region$overlap <- anno$overlap
merged$region$left <- anno$left
merged$region$right <- anno$right

region.stats <- combineTests(merged$id, tt)
peaks <- merged$region 
mcols(peaks) <- c(mcols(peaks), region.stats)


peaks$name <- ifelse(peaks$rep.logFC.WNT_vs_BMP.WNT_L > 0.5 & peaks$FDR < 0.05, 'WNT_vs_BMP.WNT_L_up',
                     ifelse(peaks$rep.logFC.WNT_vs_BMP.WNT_L < -0.5 & peaks$FDR < 0.05, 'BMP.WNT_L_vs_WNT_down',
                            ifelse(peaks$rep.logFC.WNT_vs_BMP.WNT_LI > 0.5 & peaks$FDR < 0.05, 'WNT_vs_BMP.WNT_LI_up',
                                   ifelse(peaks$rep.logFC.WNT_vs_BMP.WNT_LI < -0.5 & peaks$FDR < 0.05, 'BMP.WNT_LI_vs_WNT_down',
                                          ifelse(peaks$rep.logFC.BMP.WNT_L_vs_LI > 0.5 & peaks$FDR < 0.05, 'BMP.WNT_L_vs_LI_up',
                                                 ifelse(peaks$rep.logFC.BMP.WNT_L_vs_LI < -0.5 & peaks$FDR < 0.05, 'BMP.WNT_LI_vs_L_down', 'Common'))))))
summary(peaks$name)
ggplot(data = data,
       mapping = aes(x = rep.logFC.WNT_vs_BMP.WNT_L,
                     y = -log10(PValue),
                     color = Significant,
                     shape = Annotated,
                     label = overlap)) +
  rasterize(geom_point()) +
  geom_text_repel(data = data |> filter(Annotated & Significant & abs(rep.logFC.WNT_vs_BMP.WNT_L) > 0.5)) +
  geom_vline(xintercept = 0, lty = 2) +
  theme_bw() +
  ggtitle("Differential Accessibility: WNT vs BMP.WNT_L")

ggplot(data = data,
       mapping = aes(x = rep.logFC.WNT_vs_BMP.WNT_LI,
                     y = -log10(PValue),
                     color = Significant,
                     shape = Annotated,
                     label = overlap)) +
  rasterize(geom_point()) +
  geom_text_repel(data = data |> filter(Annotated & Significant & abs(rep.logFC.WNT_vs_BMP.WNT_LI) > 0.5)) +
  geom_vline(xintercept = 0, lty = 2) +
  theme_bw() +
  ggtitle("Differential Accessibility: WNT vs BMP.WNT_LI")

ggplot(data = data,
       mapping = aes(x = rep.logFC.BMP.WNT_L_vs_LI,
                     y = -log10(PValue),
                     color = Significant,
                     shape = Annotated,
                     label = overlap)) +
  rasterize(geom_point()) +
  geom_text_repel(data = data |> filter(Annotated & Significant & abs(rep.logFC.BMP.WNT_L_vs_LI) > 0.5)) +
  geom_vline(xintercept = 0, lty = 2) +
  theme_bw() +
  ggtitle("Differential Accessibility: BMP.WNT_L vs BMP.WNT_LI")

# Export the results
write.table(as.data.frame(peaks), file = gzfile('/collignon/Tawdarous/atac_seq_lab/results/csaw_result/peaks_three_conditions.tsv.gz'),
            row.names = FALSE, sep = '\t')
rtracklayer::export(peaks, '/collignon/Tawdarous/atac_seq_lab/results/csaw_result/peaks_three_conditions.bed')

# Display significant genes for each contrast
display_significant_genes <- function(data, logFC_column, PValue_column, condition) {
  significant_genes <- data %>%
    filter(Significant & Annotated & abs(.data[[logFC_column]]) > 4) %>%
    arrange(desc(abs(.data[[logFC_column]])))
  
  cat(sprintf("Top upregulated genes in %s:\n", condition))
  significant_genes %>%
    filter(.data[[logFC_column]] > 0) %>%
    select(overlap, logFC = .data[[logFC_column]], PValue = .data[[PValue_column]]) %>%
    head(10) %>%
    print(n = 10)
  
  cat(sprintf("\nTop downregulated genes in %s:\n", condition))
  significant_genes %>%
    filter(.data[[logFC_column]] < 0) %>%
    select(overlap, logFC = .data[[logFC_column]], PValue = .data[[PValue_column]]) %>%
    head(10) %>%
    print(n = 10)
}

# Display significant genes for each pairwise comparison
display_significant_genes(data, "rep.logFC.WNT_vs_BMP.WNT_L", "PValue", "WNT vs BMP.WNT_L")
display_significant_genes(data, "rep.logFC.WNT_vs_BMP.WNT_LI", "PValue", "WNT vs BMP.WNT_LI")
display_significant_genes(data, "rep.logFC.BMP.WNT_L_vs_LI", "PValue", "BMP.WNT_L vs BMP.WNT_LI")

#__________________________t06_BMP.WNT VS t06_WNT _______________________
contr2 <- makeContrasts(t06_BMP.WNT - t06_WNT, levels = colnames(moma))
res2 <- glmQLFTest(fit, contrast = contr2)
table(decideTests(res2, p.value = 0.05, lfc = 1.5))
tt2 <- topTags(res2, n=Inf, sort.by="none")$table[,c('logFC','logCPM','F','PValue','FDR')]

tt2 |>
  as_tibble() |>
  mutate(Windows = rownames(dgel), .before=1) |>
  arrange(FDR) |>
  head()


rowData(filtered.processed_peak) <- cbind(rowData(filtered.processed_peak), tt2)
merged <- mergeWindows(rowRanges(filtered.processed_peak), tol=300L)
str(merged,1)
merged$region
summary(width(merged$region))



anno <- detailRanges(merged$region,
                     txdb=TxDb.Mmusculus.UCSC.mm39.knownGene,
                     orgdb=org.Mm.eg.db, promoter=c(1500, 500), dist=10*1e3L)

merged$region$overlap <- anno$overlap
merged$region$left <- anno$left
merged$region$right <- anno$right
merged$region
head(merged$region$overlap)

region.stats <- combineTests(merged$id, tt2) |>
  as_tibble() |>
  mutate(rep.logCPM = tt2[rep.test,'logCPM'],
         Contrast = colnames(contr))
region.stats |>
  head()
peaks <- merged$region 

mcols(peaks) <- c(mcols(peaks), region.stats)

peaks$name <- ifelse(peaks$rep.logFC > 0.5 & peaks$FDR < 0.05, ' t06_BMP.WNT',
                     ifelse(peaks$rep.logFC < (-0.5) & peaks$FDR < 0.05, 't06_WNT', 'Common'))

table(peaks$name, peaks$direction)

data <- mcols(peaks) |>
  as_tibble() |>
  mutate(logFC = rep.logFC,
         Significant = FDR < 0.05,
         Annotated = overlap != "")


ggplot(data = data,
       mapping = aes(x = logFC,
                     y = -log10(PValue),
                     color = Significant,
                     shape = Annotated,
                     label = overlap)) +
  rasterize( geom_point() ) +
  geom_text_repel(data = data |> filter(Annotated & Significant & abs(logFC) > 0.5)) +
  geom_vline(xintercept = 0, lty = 2) +
  theme_bw()


write.table(as.data.frame(peaks), file=gzfile('/collignon/Tawdarous/atac_seq_lab/results/csaw_result/peaks_t06_WNT_VS_t24_WNT-06h_L.tsv.gz'),
            row.names=FALSE, sep='\t')
rtracklayer::export(peaks, '/collignon/Tawdarous/atac_seq_lab/results/csaw_result/peaks_t06_BMP.WNT_VS_t06_WNT.bed')
saveRDS(peaks, file='/collignon/Tawdarous/atac_seq_lab/results/csaw_result/peaks_t06_BMP.WNT_VS_t06_WNT.rds')

significant_genes <- data %>%
  filter(Significant & Annotated & abs(logFC) > 0.5) %>%
  arrange(desc(abs(logFC)))
significant_genes %>%
  select(overlap, logFC, PValue) %>%
  head(20) %>%
  print(n = 20)

# Upregulated genes
cat("Top upregulated genes:\n")
significant_genes %>%
  filter(logFC > 0.5) %>%
  select(overlap, logFC, PValue) %>%
  head(10) %>%
  print(n = 10)

# Downregulated genes
cat("\nTop downregulated genes:\n")
significant_genes %>%
  filter(logFC < 0.5) %>%
  select(overlap, logFC, PValue) %>%
  head(10) %>%
  print(n = 10)


upregulated_genes <- significant_genes %>%
  filter(logFC > 0) %>%
  select(overlap, logFC, PValue)
write.table(upregulated_genes, file = "/collignon/Tawdarous/atac_seq_lab/results/csaw_result/upregulated_genes_t06_BMP.WNT_VS_t06_WNT_1_5.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

downregulated_genes <- significant_genes %>%
  filter(logFC < 0) %>%
  select(overlap, logFC, PValue)
write.table(downregulated_genes, file = "/collignon/Tawdarous/atac_seq_lab/results/csaw_result/downregulated_genes_t06_BMP.WNT_VS_t06_WNT_1_5.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
#________tf
#Extract gene names from the overlap column
extract_gene_names <- function(overlap_column) {
  sapply(strsplit(overlap_column, ":"),
         function(x) x[1])
}

# Apply the function to extract gene names
upregulated_genes <- significant_genes %>%
  filter(logFC > 0) %>%
  mutate(GeneName = extract_gene_names(overlap))

downregulated_genes <- significant_genes %>%
  filter(logFC < 0) %>%
  mutate(GeneName = extract_gene_names(overlap))

# Check the first few extracted gene names
head(upregulated_genes$GeneName)
head(downregulated_genes$GeneName)

# Perform enrichment analysis using ChEA 2016 database
up_enrichment <- enrichr(upregulated_genes$GeneName, databases = "ChEA_2016")
down_enrichment <- enrichr(downregulated_genes$GeneName, databases = "ChEA_2016")


# Plot top enriched transcription factors for upregulated genes
up_enrichment_df <- up_enrichment[[1]]
top_up_tfs <- head(up_enrichment_df[order(up_enrichment_df$Adjusted.P.value), ], 10)
ggplot(top_up_tfs, aes(x = reorder(Term, -Adjusted.P.value), y = -log10(Adjusted.P.value))) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Top Enriched Transcription Factors (Upregulated)", x = "Transcription Factor", y = "-log10(Adjusted P-value)")

# Plot top enriched transcription factors for downregulated genes
down_enrichment_df <- down_enrichment[[1]]
top_down_tfs <- head(down_enrichment_df[order(down_enrichment_df$Adjusted.P.value), ], 10)
ggplot(top_down_tfs, aes(x = reorder(Term, -Adjusted.P.value), y = -log10(Adjusted.P.value))) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Top Enriched Transcription Factors (Downregulated)", x = "Transcription Factor", y = "-log10(Adjusted P-value)")

#__________________________t06_BMP.WNT VS t24_BMP.WNT_06h_L _______________________
contr2 <- makeContrasts(t06_BMP.WNT - t24_BMP.WNT_06h_L, levels = colnames(moma))
res2 <- glmQLFTest(fit, contrast = contr2)
table(decideTests(res2, p.value = 0.05, lfc = 1.5))
tt2 <- topTags(res2, n=Inf, sort.by="none")$table[,c('logFC','logCPM','F','PValue','FDR')]

tt2 |>
  as_tibble() |>
  mutate(Windows = rownames(dgel), .before=1) |>
  arrange(FDR) |>
  head()


rowData(filtered.processed_peak) <- cbind(rowData(filtered.processed_peak), tt2)
merged <- mergeWindows(rowRanges(filtered.processed_peak), tol=300L)
str(merged,1)
merged$region
summary(width(merged$region))



anno <- detailRanges(merged$region,
                     txdb=TxDb.Mmusculus.UCSC.mm39.knownGene,
                     orgdb=org.Mm.eg.db, promoter=c(1500, 500), dist=10*1e3L)

merged$region$overlap <- anno$overlap
merged$region$left <- anno$left
merged$region$right <- anno$right
merged$region
head(merged$region$overlap)

region.stats <- combineTests(merged$id, tt2) |>
  as_tibble() |>
  mutate(rep.logCPM = tt2[rep.test,'logCPM'],
         Contrast = colnames(contr2))
region.stats |>
  head()
peaks <- merged$region 

mcols(peaks) <- c(mcols(peaks), region.stats)

peaks$name <- ifelse(peaks$rep.logFC > 2.5 & peaks$FDR < 0.05, ' t06_BMP.WNT',
                     ifelse(peaks$rep.logFC < (-2.5) & peaks$FDR < 0.05, 't24_BMP.WNT_06h_L', 'Common'))

table(peaks$name, peaks$direction)

data <- mcols(peaks) |>
  as_tibble() |>
  mutate(logFC = rep.logFC,
         Significant = FDR < 0.05,
         Annotated = overlap != "")


ggplot(data = data,
       mapping = aes(x = logFC,
                     y = -log10(PValue),
                     color = Significant,
                     shape = Annotated,
                     label = overlap)) +
  rasterize( geom_point() ) +
  geom_text_repel(data = data |> filter(Annotated & Significant & abs(logFC) > 0.5)) +
  geom_vline(xintercept = 0, lty = 2) +
  theme_bw()


write.table(as.data.frame(peaks), file=gzfile('/collignon/Tawdarous/atac_seq_lab/results/csaw_result/peaks_t06_BMP.WNT_VS_t24_BMP.WNT_06h_L_logFC_05.tsv.gz'),
            row.names=FALSE, sep='\t')
rtracklayer::export(peaks, '/collignon/Tawdarous/atac_seq_lab/results/csaw_result/peaks_t06_BMP.WNT_VS_t24_BMP.WNT_06h_L_logFC_05.bed')
saveRDS(peaks, file='/collignon/Tawdarous/atac_seq_lab/results/csaw_result/peaks_t06_BMP.WNT_VS_t24_BMP.WNT_06h_L_logFC_05.rds')



#-------------t24_BMP.WNT_06h_L vs t24_WNT_06h_L --------
  
  
contr2 <- makeContrasts(t24_BMP.WNT_06h_L - t24_WNT_06h_L, levels = colnames(moma))
res2 <- glmQLFTest(fit, contrast = contr2)
table(decideTests(res2, p.value = 0.01, lfc = 1))
tt2 <- topTags(res2, n=Inf, sort.by="none")$table[,c('logFC','logCPM','F','PValue','FDR')]

tt2 |>
  as_tibble() |>
  mutate(Windows = rownames(dgel), .before=1) |>
  arrange(FDR) |>
  head()


rowData(filtered.processed_peak) <- cbind(rowData(filtered.processed_peak), tt2)
merged <- mergeWindows(rowRanges(filtered.processed_peak), tol=300L)
str(merged,1)
merged$region
summary(width(merged$region))



anno <- detailRanges(merged$region,
                     txdb=TxDb.Mmusculus.UCSC.mm39.knownGene,
                     orgdb=org.Mm.eg.db, promoter=c(1500, 500), dist=10*1e3L)

merged$region$overlap <- anno$overlap
merged$region$left <- anno$left
merged$region$right <- anno$right
merged$region
head(merged$region$overlap)

region.stats <- combineTests(merged$id, tt2) |>
  as_tibble() |>
  mutate(rep.logCPM = tt2[rep.test,'logCPM'],
         Contrast = colnames(contr2))
region.stats |>
  head()
peaks <- merged$region 

mcols(peaks) <- c(mcols(peaks), region.stats)

peaks$name <- ifelse(peaks$rep.logFC > 1 & peaks$FDR < 0.01 & peaks$PValue<0.01, ' t24_BMP.WNT_06h_L',
                     ifelse(peaks$rep.logFC < (-1) & peaks$FDR < 0.01 & peaks$PValue<0.01, 't24_WNT_06h_L', 'Common'))

table(peaks$name, peaks$direction)


data <- mcols(peaks) |>
  as_tibble() |>
  mutate(logFC = rep.logFC,
         Significant = abs(logFC) > 1 & PValue < 0.01,
         Annotated = overlap != "")

# Generate the plot
ggplot(data = data,
       mapping = aes(x = logFC,
                     y = -log10(PValue),
                     color = Significant,
                     shape = Annotated,
                     label = overlap)) +
  rasterize(geom_point()) +
  geom_text_repel(data = data |> filter(Annotated & Significant)) +
  geom_vline(xintercept = 1, lty = 2) +
  geom_vline(xintercept = -1, lty = 2) +
  geom_hline(yintercept = -log10(0.01), lty = 2) +
  theme_bw()





write.table(as.data.frame(peaks), file=gzfile('/collignon/Tawdarous/atac_seq_lab/results/csaw_result/peaks_t24_BMP.WNT_06h_L_vs_t24_WNT_06h_L.tsv.gz'),
            row.names=FALSE, sep='\t')
rtracklayer::export(peaks, '/collignon/Tawdarous/atac_seq_lab/results/csaw_result/peaks_t24_BMP.WNT_06h_L_vs_t24_WNT_06h_L.bed')
saveRDS(peaks, file='/collignon/Tawdarous/atac_seq_lab/results/csaw_result/peaks_t24_BMP.WNT_06h_L_vs_t24_WNT_06h_L.rds')

significant_genes_t24_BMP.WNT_06h_L_vs_t24_WNT_06h_L <- data %>%
  filter(Significant & Annotated ) %>%
  arrange(desc(abs(logFC)))
significant_genes_t24_BMP.WNT_06h_L_vs_t24_WNT_06h_L %>%
  select(overlap, logFC, PValue) %>%
  head(20) %>%
  print(n = 20)

upregulated_genes_t24_BMP.WNT_06h_L_vs_t24_WNT_06h_L <- significant_genes_t24_BMP.WNT_06h_L_vs_t24_WNT_06h_L %>%
  filter(logFC > 1) %>%
  select(overlap, logFC, PValue)
write.table(upregulated_genes, file = "/collignon/Tawdarous/atac_seq_lab/results/csaw_result/upregulated_genes_t24_BMP.WNT_06h_L_vs_t24_WNT_06h_L.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

downregulated_genes <- significant_genes_t24_BMP.WNT_06h_L_vs_t24_WNT_06h_L %>%
  filter(logFC < -1) %>%
  select(overlap, logFC, PValue)
write.table(downregulated_genes, file = "/collignon/Tawdarous/atac_seq_lab/results/csaw_result/downregulated_genes_t24_BMP.WNT_06h_L_vs_t24_WNT_06h_L.tsv", sep = "\t", row.names = FALSE, quote = FALSE)



extract_gene_names <- function(overlap_column) {
  sapply(strsplit(overlap_column, ":"),
         function(x) x[1])
}

# Apply the function to extract gene names
upregulated_genes <- significant_genes_t24_BMP.WNT_06h_L_vs_t24_WNT_06h_L %>%
  filter(logFC > 1) %>%
  mutate(GeneName = extract_gene_names(overlap))

downregulated_genes <- significant_genes_t24_BMP.WNT_06h_L_vs_t24_WNT_06h_L %>%
  filter(logFC < -1) %>%
  mutate(GeneName = extract_gene_names(overlap))

# Check the first few extracted gene names
head(upregulated_genes$GeneName)
head(downregulated_genes$GeneName)


gene_list_upregulated_genes <- bitr(upregulated_genes$GeneName, fromType = "SYMBOL", 
                                    toType = "ENTREZID", OrgDb = org.Mm.eg.db)
gene_list_downregulated_genes <- bitr(downregulated_genes$GeneName, fromType = "SYMBOL", 
                                    toType = "ENTREZID", OrgDb = org.Mm.eg.db)
ego_BP_increased <- enrichGO(gene = gene_list_upregulated_genes$ENTREZID,
                             OrgDb = org.Mm.eg.db,
                             ont = "BP", #or BP
                             pAdjustMethod = "BH",
                             pvalueCutoff = 0.05,
                             qvalueCutoff = 0.05)
dotplot(ego_BP_increased, title = "ontologie: +++ dans t24_BMP.WNT_06h_L")

ego_BP_decreased <- enrichGO(gene = gene_list_downregulated_genes$ENTREZID,
                             OrgDb = org.Mm.eg.db,
                             ont = "BP",#or BP or CC
                             pAdjustMethod = "BH",
                             pvalueCutoff = 0.05,
                             qvalueCutoff = 0.05)
dotplot(ego_BP_decreased, title = "ontologie: +++ dans t24_WNT_06h_L")
#____________________KEGG_____________
# Upregulated genes
kegg_up <- enrichKEGG(gene = gene_list_upregulated_genes$ENTREZID,
                      organism = 'mmu', # mmu for mouse
                      pvalueCutoff = 0.05)

# Downregulated genes
kegg_down <- enrichKEGG(gene = gene_list_downregulated_genes$ENTREZID,
                        organism = 'mmu', # mmu for mouse
                        pvalueCutoff = 0.01)

# Plot results
dotplot(kegg_up, title = "KEGG Pathways: Upregulated Genes, t24_BMP.WNT_06h_L")
dotplot(kegg_down, title = "KEGG Pathways: Downregulated Genes, t24_WNT_06h_L")



# Upregulated genes
reactome_up <- enrichPathway(gene = gene_list_upregulated_genes$ENTREZID,
                             organism = 'mouse', # specify organism
                             pvalueCutoff = 0.05)

# Downregulated genes
reactome_down <- enrichPathway(gene = gene_list_downregulated_genes$ENTREZID,
                               organism = 'mouse', # specify organism
                               pvalueCutoff = 0.05)

# Plot results
dotplot(reactome_up, title = "Reactome Pathways: Upregulated Genes")
dotplot(reactome_down, title = "Reactome Pathways: Downregulated Genes")


#-------------t24_BMP.WNT_06h_L vs t24_BMP.WNT_06h_LI --------



contr2 <- makeContrasts(t24_BMP.WNT_06h_L - t24_BMP.WNT_06h_LI, levels = colnames(moma))
res2 <- glmQLFTest(fit, contrast = contr2)
table(decideTests(res2, p.value = 0.01, lfc = 1))
tt2 <- topTags(res2, n=Inf, sort.by="none")$table[,c('logFC','logCPM','F','PValue','FDR')]

tt2 |>
  as_tibble() |>
  mutate(Windows = rownames(dgel), .before=1) |>
  arrange(FDR) |>
  head()


rowData(filtered.processed_peak) <- cbind(rowData(filtered.processed_peak), tt2)
merged <- mergeWindows(rowRanges(filtered.processed_peak), tol=300L)
str(merged,1)
merged$region
summary(width(merged$region))



anno <- detailRanges(merged$region,
                     txdb=TxDb.Mmusculus.UCSC.mm39.knownGene,
                     orgdb=org.Mm.eg.db, promoter=c(1500, 500), dist=10*1e3L)

merged$region$overlap <- anno$overlap
merged$region$left <- anno$left
merged$region$right <- anno$right
merged$region
head(merged$region$overlap)

region.stats <- combineTests(merged$id, tt2) |>
  as_tibble() |>
  mutate(rep.logCPM = tt2[rep.test,'logCPM'],
         Contrast = colnames(contr2))
region.stats |>
  head()
peaks <- merged$region 

mcols(peaks) <- c(mcols(peaks), region.stats)

peaks$name <- ifelse(peaks$rep.logFC > 1 & peaks$FDR < 0.01 & peaks$PValue<0.01, ' t24_BMP.t24_BMP.WNT_06h_L',
                     ifelse(peaks$rep.logFC < (-1) & peaks$FDR < 0.01 & peaks$PValue<0.01, 't24_BMP.WNT_06h_LI', 'Common'))

table(peaks$name, peaks$direction)


data <- mcols(peaks) |>
  as_tibble() |>
  mutate(logFC = rep.logFC,
         Significant = abs(logFC) > 1 & PValue < 0.01,
         Annotated = overlap != "")

# Generate the plot
ggplot(data = data,
       mapping = aes(x = logFC,
                     y = -log10(PValue),
                     color = Significant,
                     shape = Annotated,
                     label = overlap)) +
  rasterize(geom_point()) +
  geom_text_repel(data = data |> filter(Annotated & Significant), , max.overlaps = 10) +
  geom_vline(xintercept = 1, lty = 2) +
  geom_vline(xintercept = -1, lty = 2) +
  geom_hline(yintercept = -log10(0.01), lty = 2) +
  theme_bw()





write.table(as.data.frame(peaks), file=gzfile('/collignon/Tawdarous/atac_seq_lab/results/csaw_result/peaks_t24_BMP.WNT_06h_L_vs_t24_BMP.WNT_06h_LI.tsv.gz'),
            row.names=FALSE, sep='\t')
rtracklayer::export(peaks, '/collignon/Tawdarous/atac_seq_lab/results/csaw_result/peaks_t24_BMP.WNT_06h_L_vs_t24_BMP.WNT_06h_LI.bed')
saveRDS(peaks, file='/collignon/Tawdarous/atac_seq_lab/results/csaw_result/peaks_t24_BMP.WNT_06h_L_vs_t24_BMP.WNT_06h_LI.rds')

significant_genes_t24_BMP.WNT_06h_L_vs_t24_BMP.WNT_06h_LI<- data %>%
  filter(Significant & Annotated ) %>%
  arrange(desc(abs(logFC)))
significant_genes_t24_BMP.WNT_06h_L_vs_t24_BMP.WNT_06h_LI %>%
  select(overlap, logFC, PValue) %>%
  head(20) %>%
  print(n = 20)

upregulated_genes_t24_BMP.WNT_06h_L_vs_t24_BMP.WNT_06h_LI  <- significant_genes_t24_BMP.WNT_06h_L_vs_t24_BMP.WNT_06h_LI %>%
  filter(logFC > 1) %>%
  select(overlap, logFC, PValue)
write.table(upregulated_genes_t24_BMP.WNT_06h_L_vs_t24_BMP.WNT_06h_LI, file = "/collignon/Tawdarous/atac_seq_lab/results/csaw_result/upregulated_genes_t24_BMP.WNT_06h_L_vs_t24_BMP.WNT_06h_LI.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

downregulated_genes_t24_BMP.WNT_06h_L_vs_t24_BMP.WNT_06h_LI  <- significant_genes_t24_BMP.WNT_06h_L_vs_t24_BMP.WNT_06h_LI %>%
  filter(logFC < -1) %>%
  select(overlap, logFC, PValue)
write.table(downregulated_genes_t24_BMP.WNT_06h_L_vs_t24_BMP.WNT_06h_LI, file = "/collignon/Tawdarous/atac_seq_lab/results/csaw_result/downregulated_genes_t24_BMP.WNT_06h_L_vs_t24_BMP.WNT_06h_LI.tsv", sep = "\t", row.names = FALSE, quote = FALSE)



extract_gene_names <- function(overlap_column) {
  sapply(strsplit(overlap_column, ":"),
         function(x) x[1])
}

# Apply the function to extract gene names
upregulated_genes_t24_BMP.WNT_06h_L_vs_t24_BMP.WNT_06h_LI <- significant_genes_t24_BMP.WNT_06h_L_vs_t24_BMP.WNT_06h_LI %>%
  filter(logFC > 1) %>%
  mutate(GeneName = extract_gene_names(overlap))

down_genes_t24_BMP.WNT_06h_L_vs_t24_BMP.WNT_06h_LI <- significant_genes_t24_BMP.WNT_06h_L_vs_t24_BMP.WNT_06h_LI %>%
  filter(logFC < -1) %>%
  mutate(GeneName = extract_gene_names(overlap))

# Check the first few extracted gene names
head(upregulated_genes_t24_BMP.WNT_06h_L_vs_t24_BMP.WNT_06h_LI$GeneName)
head(down_genes_t24_BMP.WNT_06h_L_vs_t24_BMP.WNT_06h_LI$GeneName)


gene_list_upregulated_genes <- bitr(upregulated_genes_t24_BMP.WNT_06h_L_vs_t24_BMP.WNT_06h_LI$GeneName, fromType = "SYMBOL", 
                                    toType = "ENTREZID", OrgDb = org.Mm.eg.db)
gene_list_downregulated_genes <- bitr(down_genes_t24_BMP.WNT_06h_L_vs_t24_BMP.WNT_06h_LI$GeneName, fromType = "SYMBOL", 
                                      toType = "ENTREZID", OrgDb = org.Mm.eg.db)
ego_BP_increased <- enrichGO(gene = gene_list_upregulated_genes$ENTREZID,
                             OrgDb = org.Mm.eg.db,
                             ont = "CC", #or BP
                             pAdjustMethod = "BH",
                             pvalueCutoff = 0.05,
                             qvalueCutoff = 0.05)
dotplot(ego_BP_increased, title = "upregulated_genes_t24_BMP.WNT_06h_L")

ego_BP_decreased <- enrichGO(gene = gene_list_downregulated_genes$ENTREZID,
                             OrgDb = org.Mm.eg.db,
                             ont = "CC",#or BP or CC
                             pAdjustMethod = "BH",
                             pvalueCutoff = 0.05,
                             qvalueCutoff = 0.05)
dotplot(ego_BP_decreased, title = "upregulated_genes_t24_BMP.WNT_06h_LI")
#____________________KEGG_____________
# Upregulated genes
kegg_up <- enrichKEGG(gene = gene_list_upregulated_genes$ENTREZID,
                      organism = 'mmu', # mmu for mouse
                      pvalueCutoff = 0.01)

# Downregulated genes
kegg_down <- enrichKEGG(gene = gene_list_downregulated_genes$ENTREZID,
                        organism = 'mmu', # mmu for mouse
                        pvalueCutoff = 0.01)

# Plot results
dotplot(kegg_up, title = "Upregulated Genes, t24_BMP.WNT_06h_L")
dotplot(kegg_down, title = "upregulated_genes_t24_BMP.WNT_06h_LI")



# Upregulated genes
reactome_up <- enrichPathway(gene = gene_list_upregulated_genes$ENTREZID,
                             organism = 'mouse', # specify organism
                             pvalueCutoff = 0.01)

# Downregulated genes
reactome_down <- enrichPathway(gene = gene_list_downregulated_genes$ENTREZID,
                               organism = 'mouse', # specify organism
                               pvalueCutoff = 0.01)

# Plot results
dotplot(reactome_up, title = "Reactome Pathways: Upregulated Genes")
dotplot(reactome_down, title = "Reactome Pathways: Downregulated Genes")

#----------------------t24_BMP+WNT-06h_LI_VS_t24_WNT-06h_L------------------------



contr2 <- makeContrasts(t24_BMP.WNT_06h_LI - t24_WNT_06h_L, levels = colnames(moma))
res2 <- glmQLFTest(fit, contrast = contr2)
table(decideTests(res2, p.value = 0.01, lfc = 1))
tt2 <- topTags(res2, n=Inf, sort.by="none")$table[,c('logFC','logCPM','F','PValue','FDR')]

tt2 |>
  as_tibble() |>
  mutate(Windows = rownames(dgel), .before=1) |>
  arrange(FDR) |>
  head()


rowData(filtered.processed_peak) <- cbind(rowData(filtered.processed_peak), tt2)
merged <- mergeWindows(rowRanges(filtered.processed_peak), tol=300L)
str(merged,1)
merged$region
summary(width(merged$region))



anno <- detailRanges(merged$region,
                     txdb=TxDb.Mmusculus.UCSC.mm39.knownGene,
                     orgdb=org.Mm.eg.db, promoter=c(3000, 3000), dist=10*1e3L)

merged$region$overlap <- anno$overlap
merged$region$left <- anno$left
merged$region$right <- anno$right
merged$region
head(merged$region$overlap)

region.stats <- combineTests(merged$id, tt2) |>
  as_tibble() |>
  mutate(rep.logCPM = tt2[rep.test,'logCPM'],
         Contrast = colnames(contr2))
region.stats |>
  head()
peaks <- merged$region 

mcols(peaks) <- c(mcols(peaks), region.stats)

peaks$name <- ifelse(peaks$rep.logFC > 1 & peaks$FDR < 0.01 & peaks$PValue<0.01, ' t24_BMP.WNT_06h_LI',
                     ifelse(peaks$rep.logFC < (-1) & peaks$FDR < 0.01 & peaks$PValue<0.01, 't24_WNT_06h_L', 'Common'))

table(peaks$name, peaks$direction)


data <- mcols(peaks) |>
  as_tibble() |>
  mutate(logFC = rep.logFC,
         Significant = abs(logFC) > 1 & PValue < 0.01,
         Annotated = overlap != "")

# Generate the plot
ggplot(data = data,
       mapping = aes(x = logFC,
                     y = -log10(PValue),
                     color = Significant,
                     shape = Annotated,
                     label = overlap)) +
  rasterize(geom_point()) +
  geom_text_repel(data = data |> filter(Annotated & Significant), , max.overlaps = 10) +
  geom_vline(xintercept = 1, lty = 2) +
  geom_vline(xintercept = -1, lty = 2) +
  geom_hline(yintercept = -log10(0.01), lty = 2) +
  theme_bw()





write.table(as.data.frame(peaks), file=gzfile('/collignon/Tawdarous/atac_seq_lab/results/csaw_result/peaks_t24_BMP+WNT-06h_LI_VS_t24_WNT-06h_Ltsv.gz'),
            row.names=FALSE, sep='\t')
rtracklayer::export(peaks, '/collignon/Tawdarous/atac_seq_lab/results/csaw_result/peaks_t24_BMP+WNT-06h_LI_VS_t24_WNT-06h_L.bed')
saveRDS(peaks, file='/collignon/Tawdarous/atac_seq_lab/results/csaw_result/peaks_t24_BMP+WNT-06h_LI_VS_t24_WNT-06h_L.rds')

significant_genes_t24_BMP_WNT_06h_LI_VS_t24_WNT_06h_L<- data %>%
  filter(Significant & Annotated ) %>%
  arrange(desc(abs(logFC)))
significant_genes_t24_BMP_WNT_06h_LI_VS_t24_WNT_06h_L %>%
  select(overlap, logFC, PValue) %>%
  head(20) %>%
  print(n = 20)

upregulated_genes_t24_BMP_WNT_06h_LI_VS_t24_WNT_06h_L  <- significant_genes_t24_BMP_WNT_06h_LI_VS_t24_WNT_06h_L %>%
  filter(logFC > 1) %>%
  select(overlap, logFC, PValue)
write.table(upregulated_genes_t24_BMP_WNT_06h_LI_VS_t24_WNT_06h_L, file = "/collignon/Tawdarous/atac_seq_lab/results/csaw_result/upregulated_genes_t24_BMP+WNT-06h_LI_VS_t24_WNT-06h_L.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

downregulated_genes_t24_BMP_WNT_06h_LI_VS_t24_WNT_06h_L  <- significant_genes_t24_BMP_WNT_06h_LI_VS_t24_WNT_06h_L %>%
  filter(logFC < -1) %>%
  select(overlap, logFC, PValue)
write.table(downregulated_genes_t24_BMP_WNT_06h_LI_VS_t24_WNT_06h_L, file = "/collignon/Tawdarous/atac_seq_lab/results/csaw_result/downregulated_genes_t24_BMP+WNT-06h_LI_VS_t24_WNT-06h_L.tsv", sep = "\t", row.names = FALSE, quote = FALSE)




extract_gene_names <- function(overlap_column) {
  sapply(strsplit(overlap_column, ":"),
         function(x) x[1])
}

# Apply the function to extract gene names
upregulated_genes_t24_BMP_WNT_06h_LI_VS_t24_WNT_06h_L <- significant_genes_t24_BMP_WNT_06h_LI_VS_t24_WNT_06h_L %>%
  filter(logFC > 1) %>%
  mutate(GeneName = extract_gene_names(overlap))

downregulated_genes_t24_BMP_WNT_06h_LI_VS_t24_WNT_06h_L <- significant_genes_t24_BMP_WNT_06h_LI_VS_t24_WNT_06h_L %>%
  filter(logFC < -1) %>%
  mutate(GeneName = extract_gene_names(overlap))

# Check the first few extracted gene names
head(upregulated_genes_t24_BMP_WNT_06h_LI_VS_t24_WNT_06h_L$GeneName)
head(downregulated_genes_t24_BMP_WNT_06h_LI_VS_t24_WNT_06h_L$GeneName)



gene_list_upregulated_genes <- bitr(upregulated_genes_t24_BMP_WNT_06h_LI_VS_t24_WNT_06h_L$GeneName, fromType = "SYMBOL", 
                                    toType = "ENTREZID", OrgDb = org.Mm.eg.db)
gene_list_downregulated_genes <- bitr(downregulated_genes_t24_BMP_WNT_06h_LI_VS_t24_WNT_06h_L$GeneName, fromType = "SYMBOL", 
                                      toType = "ENTREZID", OrgDb = org.Mm.eg.db)
ego_BP_increased <- enrichGO(gene = gene_list_upregulated_genes$ENTREZID,
                             OrgDb = org.Mm.eg.db,
                             ont = "CC", #or BP or CC
                             pAdjustMethod = "BH",
                             pvalueCutoff = 0.05,
                             qvalueCutoff = 0.05)
dotplot(ego_BP_increased, title = "upregulated_genes_t24_BMP.WNT_06h_LI")

ego_BP_decreased <- enrichGO(gene = gene_list_downregulated_genes$ENTREZID,
                             OrgDb = org.Mm.eg.db,
                             ont = "CC",#or BP or CC
                             pAdjustMethod = "BH",
                             pvalueCutoff = 0.05,
                             qvalueCutoff = 0.05)
dotplot(ego_BP_decreased, title = "upregulated_genes_t24_WNT_06h_L")
#____________________KEGG_____________
# Upregulated genes
kegg_up <- enrichKEGG(gene = gene_list_upregulated_genes$ENTREZID,
                      organism = 'mmu', # mmu for mouse
                      pvalueCutoff = 0.01)

# Downregulated genes
kegg_down <- enrichKEGG(gene = gene_list_downregulated_genes$ENTREZID,
                        organism = 'mmu', # mmu for mouse
                        pvalueCutoff = 0.01)

# Plot results
dotplot(kegg_up, title = "upregulated_genes_t24_BMP.WNT_06h_LI")
dotplot(kegg_down, title = "upregulated_genes_t24_WNT_06h_L")



# Upregulated genes
reactome_up <- enrichPathway(gene = gene_list_upregulated_genes$ENTREZID,
                             organism = 'mouse', # specify organism
                             pvalueCutoff = 0.01)

# Downregulated genes
reactome_down <- enrichPathway(gene = gene_list_downregulated_genes$ENTREZID,
                               organism = 'mouse', # specify organism
                               pvalueCutoff = 0.01)

# Plot results
dotplot(reactome_up, title = "Reactome Pathways: Upregulated Genes")
dotplot(reactome_down, title = "Reactome Pathways: Downregulated Genes")

#__________________________t24_WNT_06h_L VS t24_BMP.WNT_06h_LX1 _______________________
contr <- makeContrasts(t24_WNT_06h_L - t24_BMP.WNT_06h_LX1, levels = colnames(moma))
res <- glmQLFTest(fit, contrast = contr)
table(decideTests(res, p.value = 0.05, lfc = 1.5))
tt <- topTags(res, n=Inf, sort.by="none")$table[,c('logFC','logCPM','F','PValue','FDR')]

tt |>
  as_tibble() |>
  mutate(Windows = rownames(dgel), .before=1) |>
  arrange(FDR) |>
  head()


rowData(filtered.processed_peak) <- cbind(rowData(filtered.processed_peak), tt)
merged <- mergeWindows(rowRanges(filtered.processed_peak), tol=300L)
str(merged,1)
merged$region
summary(width(merged$region))



anno <- detailRanges(merged$region,
                     txdb=TxDb.Mmusculus.UCSC.mm39.knownGene,
                     orgdb=org.Mm.eg.db, promoter=c(1500, 500), dist=10*1e3L)

merged$region$overlap <- anno$overlap
merged$region$left <- anno$left
merged$region$right <- anno$right
merged$region
head(merged$region$overlap)

region.stats <- combineTests(merged$id, tt) |>
  as_tibble() |>
  mutate(rep.logCPM = tt[rep.test,'logCPM'],
         Contrast = colnames(contr))
region.stats |>
  head()
peaks <- merged$region 

mcols(peaks) <- c(mcols(peaks), region.stats)

peaks$name <- ifelse(peaks$rep.logFC > 1.5 & peaks$FDR < 0.05, 't24_WNT_06h_L',
                     ifelse(peaks$rep.logFC < (-1.5) & peaks$FDR < 0.05, 't24_BMP.WNT_06h_LX1', 'Common'))

table(peaks$name, peaks$direction)

data <- mcols(peaks) |>
  as_tibble() |>
  mutate(logFC = rep.logFC,
         Significant = FDR < 0.05,
         Annotated = overlap != "")


ggplot(data = data,
       mapping = aes(x = logFC,
                     y = -log10(PValue),
                     color = Significant,
                     shape = Annotated,
                     label = overlap)) +
  rasterize( geom_point() ) +
  geom_text_repel(data = data |> filter(Annotated & Significant & abs(logFC) > 1.5)) +
  geom_vline(xintercept = 0, lty = 2) +
  theme_bw()


write.table(as.data.frame(peaks), file=gzfile('/collignon/Tawdarous/atac_seq_lab/results/csaw_result/peaks_t24_WNT_06h_L_VS_t24_BMP.WNT_06h_LX1.tsv.gz'),
            row.names=FALSE, sep='\t')
rtracklayer::export(peaks, '/collignon/Tawdarous/atac_seq_lab/results/csaw_result/peaks_t24_WNT_06h_L_VS_t24_BMP.WNT_06h_LX1.bed')

significant_genes <- data %>%
  filter(Significant & Annotated & abs(logFC) > 2.5) %>%
  arrange(desc(abs(logFC)))
significant_genes %>%
  select(overlap, logFC, PValue) %>%
  head(20) %>%
  print(n = 20)

# Upregulated genes
cat("Top upregulated genes:\n")
significant_genes %>%
  filter(logFC > 0) %>%
  select(overlap, logFC, PValue) %>%
  head(10) %>%
  print(n = 10)

# Downregulated genes
cat("\nTop downregulated genes:\n")
significant_genes %>%
  filter(logFC < 0) %>%
  select(overlap, logFC, PValue) %>%
  head(10) %>%
  print(n = 10)


#-----------------------------------CHIPSEEKER-------------------------

files<- c(t06_WNT.rep1,t06_WNT.rep2,t06_WNT.rep3,
              t24_WNT_06h_L_rep1,t24_WNT_06h_L_rep2,t24_WNT_06h_L_rep3)
txdb <- TxDb.Mmusculus.UCSC.mm39.knownGene
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrixList <- lapply(files, getTagMatrix, windows=promoter)

#-----------------------------------Peak calling 2-------------------------

filterStat <- filterWindowsGlobal(processed_peak, processed_backgroud)
str(filterStat)
globalBG <- median( filterStat$back.abundances )
globalBG
head(filterStat$abundances - globalBG)
keep <- filterStat$filter > log2(3)
sum(keep)
mean(keep)
min(filterStat$abundances[keep])
globalBG + log2(3)



ggplot(data = tibble(x=filterStat$abundances),
       mapping = aes(x)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = globalBG, color='red') +
  geom_vline(xintercept = globalBG + log2(3), color='blue') +
  labs(x="Abundance", y="Frequency", title="Histogram of Window Abundances") +
  theme_bw(base_size=15)


rowData(processed_peak)$enrichment <- filterStat$filter
filtered.processed_peak <- processed_peak[keep,]
rtracklayer::export(rowRanges(filtered.processed_peak),
                    '/collignon/Tawdarous/atac_seq_lab/results/csaw_result/filtered_global.bed')
###############
filtered.processed_peak <- normOffsets(filtered.processed_peak)
logCPM <- csaw::calculateCPM(filtered.processed_peak, log=TRUE, use.offsets = TRUE)
colnames(filtered.processed_peak) <- sample_names
colData(filtered.processed_peak)$CellType <- cell_types

colnames(logCPM) <- colnames(filtered.processed_peak)

ct <- split(colnames(logCPM), colData(filtered.processed_peak)$CellType)
combi <- list(
  c('t06_BMP+WNT', 't06_WNT'),
  c('t24_BMP+WNT-06h_LI', 't24_BMP+WNT-06h_L'),
  c('t24_BMP+WNT-06h_L', 't24_WNT-06h_L')
)


df <- vapply(ct, function(x) {
  rowMeans(logCPM[,x,drop=FALSE])
}, rep(0, nrow(filtered.processed_peak))) |>
  as_tibble() |>
  mutate(overlapTSS = rowData(filtered.processed_peak)$overlapTSS,
         GC = rowData(filtered.processed_peak)[,'G|C'])

patchwork::wrap_plots(lapply(combi, function(x) plotMA(df, x[1], x[2])), ncol=2)


# Convert to tibble and add overlapTSS and GC columns
df <- as_tibble(df) %>%
  mutate(overlapTSS = rowData(filtered.processed_peak)$overlapTSS,
         GC = rowData(filtered.processed_peak)[, 'G|C'])

# Verify the columns in the dataframe
print(colnames(df))

# Define combinations for comparisons
combi <- list(
  c('t06_BMP+WNT', 't06_WNT'),
  c('t24_BMP+WNT-06h_LI', 't24_BMP+WNT-06h_L'),
  c('t24_BMP+WNT-06h_L', 't24_WNT-06h_L')
)

# Function to plot MA plots
plotMA <- function(df, s1, s2) {
  ggplot(df, mapping = aes(x = .data[[s1]] + .data[[s2]], y = .data[[s1]] - .data[[s2]])) +
    geom_hex(bins = 100, aes(fill = after_stat(density)^(1/16)), show.legend = FALSE) +
    rasterize(geom_point(data = subset(df, overlapTSS), col = 'red', pch = '.')) +
    rasterize(geom_smooth(data = subset(df, overlapTSS), col = 'plum', se = FALSE)) +
    geom_hline(yintercept = 0, lty = 2) +
    labs(x = "A", y = "M", title = sprintf("%s vs. %s", s1, s2)) +
    theme_bw(base_size = 15)
}

wrap_plots(lapply(combi, function(x) plotMA(df, x[1], x[2])), ncol=2)

########
logCPM <- csaw::calculateCPM(filtered.processed_peak, log=TRUE, use.offsets = TRUE)

# Ensure the correct sample group assignment
colData(filtered.processed_peak)$SampleGroup <- gsub("_rep\\d+.*", "", colnames(filtered.processed_peak))

# Split samples by SampleGroup
ct <- split(colnames(logCPM), colData(filtered.processed_peak)$SampleGroup)

# Create a dataframe for MA plot
df <- vapply(ct, function(x) {
  if (length(x) > 0) {
    rowMeans(logCPM[, x, drop = FALSE])
  } else {
    rep(NA, nrow(logCPM))
  }
}, numeric(nrow(logCPM)))

# Convert to tibble and add overlapTSS and GC columns
df <- as_tibble(df) %>%
  mutate(overlapTSS = rowData(filtered.processed_peak)$overlapTSS,
         GC = rowData(filtered.processed_peak)[, 'G|C'])

# Verify the columns in the dataframe
print(colnames(df))

# Define combinations for comparisons
combi <- list(
  c('t06_BMP_WNT', 't06_WNT'),
  c('t24_BMP_WNT_06h_LI', 't24_BMP_WNT_06h_L'),
  c('t24_BMP_WNT_06h_L', 't24_WNT_06h_L')
)

# Function to plot MA plots
plotMA <- function(df, s1, s2) {
  ggplot(df, mapping = aes(x = .data[[s1]] + .data[[s2]], y = .data[[s1]] - .data[[s2]])) +
    geom_hex(bins = 100, aes(fill = after_stat(density)^(1/16)), show.legend = FALSE) +
    rasterize(geom_point(data = subset(df, overlapTSS), col = 'red', pch = '.')) +
    rasterize(geom_smooth(data = subset(df, overlapTSS), col = 'plum', se = FALSE)) +
    geom_hline(yintercept = 0, lty = 2) +
    labs(x = "A", y = "M", title = sprintf("%s vs. %s", s1, s2)) +
    theme_bw(base_size = 15)
}

# Plot MA plots for defined combinations
patchwork::wrap_plots(lapply(combi, function(x) plotMA(df, x[1], x[2])), ncol=2)

##########

#-----------------------------------Add sample annotation-------------------------
#coldata
sample_names <- basename(bam.files)
cell_types <- sub("_rep\\d+.*", "", sample_names)

colData <- DataFrame(
  bam.files = bam.files,
  sample_names = sample_names,
  cell_types = cell_types
)
colData(processed_peak) <- colData
#head(colData(processed_peak))

#Row data
seqs <- getSeq(BSgenome.Mmusculus.UCSC.mm39, rowRanges(processed_peak))
rowData(processed_peak) <- letterFrequency(seqs, letters=c('GC', 'N'), as.prob=TRUE) # G or C, N frequency
#head(rowData(processed_peak))

is.agap <- rowData(processed_peak)[,'N'] > 2/3 # i.e. more than 100 of 150 bases
table(is.agap)/nrow(processed_peak)

is.smaller <- width(rowRanges(processed_peak)) != 150
table(is.smaller)

processed_peak <- processed_peak[!is.agap & !is.smaller,]

saveRDS(processed_peak, file='/collignon/Tawdarous/atac_seq_lab/results/csaw_result/processed_peak.rds') #save the windows


#-------------------------Export coverage tracks for IGV------------------------
colData(processed_peak)$totals <- colData(windows)$totals
gr <- rowRanges(processed_peak)

CPM <-  csaw::calculateCPM(processed_peak, log=FALSE)

sample_names <- colData(processed_peak)$sample_names
colnames(processed_peak) <- sample_names
colnames(CPM) <- colnames(processed_peak)
saveRDS(CPM, file='/collignon/Tawdarous/atac_seq_lab/results/csaw_result/CPM.rds') #save the windows

# Ensure the output directory exists
output_dir <- "/collignon/Tawdarous/atac_seq_lab/results/csaw_result/bw"


# Iterate over each sample and export coverage tracks
for (cn in colnames(CPM)) {
  gr$score <- CPM[, cn]
  bigwig_file <- file.path(output_dir, sprintf('%s.bigWig', cn))
  print(paste("Exporting:", bigwig_file))  # Debugging line
  rtracklayer::export.bw(gr, bigwig_file, format = "bigWig")
}

print("Export completed.")



#Read enrichment plot
fracLib <- matrix(0, nrow = nrow(processed_peak), ncol = ncol(processed_peak),
                  dimnames = dimnames(processed_peak))
for (cn in colnames(processed_peak)) {
  counts <- as.numeric(assays(processed_peak)$counts[,cn])
  counts <- sort(counts, decreasing=FALSE)
  fracLib[,cn] <- cumsum(counts)/sum(counts)
}
fracGenome <- seq(0, 1, length.out=nrow(processed_peak)) # since all windows are of the same size


data <- tibble(x=fracGenome) |>
  bind_cols(as_tibble(fracLib)) |>
  pivot_longer(-x) |>
  mutate(SampleGroup=str_replace(name, "_rep\\d+$", ""), 
         Replicate=str_replace(name, ".*_rep(\\d+)$", "R\\1"))

ggplot(data = data,
       mapping = aes(x=x, y=value, colour=SampleGroup, linetype=Replicate)) + 
  rasterize( geom_line() ) + 
  geom_abline(slope=1, intercept=0, linetype="dotted", colour="black") +
  scale_color_brewer(palette = "Set1") +
  labs(x="Fraction Genome", y="Fraction Library") +
  theme_bw(base_size=15)










# Define transcription start sites (TSS)
tss <- genes(TxDb.Mmusculus.UCSC.mm39.knownGene)
tss <- resize(tss, width = 1, fix = "start")
overlapTSS <- rowRanges(processed_peak) %over% tss
table(overlapTSS)

# Add overlap information to rowData
rowData(processed_peak)$overlapTSS <- overlapTSS

# Calculate logCPM and prepare data frame
colData(processed_peak)$totals <- colData(windows)$totals

logCPM <- csaw::calculateCPM(processed_peak, log = TRUE)

# Ensure the column names are assigned
colnames(logCPM) <- colnames(processed_peak)
colnames(logCPM) <- colnames(processed_peak)

# Verify the column names
print(colnames(logCPM))

# Merge cell types and calculate mean logCPM
ct <- split(colnames(logCPM), colData(processed_peak)$CellType)

df <- vapply(ct, function(x) {
  rowMeans(logCPM[, x, drop = FALSE])
}, rep(0, nrow(processed_peak))) %>%
  as_tibble() %>%
  mutate(overlapTSS = rowData(processed_peak)$overlapTSS,
         GC = rowData(processed_peak)$`G|C`)

# Define the plotMA function
plotMA <- function(df, s1, s2) {
  ggplot(df, aes_string(x = paste0(s1, "+", s2), y = paste0(s1, "-", s2))) +
    geom_hex(bins = 100, aes(fill = after_stat(density)^(1/16)), show.legend = FALSE) +
    geom_point(data = df[df$overlapTSS, ], aes_string(x = paste0(s1, "+", s2), y = paste0(s1, "-", s2)), col = 'red', pch = '.') +
    geom_smooth(data = df[df$overlapTSS, ], aes_string(x = paste0(s1, "+", s2), y = paste0(s1, "-", s2)), col = 'plum', se = FALSE) +
    geom_hline(yintercept = 0, lty = 2) +
    labs(x = "A", y = "M", title = sprintf("%s vs. %s", s1, s2)) +
    theme_bw(base_size = 15)
}

# Define combinations for comparison
combi <- list(c('t06_BMP+WNT', 't06_WNT'), c('t06_BMP+WNT', 't24_BMP+WNT-06h_LI'), c('t24_BMP+WNT-06h_L', 't24_WNT-06h_L'))

# Plot and combine the plots
gl <- lapply(combi, function(x) plotMA(df, x[1], x[2]))

patchwork::wrap_plots(gl, ncol = 2)



