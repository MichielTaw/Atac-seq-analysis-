setwd("/collignon/Tawdarous/atac_seq_lab/results/csaw_result/new_csaw/comparaison_paper/t06_BMP.WNT_vs_t06_WNT")

significant_peaks_m1 <- readRDS("/collignon/Tawdarous/atac_seq_lab/results/csaw_result/new_csaw/comparaison_paper/t06_BMP.WNT_vs_t06_WNT/significant_peaks_m1.rds")
significant_peaks_m2 <- readRDS("/collignon/Tawdarous/atac_seq_lab/results/csaw_result/new_csaw/comparaison_paper/t06_BMP.WNT_vs_t06_WNT/significant_peaks_m2.rds")
significant_peaks_m3 <- readRDS("/collignon/Tawdarous/atac_seq_lab/results/csaw_result/new_csaw/comparaison_paper/t06_BMP.WNT_vs_t06_WNT/significant_peaks_m3.rds")
significant_peaks_m4 <- readRDS("/collignon/Tawdarous/atac_seq_lab/results/csaw_result/new_csaw/comparaison_paper/t06_BMP.WNT_vs_t06_WNT/significant_peaks_m4.rds")

up_significant_peaks_m1 <- significant_peaks_m1[significant_peaks_m1 $direction=="up"]
down_significant_peaks_m1 <- significant_peaks_m1[significant_peaks_m1 $direction=="down"]

up_significant_peaks_m2 <- significant_peaks_m2[significant_peaks_m2 $direction=="up"]
down_significant_peaks_m2 <- significant_peaks_m2[significant_peaks_m2 $direction=="down"]

up_significant_peaks_m3 <- significant_peaks_m3[significant_peaks_m3 $direction=="up"]
down_significant_peaks_m3 <- significant_peaks_m3[significant_peaks_m3 $direction=="down"]

up_significant_peaks_m4 <- significant_peaks_m4[significant_peaks_m4 $direction=="up"]
down_significant_peaks_m4 <- significant_peaks_m4[significant_peaks_m4 $direction=="down"]


output <- "/collignon/Tawdarous/atac_seq_lab/results/csaw_result/new_csaw/comparaison_paper/t06_BMP.WNT_vs_t06_WNT"

# 1. Create a standard BED file
# First, let's add a name column to our GRanges object
up_significant_peaks_m1$name <- paste0("peak_", 1:length(up_significant_peaks_m1))
down_significant_peaks_m1$name <- paste0("peak_", 1:length(down_significant_peaks_m1))
# Export as BED
#export.bed(significant_peaks_m1, "significant_peaks_m1.bed")


up_significant_peaks_m1_homer_df <- data.frame(
  chr = seqnames(up_significant_peaks_m1),
  start = start(up_significant_peaks_m1) - 1,  # BED is 0-based, GRanges is 1-based
  end = end(up_significant_peaks_m1),
  name = up_significant_peaks_m1$name,
  score = up_significant_peaks_m1$logFC,  # Using logFC as score, you can change this
  strand = strand(up_significant_peaks_m1)
)

down_significant_peaks_m1_homer_df <- data.frame(
  chr = seqnames(down_significant_peaks_m1),
  start = start(down_significant_peaks_m1) - 1,  # BED is 0-based, GRanges is 1-based
  end = end(down_significant_peaks_m1),
  name = down_significant_peaks_m1$name,
  score = down_significant_peaks_m1$logFC,  # Using logFC as score, you can change this
  strand = strand(down_significant_peaks_m1)
)

# Replace '*' with '.' for strand (HOMER preference)
up_significant_peaks_m1_homer_df$strand <- ifelse(up_significant_peaks_m1_homer_df$strand == "*", ".", as.character(up_significant_peaks_m1_homer_df$strand))
down_significant_peaks_m1_homer_df$strand <- ifelse(down_significant_peaks_m1_homer_df$strand == "*", ".", as.character(down_significant_peaks_m1_homer_df$strand))



# If you want to include additional columns for HOMER (optional)
up_significant_peaks_m1_homer_df <- up_significant_peaks_m1_homer_df %>%
  mutate(
    logCPM = up_significant_peaks_m1$logCPM,
    F_statistic = up_significant_peaks_m1$F,
    PValue = up_significant_peaks_m1$PValue,
    FDR = up_significant_peaks_m1$FDR
  )

write.table(up_significant_peaks_m1_homer_df, "up_significant_peaks_m1_homer_df.txt", 
            sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

down_significant_peaks_m1_homer_df <- down_significant_peaks_m1_homer_df %>%
  mutate(
    logCPM = down_significant_peaks_m1$logCPM,
    F_statistic = down_significant_peaks_m1$F,
    PValue = down_significant_peaks_m1$PValue,
    FDR = down_significant_peaks_m1$FDR
  )


write.table(down_significant_peaks_m1_homer_df, "down_significant_peaks_m1_homer_df.txt", 
            sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

up_significant_peaks_m1_homer_annotation <- read.delim("up_significant_peaks_m1_homer_df_annotated.txt", stringsAsFactors = FALSE)
down_significant_peaks_m1_homer_annotation <- read.delim("down_significant_peaks_m1_homer_df_annotated_annotated.txt", stringsAsFactors = FALSE)



head(up_significant_peaks_m1_homer_annotation)
head(down_significant_peaks_m1_homer_annotation)

colnames(up_significant_peaks_m1_homer_annotation) <- make.names(colnames(up_significant_peaks_m1_homer_annotation))
colnames(down_significant_peaks_m1_homer_annotation) <- make.names(colnames(down_significant_peaks_m1_homer_annotation))



up_tidy_m1 <- up_significant_peaks_m1_homer_annotation %>%
  as_tibble() %>%
  rename(seqnames = PeakID..cmd.annotatePeaks.pl.up_significant_peaks_m1_homer_df.txt.mm39., start = Start, end = End)

down_tidy_m1 <- down_significant_peaks_m1_homer_annotation %>%
  as_tibble() %>%
  rename(seqnames = PeakID..cmd.annotatePeaks.pl.down_significant_peaks_m1_homer_df.txt.mm39., start = Start, end = End)


up_tidy_m1_gr <- GRanges(
  seqnames = up_tidy_m1$seqnames,
  ranges = IRanges(start = up_tidy_m1$start, end = up_tidy_m1$end),
  strand = up_tidy_m1$Strand
)

down_tidy_m1_gr <- GRanges(
  seqnames = down_tidy_m1$seqnames,
  ranges = IRanges(start = down_tidy_m1$start, end = down_tidy_m1$end),
  strand = down_tidy_m1$Strand
)

metadata_up_tidy_m1 <- setdiff(names(up_tidy_m1), c("seqnames", "start", "end", "Strand"))
mcols(up_tidy_m1_gr) <- up_tidy_m1[, metadata_up_tidy_m1]
up_tidy_m1_gr


metadata_down_tidy_m1 <- setdiff(names(down_tidy_m1), c("seqnames", "start", "end", "Strand"))
mcols(down_tidy_m1_gr) <- down_tidy_m1[, metadata_down_tidy_m1]

# Basic summary of annotations
table(up_tidy_m1$Annotation)

# Distribution of distances to TSS
hist(up_tidy_m1$Distance.to.TSS, breaks = 50, main = "Distribution of Distance to TSS")
hist(down_tidy_m1$Distance.to.TSS, breaks = 50, main = "Distribution of Distance to TSS")

# Genes with the most peaks
up_gene_peak_counts <- table(up_tidy_m1$Gene.Name)
up_top_genes <- sort(up_gene_peak_counts, decreasing = TRUE)[1:20]
barplot(up_top_genes, las = 2, main = "Top 20 Genes up_top_genes")

down_gene_peak_counts <- table(down_tidy_m1$Gene.Name)
down_top_genes <- sort(down_gene_peak_counts, decreasing = TRUE)[1:20]
barplot(down_top_genes, las = 2, main = "Top 20 Genes down_top_genes")


genes_with_peaks <- unique(up_tidy_m1$Gene.Name)
ego <- enrichGO(gene = genes_with_peaks,
                OrgDb = org.Mm.eg.db,
                keyType = "SYMBOL",
                ont = "BP",
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05,
                readable = TRUE)

# Plot top GO terms
dotplot(ego, showCategory = 10)

m_df <- msigdbr(species = "Mus musculus", category = "H")
hallmark_gene_sets <- split(m_df$gene_symbol, m_df$gs_name)


nPerm <- 100000


up_gene_scores <- setNames(up_tidy_m1_gr$Peak.Score, up_tidy_m1_gr$Gene.Name)
down_gene_scores <- setNames(down_tidy_m1_gr$Peak.Score, down_tidy_m1_gr$Gene.Name)

up_gene_scores <- up_gene_scores[!duplicated(names(up_gene_scores))]
up_gene_scores <- sort(up_gene_scores, decreasing = TRUE)

down_gene_scores <- down_gene_scores[!duplicated(names(down_gene_scores))]
down_gene_scores <- sort(down_gene_scores, decreasing = TRUE)

fgseaRes_down <- fgsea(pathways = hallmark_gene_sets, 
                       stats = down_gene_scores,
                       minSize = 15,
                       maxSize = 500,
                       nPermSimple = nPerm)
fgseaRes_up <- fgsea(pathways = hallmark_gene_sets, 
                       stats = up_gene_scores,
                       minSize = 15,
                       maxSize = 500,
                       nPermSimple = nPerm)
fgseaResTidy_up <- fgseaRes_up %>%
  as_tibble() %>%
  arrange(pval)

fgseaResTidy_down <- fgseaRes_down %>%
  as_tibble() %>%
  arrange(pval)

ggplot(fgseaResTidy_down %>% head(n = 5), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = pval < 0.05)) +
  coord_flip() +
  labs(x = "Pathway", y = "Normalized Enrichment Score",
       title = "Top 5 Hallmark Pathways t06 WNT") +
  theme_minimal()

pathway_of_interest <- fgseaResTidy_up$pathway[2]  # First pathway in the results

plotEnrichment(hallmark_gene_sets[[pathway_of_interest]], gene_scores) +
  labs(title = pathway_of_interest)

gene_to_pathway <- fgseaRes_up %>%
  unnest(leadingEdge) %>%
  group_by(leadingEdge) %>%
  slice_min(order_by = padj, n = 1) %>%
  select(gene = leadingEdge, top_hallmark_pathway = pathway)

up_tidy_m1_gr$top_hallmark_pathway <- gene_to_pathway$top_hallmark_pathway[match(up_tidy_m1_gr$Nearest.PromoterID, gene_to_pathway$gene)]
