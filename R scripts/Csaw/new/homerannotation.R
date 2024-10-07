#t06_BMP.WNT_vs_t06_WNT:
#setwd("/collignon/Tawdarous/atac_seq_lab/results/csaw_result/new_csaw/comparaison_paper/t06_BMP.WNT_vs_t06_WNT/homer")


# t24_BMP+WNT_L VS   t06_BMP+WNT 
setwd("/collignon/Tawdarous/atac_seq_lab/results/csaw_result/new_csaw/comparaison_paper/t24_BMP.WNT_06h_L_VS_t06_BMP.WNT/homer")

#t24_BMP.WNT_06h_LI  vs  t24_BMP.WNT_06h_L
setwd("/collignon/Tawdarous/atac_seq_lab/results/csaw_result/new_csaw/comparaison_paper/t24_BMP.WNT_06h_LI_vs_t24_BMP.WNT_06h_L/homer")


#t24_BMP.WNT_06h_LX1  vs  t24_BMP.WNT_06h_L
setwd("/collignon/Tawdarous/atac_seq_lab/results/csaw_result/new_csaw/comparaison_paper/t24_BMP.WNT_06h_LX1_vs_t24_BMP.WNT_06h_L/homer")


#t24_BMP+WNT_L vs t24_WNT_L

setwd("/collignon/Tawdarous/atac_seq_lab/results/csaw_result/new_csaw/comparaison_paper/t24_BMP.WNT_06h_L_vs_t24_WNT_06h_L/homer")


annotatePeaks.pl homer_down_m1.txt mm39 > homer_down_m1_annotated.txt
annotatePeaks.pl homer_down_m2.txt mm39 > homer_down_m2_annotated.txt
annotatePeaks.pl homer_down_m3.txt mm39 >homer_down_m3_annotated.txt
annotatePeaks.pl homer_down_m4.txt mm39 > homer_down_m4_annotated.txt
annotatePeaks.pl homer_up_m1.txt mm39 > homer_up_m1_annoteted.txt
annotatePeaks.pl homer_up_m2.txt mm39 > homer_up_m2_annotated.txt
annotatePeaks.pl homer_up_m3.txt mm39 > homer_up_m3_annotated.txt
annotatePeaks.pl homer_up_m4.txt mm39 > homer_up_m4_annotated.txt




homer_up_m1_annotated <- read.delim("homer_up_m1_annoteted.txt", stringsAsFactors = FALSE)
homer_down_m1_annotated <- read.delim("homer_down_m1_annotated.txt", stringsAsFactors = FALSE)

homer_up_m2_annotated <- read.delim("homer_up_m2_annotated.txt", stringsAsFactors = FALSE)
homer_down_m2_annotated <- read.delim("homer_down_m2_annotated.txt", stringsAsFactors = FALSE)


homer_up_m3_annotated <- read.delim("homer_up_m3_annotated.txt", stringsAsFactors = FALSE)
homer_down_m3_annotated <- read.delim("homer_down_m3_annotated.txt", stringsAsFactors = FALSE)


homer_up_m4_annotated <- read.delim("homer_up_m4_annotated.txt", stringsAsFactors = FALSE)
homer_down_m4_annotated <- read.delim("homer_down_m4_annotated.txt", stringsAsFactors = FALSE)


#______HELPER FUNCTION_______
process_homer_annotation <- function(homer_annotation, output_name_prefix) {
  # Dynamically get the column name for PeakID
  peakid_col <- names(homer_annotation)[grepl("PeakID", names(homer_annotation))][1]
  
  # Step 1: Convert the data frame to a tibble and rename columns dynamically
  homer_tidy <- homer_annotation %>%
    as_tibble() %>%
    rename(
      seqnames = !!sym(peakid_col),  # Use the dynamically identified column for seqnames
      start = Start, 
      end = End
    )
  
  # Step 2: Create GRanges object
  homer_gr <- GRanges(
    seqnames = homer_tidy$seqnames,
    ranges = IRanges(start = homer_tidy$start, end = homer_tidy$end),
    strand = homer_tidy$Strand
  )
  
  # Assign the outputs to the global environment with the specified names
  assign(paste0(output_name_prefix, "_tidy"), homer_tidy, envir = .GlobalEnv)
  assign(paste0(output_name_prefix, "_gr"), homer_gr, envir = .GlobalEnv)
}

# This will create two objects:
# - homer_up_m1_tidy: the tidied data
# - homer_up_m1_gr: the GRanges object
process_homer_annotation(homer_up_m1_annotated, "homer_up_m1")
process_homer_annotation(homer_down_m1_annotated, "homer_down_m1")
process_homer_annotation(homer_up_m2_annotated, "homer_up_m2")
process_homer_annotation(homer_down_m2_annotated, "homer_down_m2")

process_homer_annotation(homer_up_m3_annotated, "homer_up_m3")
process_homer_annotation(homer_down_m3_annotated, "homer_down_m3")

process_homer_annotation(homer_up_m4_annotated, "homer_up_m4")
process_homer_annotation(homer_down_m4_annotated, "homer_down_m4")






# Create a GRanges object
metadata_columns <- setdiff(names(homer_down_m1_tidy), c("seqnames", "start", "end", "Strand"))
mcols(homer_down_m1_gr) <- homer_down_m1_tidy[, metadata_columns]
homer_down_m1_gr


# Distribution of distances to TSS
hist(homer_down_m1_tidy$Distance.to.TSS, breaks = 50, main = "Distribution of Distance to TSS")

# Function to count peaks and plot top 20 genes by number of peaks
# Load necessary libraries
library(dplyr)

# Function to plot top 20 upregulated and downregulated genes by number of peaks
plot_top_genes_by_peaks <- function(homer_up_tidy, homer_down_tidy, prefix = "M1") {
  
  # Step 1: Process upregulated genes
  gene_peak_counts_up <- table(homer_up_tidy$Gene.Name)
  top_genes_up <- sort(gene_peak_counts_up, decreasing = TRUE)[1:20]
  
  # Plot top 20 upregulated genes
  barplot(
    top_genes_up, 
    las = 2, 
    main = paste("Top 20 Upregulated Genes by Number of Peaks (", prefix, ")", sep = ""),
    col = "skyblue", 
    cex.names = 0.7
  )
  
  # Step 2: Process downregulated genes
  gene_peak_counts_down <- table(homer_down_tidy$Gene.Name)
  top_genes_down <- sort(gene_peak_counts_down, decreasing = TRUE)[1:20]
  
  # Plot top 20 downregulated genes
  barplot(
    top_genes_down, 
    las = 2, 
    main = paste("Top 20 Downregulated Genes by Number of Peaks (", prefix, ")", sep = ""),
    col = "salmon", 
    cex.names = 0.7
  )
  
  # Optionally return the top genes for further inspection
  return(list(top_genes_up = top_genes_up, top_genes_down = top_genes_down))
}

# Example usage for M1
results_m1 <- plot_top_genes_by_peaks(homer_up_m1_tidy, homer_down_m1_tidy, prefix = "M1")

# Example usage for M2
results_m2 <- plot_top_genes_by_peaks(homer_up_m2_tidy, homer_down_m2_tidy, prefix = "M2")

# Example usage for M3
results_m3 <- plot_top_genes_by_peaks(homer_up_m3_tidy, homer_down_m3_tidy, prefix = "M3")

# Example usage for M4
results_m4 <- plot_top_genes_by_peaks(homer_up_m4_tidy, homer_down_m4_tidy, prefix = "M4")



perform_go_enrichment <- function(homer_tidy, ont = "BP", top_n = 10, org_db = org.Mm.eg.db) {
  
  # Step 1: Extract unique genes with peaks
  genes_with_peaks <- unique(homer_tidy$Gene.Name)
  
  # Step 2: Perform GO enrichment analysis
  ego <- enrichGO(gene = genes_with_peaks,
                  OrgDb = org_db,
                  keyType = "SYMBOL",
                  ont = ont,
                  pAdjustMethod = "BH",
                  qvalueCutoff = 0.05,
                  readable = TRUE)
  
  # Step 3: Plot top GO terms
  plot <- dotplot(ego, showCategory = top_n) +
    ggtitle(paste("Top", top_n, "GO terms in", ont, "ontology"))
  
  # Print the plot to ensure it displays
  print(plot)
  
  # Return the enrichment result for further inspection if needed
  return(ego)
}


ego_m1_up <- perform_go_enrichment(homer_up_m1_tidy, ont = "BP", top_n = 10)
ego_m1_down <- perform_go_enrichment(homer_down_m1_tidy, ont = "BP", top_n = 10)

ego_m2_up <- perform_go_enrichment(homer_up_m2_tidy, ont = "BP", top_n = 10)
ego_m2_down <- perform_go_enrichment(homer_down_m2_tidy, ont = "BP", top_n = 10)


ego_m3_up <- perform_go_enrichment(homer_up_m3_tidy, ont = "BP", top_n = 10)
ego_m3_down <- perform_go_enrichment(homer_down_m3_tidy, ont = "BP", top_n = 10)

ego_m4_up <- perform_go_enrichment(homer_up_m4_tidy, ont = "BP", top_n = 10)
ego_m4_down <- perform_go_enrichment(homer_down_m4_tidy, ont = "BP", top_n = 10)


#-------------



pathway_of_interest <- "positive regulation of nervous system development"
selected_pathway <- ego_m1_up %>%
  dplyr::filter(Description == pathway_of_interest) # or use the column name for the pathway ID

genes_in_pathway <- selected_pathway$geneID # Replace with the actual column name where genes are stored
genes_in_pathway_vector <- strsplit(genes_in_pathway, "/")[[1]]
print(genes_in_pathway_vector)

#-----------mol signature----------

m_df <- msigdbr(species = "Mus musculus", category = "H")
hallmark_gene_sets <- split(m_df$gene_symbol, m_df$gs_name)

#______________
# Create a GRanges object

metadata_columns <- setdiff(names(homer_up_m1_tidy), c("seqnames", "start", "end", "Strand"))
mcols(homer_up_m1_gr) <- homer_up_m1_tidy[, metadata_columns]

metadata_columns <- setdiff(names(homer_down_m1_tidy), c("seqnames", "start", "end", "Strand"))
mcols(homer_down_m1_gr) <- homer_down_m1_tidy[, metadata_columns]


metadata_columns <- setdiff(names(homer_up_m2_tidy), c("seqnames", "start", "end", "Strand"))
mcols(homer_up_m2_gr) <- homer_up_m2_tidy[, metadata_columns]

metadata_columns <- setdiff(names(homer_down_m2_tidy), c("seqnames", "start", "end", "Strand"))
mcols(homer_down_m2_gr) <- homer_down_m2_tidy[, metadata_columns]



metadata_columns <- setdiff(names(homer_up_m3_tidy), c("seqnames", "start", "end", "Strand"))
mcols(homer_up_m3_gr) <- homer_up_m3_tidy[, metadata_columns]

metadata_columns <- setdiff(names(homer_down_m3_tidy), c("seqnames", "start", "end", "Strand"))
mcols(homer_down_m3_gr) <- homer_down_m3_tidy[, metadata_columns]



metadata_columns <- setdiff(names(homer__m4_tidy), c("seqnames", "start", "end", "Strand"))
mcols(homer_up_m4_gr) <- homer_up_m4_tidy[, metadata_columns]

metadata_columns <- setdiff(names(homer_down_m4_tidy), c("seqnames", "start", "end", "Strand"))
mcols(homer_down_m4_gr) <- homer_down_m4_tidy[, metadata_columns]


#__________
gene_scores <- setNames(homer_down_m4_gr$Peak.Score, homer_down_m4_gr$Gene.Name)

gene_scores <- gene_scores[!duplicated(names(gene_scores))]
gene_scores <- sort(gene_scores, decreasing = TRUE)


fgseaRes <- fgsea(pathways = hallmark_gene_sets, 
                  stats = gene_scores,
                  minSize = 1,
                  maxSize = 500,
                  )

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(pval)

# Print top enriched pathways
print(fgseaResTidy %>% select(pathway, pval, padj, NES) %>% head(n = 10))

ggplot(fgseaResTidy %>% head(n = 10), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = pval < 0.05)) +
  coord_flip() +
  labs(x = "Pathway", y = "Normalized Enrichment Score",
       title = "Top 10 up Pathways") +
  theme_minimal()

pathway_of_interest <- fgseaResTidy$pathway[2]  # First pathway in the results

# Plot the enrichment plot for this pathway
plotEnrichment(hallmark_gene_sets[[pathway_of_interest]], gene_scores) +
  labs(title = pathway_of_interest)

gene_to_pathway <- fgseaRes %>%
  filter(padj < 0.05) %>%
  unnest(leadingEdge) %>%
  group_by(leadingEdge) %>%
  slice_min(order_by = padj, n = 1) %>%
  select(gene = leadingEdge, top_hallmark_pathway = pathway)

homer_gr$top_hallmark_pathway <- gene_to_pathway$top_hallmark_pathway[match(homer_gr$Nearest.PromoterID, gene_to_pathway$gene)]


