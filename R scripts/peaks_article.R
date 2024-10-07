setwd("/collignon/Tawdarous/atac_seq_lab/results/peak_analysis/macs3/Homer/replicated_peaks_annotate")

#Function to transform gene names to uppercase
transform_gene_names <- function(data) {
  data %>%
    mutate(Gene.Name = toupper(Gene.Name))
}

# Function to calculate mean peak score for specified genes
calculate_mean_peak_score <- function(data, gene_list) {
  data %>%
    filter(Gene.Name %in% gene_list) %>%
    group_by(Gene.Name) %>%
    summarize(mean_peak_score = mean(Peak.Score, na.rm = TRUE))
}

# List of genes to analyze (already in uppercase)
genes_of_interest <- c("TAL1", "LYL1", "LMO2", "KDR", "HHEX", "SOX17", "PRDM14", "PROM1", "DND1", "DPPA3",
                       "MSXI", "HAND1", "FOXF1", "CDX2", "BMP4", "SNAI", "RSP03", "MESP2", "MESP1", "CDH2",
                       "WNT3", "FGF8", "MIXL1", "FST", "T", "SFN", "KRT18", "KRTB", "CLDN4", "NDRG2",
                       "NDRG1", "CLDN7", "CLDN6", "CLDN3", "CLONG", "TDGF1", "POU5F1", "NODAL", "NANDG",
                       "ID1", "UTF1", "SOX2", "POU3F1", "OTX2", "CDH1")
genes_of_interest <- c("POU5F1", "CRIPTO", "NODAL", "SOX2", "OTX2",
                       "T", "AXIN2", "WNT3", "TFAP2A", "ID1", "NANOG", "FST", 
                       "HAND1", "BMP2", "GATA4", "CDX2", "TBX4")
# Read and transform dataframes
t06_BMP_WNT <- read.delim("t06_BMP+WNT.replicated_broadPeak_chrom_mapped_annotated.txt", stringsAsFactors = FALSE) %>% transform_gene_names()
t06_WNT <- read.delim("t06_WNT.replicated_broadPeak_chrom_mapped_annotated.txt", stringsAsFactors = FALSE) %>% transform_gene_names()
t24_BMP_WNT_06h_LI <- read.delim("t24_BMP+WNT-06h_LI.replicated_broadPeak_chrom_mapped_annotated.txt", stringsAsFactors = FALSE) %>% transform_gene_names()
t24_BMP_WNT_06h_L <- read.delim("t24_BMP+WNT-06h_L.replicated_broadPeak_chrom_mapped_annotated.txt", stringsAsFactors = FALSE) %>% transform_gene_names()
t24_BMP_WNT_06h_LX1 <- read.delim("t24_BMP+WNT-06h_LX1.replicated_broadPeak_chrom_mapped_annotated.txt", stringsAsFactors = FALSE) %>% transform_gene_names()
t24_WNT_06h_L <- read.delim("t24_WNT-06h_L.replicated_broadPeak_chrom_mapped_annotated.txt", stringsAsFactors = FALSE) %>% transform_gene_names()

# List of dataframes
dataframes <- list(t06_BMP_WNT, t06_WNT, t24_BMP_WNT_06h_LI, t24_BMP_WNT_06h_L, t24_BMP_WNT_06h_LX1, t24_WNT_06h_L)
names(dataframes) <- c("t06_BMP_WNT", "t06_WNT", "t24_BMP_WNT_06h_LI", "t24_BMP_WNT_06h_L", "t24_BMP_WNT_06h_LX1", "t24_WNT_06h_L")


dataframes <- list(t24_BMP_WNT_06h_LI, t24_BMP_WNT_06h_L)
names(dataframes) <- c("t24_BMP_WNT_06h_LI", "t24_BMP_WNT_06h_L")

# Calculate mean peak scores for each dataset
results <- lapply(names(dataframes), function(name) {
  df <- dataframes[[name]]
  result <- calculate_mean_peak_score(df, genes_of_interest)
  result$condition <- name
  return(result)
})

# Combine results
all_results <- bind_rows(results)

# Reshape the data for easier comparison
final_results <- all_results %>%
  pivot_wider(names_from = condition, values_from = mean_peak_score)

# Display the results
print(final_results)


long_data <- final_results %>%
  pivot_longer(cols = -Gene.Name, names_to = "Condition", values_to = "Expression")

# Calculate percent expressed (assuming Expression > 0 means the gene is expressed)
percent_expressed <- long_data %>%
  group_by(Gene.Name, Condition) %>%
  summarize(Percent_Expressed = mean(Expression > 0) * 100,
            Avg_Expression = mean(Expression, na.rm = TRUE))

# Define gene groups (adjust as needed based on your data)
gene_groups <- list(
  c("TAL1", "LYL1", "LMO2", "KDR", "HHEX"),
  c("SOX17", "PRDM14", "PRDM1", "DND1", "DPPA3"),
  c("MSX1", "HAND1", "FOXF1", "CDX2", "BMP4"),
  c("SNAI1", "RSPO3", "MESP2", "MESP1", "CDH2"),
  c("WNT3", "FGF8", "MIXL1", "FST", "T"),
  c("TFAP2A", "SFN", "KRT18", "KRT8", "CLDN4"),
  c("NDRG2", "NDRG1", "CLDN7", "CLDN6", "CLDN3"),
  c("TDGF1", "POU5F1", "NODAL", "NANOG", "ID1"),
  c("UTF1", "SOX2", "POU3F1", "OTX2", "CDH1")
)
gene_groups <- list(c("POU5F1", "CRIPTO", "NODAL", "SOX2", "OTX2"),
                      c( "T", "AXIN2", "WNT3", "TFAP2A", "ID1", "NANOG", "FST"), 
                       c("HAND1", "BMP2", "GATA4", "CDX2", "TBX4")
                        )
# Add group information to the data
percent_expressed$Group <- NA
for (i in seq_along(gene_groups)) {
  percent_expressed$Group[percent_expressed$Gene.Name %in% gene_groups[[i]]] <- i
}

# Order genes based on groups
percent_expressed$Gene.Name <- factor(percent_expressed$Gene.Name, 
                                      levels = unlist(gene_groups))

# Create the plot
ggplot(percent_expressed, aes(x = Condition, y = Gene.Name)) +
  geom_point(aes(size = Percent_Expressed, color = Avg_Expression)) +
  scale_size_continuous(range = c(0, 5), breaks = c(0, 25, 50, 75, 100)) +
  scale_color_gradient(low = "lightblue", high = "darkblue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_blank()) +
  labs(x = "", y = "", size = "Percent Expressed", color = "Average Expression") +
  geom_hline(yintercept = seq(5.5, 40.5, by = 5), color = "grey80", size = 0.5)


ggsave("atac_seq_peak_analysis_dotplot_with_separators.png", width = 15, height = 20, dpi = 300)







final_results <- read.csv("atac_seq_mean_peak_scores.csv")

# Reshape the data from wide to long format
long_data <- final_results %>%
  pivot_longer(cols = -Gene.Name, names_to = "Condition", values_to = "Mean_Peak_Score")

# Calculate mean peak score and count of peaks (non-NA values)
peak_data <- long_data %>%
  group_by(Gene.Name, Condition) %>%
  summarize(
    Mean_Peak_Score = mean(Mean_Peak_Score, na.rm = TRUE),
    Peak_Count = sum(!is.na(Mean_Peak_Score))
  )

# Define gene groups
gene_groups <- list(
  c("TAL1", "LYL1", "LMO2", "KDR", "HHEX"),
  c("SOX17", "PRDM14", "PRDM1", "DND1", "DPPA3"),
  c("MSX1", "HAND1", "FOXF1", "CDX2", "BMP4"),
  c("SNAI1", "RSPO3", "MESP2", "MESP1", "CDH2"),
  c("WNT3", "FGF8", "MIXL1", "FST", "T"),
  c("TFAP2A", "SFN", "KRT18", "KRT8", "CLDN4"),
  c("NDRG2", "NDRG1", "CLDN7", "CLDN6", "CLDN3"),
  c("TDGF1", "POU5F1", "NODAL", "NANOG", "ID1"),
  c("UTF1", "SOX2", "POU3F1", "OTX2", "CDH1")
)

# Flatten the gene groups list to create the gene order
gene_order <- unlist(gene_groups)

# Order genes and conditions
peak_data$Gene.Name <- factor(peak_data$Gene.Name, levels = rev(gene_order))
peak_data$Condition <- factor(peak_data$Condition, 
                              levels = c("t24_BMP_WNT_06h_LI", "t24_BMP_WNT_06h_L"))

# Calculate positions for separating lines
separator_positions <- cumsum(sapply(gene_groups, length))
separator_positions <- separator_positions[-length(separator_positions)] + 0.5

# Create the plot
ggplot(peak_data, aes(x = Condition, y = Gene.Name)) +
  geom_point(aes(size = Peak_Count, color = Mean_Peak_Score)) +
  geom_hline(yintercept = separator_positions, color = "grey80", size = 0.5) +
  scale_size_continuous(range = c(1, 8), limits = c(0, max(peak_data$Peak_Count)), 
                        breaks = pretty_breaks(n = 5)) +
  scale_color_gradient(low = "lightblue", high = "darkblue", 
                       limits = c(0, max(peak_data$Mean_Peak_Score)),
                       labels = comma) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 8),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 16)
  ) +
  labs(
    x = "", 
    y = "", 
    size = "Number of Peaks", 
    color = "Mean Peak Score",
    title = "ATAC-seq Peak Analysis Across Conditions"
  ) +
  guides(
    size = guide_legend(override.aes = list(color = "black")),
    color = guide_colorbar(barwidth = 1, barheight = 10)
  )



ggsave("atac_seq_peak_analysis_dotplot_with_separators.png", width = 15, height = 20, dpi = 300)

