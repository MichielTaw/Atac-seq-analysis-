# Load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)
gene_groups <- list(
  c("SLC7A3", "STXBP6", "SLC35D3", "SNCA", "UCHL1", "ENPP3", "SYT11", "DNMT3A", "LYPD6B","GM42637"),#epiblast
  c("CDKN1A", "PDLIM4", "EPHA1", "AXIN2", "LGFBP3", "APELA", "NEFM", "IRF1", "WNT5B", "SALL3"), #PS
  c("DDX5", "UTP4", "EIF1A", "STX3", "MTX2", "TMX1", "EXOC5", "GNG10", "TNPO3", "MESP1"), #nascent mesoderm
  c("XKR4", "GM1992", "GM37381", "RP1", "SOX17", "GM37323", "MRPL15", "LYPLA1", "GM37988", "TCEA1"), #rostral neuroectoderm
  c("CD59A", "GPC4", "AMOT", "SFRP1", "SMIM14", "CER1", "LBH", "GGCT", "FAM89A", "KRT8"), #visceral endoderm
  c("CYBA", "RBP4", "NPL", "APOA4", "H19", "TTR", "APOA1", "APOC2", "LGALS2", "APOE"), #EXE endoderm
  c("GJB3", "CLDN4", "WNT6", "BMP4", "GM9", "AKAP17B", "PRSS50", "DUOXA2", "INPP1", "AA467197"),#exe ectoderm
  c("SERPINH1", "CRYAB", "PTH1R", "FST", "TAX1BP3", "F3", "HS3ST1", "SPARC", "SAT1", "SRGN")#parietal endoderm
)
# Define gene groups (as in your original code)
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


gene_groups <- list(
  c("SLC7A3", "STXBP6", "SLC35D3", "SNCA", "UCHL1", "ENPP3", "SYT11", "DNMT3A", "LYPD6B","GM42637"),
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
genes_of_interest <- unlist(gene_groups)

# Function to process homer annotated data
process_homer_data <- function(data, genes_of_interest) {
  # Get the actual column names
  cols <- colnames(data)
  
  # Find the column that likely contains gene names (assuming it has "Gene.Name" in its name)
  gene_col <- cols[grep("Gene.Name", cols, ignore.case = TRUE)]
  
  # Find the column that likely contains peak scores (assuming it has "Peak.Score" in its name)
  score_col <- cols[grep("Peak.Score", cols, ignore.case = TRUE)]
  
  if(length(gene_col) == 0 || length(score_col) == 0) {
    stop("Could not find appropriate columns for gene names or peak scores")
  }
  
  # Process the data
  data %>%
    select(all_of(c(gene_col, score_col))) %>%
    rename(Gene.Name = !!gene_col, Peak.Score = !!score_col) %>%
    mutate(Gene.Name = toupper(Gene.Name)) %>%
    filter(Gene.Name %in% toupper(genes_of_interest)) 
}

# Assuming you have a list of genes_of_interest to filter on
genes_of_interest <- c("GENE1", "GENE2", "GENE3")  # Replace with your actual genes of interest

# Assuming you have already processed `homer_up_m1_tidy` and `homer_down_m1_tidy`

# Add Regulation column to up and down data
up_data <- process_homer_data(homer_up_m3_tidy, genes_of_interest) %>%
  mutate(Regulation = "Upregulated")

down_data <- process_homer_data(homer_down_m3_tidy, genes_of_interest) %>%
  mutate(Regulation = "Downregulated")

# Combine all processed data
all_data <- bind_rows(up_data, down_data)

# Calculate mean peak score for each gene in each condition
peak_data <- all_data %>%
  group_by(Gene.Name, Regulation) %>%
  summarize(
    Mean_Peak_Score = mean(Peak.Score, na.rm = TRUE),
    Peak_Count = n()
  ) %>%
  ungroup()

# Order genes based on the provided gene order
peak_data$Gene.Name <- factor(peak_data$Gene.Name, levels = rev(genes_of_interest))

# Calculate positions for separating lines between groups
separator_positions <- cumsum(sapply(gene_groups, length))
separator_positions <- separator_positions[-length(separator_positions)] + 0.5

# Create the plot
p <- ggplot(peak_data, aes(x = Regulation, y = Gene.Name)) +
  geom_point(aes(size = Peak_Count, color = Mean_Peak_Score)) +
  geom_hline(yintercept = separator_positions, color = "grey80", size = 0.5) +
  scale_size_continuous(range = c(1, 8), limits = c(0, max(peak_data$Peak_Count)), 
                        breaks = scales::pretty_breaks(n = 5)) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", 
                        midpoint = median(peak_data$Mean_Peak_Score, na.rm = TRUE),
                        limits = c(0, max(peak_data$Mean_Peak_Score, na.rm = TRUE)),
                        labels = scales::comma) +
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
    title = "ATAC-seq Peak Analysis: Upregulated vs Downregulated Genes"
  ) +
  guides(
    size = guide_legend(override.aes = list(color = "black")),
    color = guide_colorbar(barwidth = 1, barheight = 10)
  )

# Display the plot
print(p)

# Save the plot
setwd("/collignon/Tawdarous/atac_seq_lab/results/csaw_result/new_csaw/comparaison_paper/t06_BMP.WNT_vs_t06_WNT/homer")

ggsave("atac_seq_peak_analysis_up_vs_down.png", plot = p, width = 12, height = 20, dpi = 300)


gene_groups <- list(
  c("SLC7A3", "STXBP6", "SLC35D3", "SNCA", "UCHL1", "ENPP3", "SYT11", "DNMT3A", "LYPD6B","GM42637"),
  c("CDKN1A", "PDLIM4", "EPHA1", "AXIN2", "LGFBP3", "APELA", "NEFM", "IRF1", "WNT5B", "SALL3"), 
  c("DDX5", "UTP4", "EIF1A", "STX3", "MTX2", "TMX1", "EXOC5", "GNG10", "TNPO3", "MESP1"),
  c("XKR4", "GM1992", "GM37381", "RP1", "SOX17", "GM37323", "MRPL15", "LYPLA1", "GM37988", "TCEA1")
)


# Function to process homer annotated data with an option to take absolute value
process_homer_data <- function(data, genes_of_interest, abs_peak_score = FALSE) {
  # Get the actual column names
  cols <- colnames(data)
  
  # Find the column that likely contains gene names (assuming it has "Gene.Name" in its name)
  gene_col <- cols[grep("Gene.Name", cols, ignore.case = TRUE)]
  
  # Find the column that likely contains peak scores (assuming it has "Peak.Score" in its name)
  score_col <- cols[grep("Peak.Score", cols, ignore.case = TRUE)]
  
  if(length(gene_col) == 0 || length(score_col) == 0) {
    stop("Could not find appropriate columns for gene names or peak scores")
  }
  
  # Process the data
  processed_data <- data %>%
    select(all_of(c(gene_col, score_col))) %>%
    rename(Gene.Name = !!gene_col, Peak.Score = !!score_col) %>%
    mutate(Gene.Name = toupper(Gene.Name)) %>%
    filter(Gene.Name %in% toupper(genes_of_interest))
  
  # Apply absolute value to Peak.Score if specified
  if (abs_peak_score) {
    processed_data <- processed_data %>%
      mutate(Peak.Score = abs(Peak.Score))
  }
  
  return(processed_data)
}

# Define your genes of interest (already done previously)
genes_of_interest <- unlist(gene_groups)

# Process upregulated data (no absolute value transformation)
up_data <- process_homer_data(homer_up_m1_annotated, genes_of_interest) %>%
  mutate(Regulation = "Upregulated")

# Process downregulated data with absolute value for Peak.Score
down_data <- process_homer_data(homer_down_m1_annotated, genes_of_interest, abs_peak_score = TRUE) %>%
  mutate(Regulation = "Downregulated")

# Combine all processed data
all_data <- bind_rows(up_data, down_data)

# Calculate mean peak score for each gene in each condition
peak_data <- all_data %>%
  group_by(Gene.Name, Regulation) %>%
  summarize(
    Mean_Peak_Score = mean(Peak.Score, na.rm = TRUE),
    Peak_Count = n()
  ) %>%
  ungroup()

# Order genes based on the provided gene order
peak_data$Gene.Name <- factor(peak_data$Gene.Name, levels = rev(genes_of_interest))

# Calculate positions for separating lines between groups
separator_positions <- cumsum(sapply(gene_groups, length))
separator_positions <- separator_positions[-length(separator_positions)] + 0.5

# Create the plot
p <- ggplot(peak_data, aes(x = Regulation, y = Gene.Name)) +
  geom_point(aes(size = Peak_Count, color = Mean_Peak_Score)) +
  geom_hline(yintercept = separator_positions, color = "grey80", size = 0.5) +
  scale_size_continuous(range = c(1, 8), limits = c(0, max(peak_data$Peak_Count)), 
                        breaks = scales::pretty_breaks(n = 5)) +
  scale_color_gradient2(low = "white", mid = "blue", high = "red", 
                        midpoint = median(peak_data$Mean_Peak_Score, na.rm = TRUE),
                        limits = c(0, max(peak_data$Mean_Peak_Score, na.rm = TRUE)),
                        labels = scales::comma) +
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
    title = "ATAC-seq Peak Analysis: Upregulated vs Downregulated Genes"
  ) +
  guides(
    size = guide_legend(override.aes = list(color = "black")),
    color = guide_colorbar(barwidth = 1, barheight = 10)
  )

# Display the plot
print(p)
ggsave("atac_seq_peak_analysis_up_vs_down.png", plot = p, width = 12, height = 20, dpi = 300)



# Calculate fold change and statistical significance
comparison_data <- peak_data %>%
  pivot_wider(names_from = Regulation, values_from = c(Mean_Peak_Score, Peak_Count)) %>%
  mutate(
    Fold_Change = log2(Mean_Peak_Score_Upregulated / Mean_Peak_Score_Downregulated),
    Total_Peaks = Peak_Count_Upregulated + Peak_Count_Downregulated
  ) %>%
  arrange(desc(abs(Fold_Change))) %>%
  mutate(Significant = abs(Fold_Change) > 1)  # You may want to adjust this threshold

# Print the comparison results
print(comparison_data)

# Save the comparison results to a CSV file
write.csv(comparison_data, "gene_comparison_results.csv", row.names = FALSE)

# Create a volcano plot
volcano_plot <- ggplot(comparison_data, aes(x = Fold_Change, y = Total_Peaks)) +
  geom_point(aes(color = Significant)) +
  geom_text(aes(label = ifelse(Significant, as.character(Gene.Name), '')), 
            vjust = 1.5, hjust = 0.5, size = 3) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  theme_minimal() +
  labs(
    x = "Log2 Fold Change (Upregulated / Downregulated)",
    y = "Total Number of Peaks",
    title = "Volcano Plot: Gene Comparison Upregulated vs Downregulated"
  )

# Display the volcano plot
list
# Save the volcano plot
ggsave("volcano_plot_up_vs_down_list1_M3.png", plot = volcano_plot, width = 12, height = 10, dpi = 300)