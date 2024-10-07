# Load required libraries
library(dplyr)
library(ggplot2)
gene_groups <- list(
  epiblast = c("SLC7A3", "UTF1", "PLA2G10", "PIM2", "L1TD1", "GRP", "DNMT3A", "TNFSF11", "DNMT3B", "PDZD4"),
  PS = c("NANOG", "CACNA1A", "POU5F1", "EOMES", "ISM2", "ZFP184", "B3GAT1", "PRPS1", "DRD2", "GABRG2"),
  `caudal epiblast` = c("TCEAL8", "ETNK2", "GM9885", "FOXB1", "OTOP3", "CXX1B", "NFE2L3", "USH1G", "PSORS1C2", "EPHA1"),
  pgc = c("FGF8", "IFITM3", "PSMB8", "SPRR2A3", "ID1", "DND1", "PSME2", "IGFBP3", "MYCL", "PSME1"),
  `anterior PS` = c("TDGF1", "LHX1", "GSC", "FGF5", "HHEX", "GRRP1", "CDH2", "UPP1", "OTX2", "NPM3"),
  NOTOCHORD = c("FOXJ1", "PIFO", "CFAP126", "1110017D15RIK", "TEKT1", "MLF1", "2410004P03RIK", "RIIAD1", "ARL4A", "CCPG1OS"),
  `DEFINITIVE ENDODERM` = c("ZFP974", "AC121871.1", "2010300C02RIK", "CADPS2", "TRAF3IP3", "GM10767", "GM28172", "GM26794", "M5C1000I18RIK", "INPPL1"),
  `NASCENT MESODERM` = c("LEFTY2", "MESP1", "MESP2", "CBLN1", "PMAIP1", "PRSS22", "TMEM229A", "ARL4D", "GM2694", "RIMBP2"),
  `surface ectoderm` = c("TFAP2A", "GAP43", "BCAM", "RBBP7", "DLX5", "GATA3", "MED7", "PLATR22", "OLIG1", "DLX2"),
  `rostral neuroectoderm` = c("NEUROG1", "PCSK2", "PDZPH1", "MIF", "MIP", "IRX3", "LY6H", "CRYBA1", "FOXI2", "SLC16A3"),
  `caudal neuroectoderm` = c("PTN", "HES3", "USP3", "AKR1B3", "DLGAP3", "TMEM74B", "GBX2", "CASP3", "MCM3", "PLPP4"),
  mesenchyme = c("TDO2", "ACTA2", "TNNT2", "TAGLN", "COLEC11", "AHNAK", "TMEM108", "GDF6", "SNAI2", "DAAM1")
)
# E6.5: list1
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

# E7.5

gene_groups <- list(
  c("SLC7A3", "UTF1", "PLA2G10", "PIM2", "L1TD1", "GRP", "DNMT3A", "TNFSF11", "DNMT3B", "PDZD4"),#epiblast
  c("NANOG", "CACNA1A", "POU5F1", "EOMES", "ISM2", "ZFP184", "B3GAT1", "PRPS1", "DRD2", "GABRG2"),#PS
  c("TCEAL8", "ETNK2", "GM9885", "FOXB1", "OTOP3", "CXX1B", "NFE2L3", "USH1G", "PSORS1C2", "EPHA1"),#caudal epiblast
  c("FGF8", "IFITM3", "PSMB8", "SPRR2A3", "ID1", "DND1", "PSME2", "IGFBP3", "MYCL", "PSME1"),#pgc
  c("TDGF1", "LHX1", "GSC", "FGF5", "HHEX", "GRRP1", "CDH2", "UPP1", "OTX2", "NPM3"),#anterior PS
  c("FOXJ1", "PIFO", "CFAP126", "1110017D15RIK", "TEKT1", "MLF1", "2410004P03RIK", "RIIAD1", "ARL4A", "CCPG1OS"),#NOTOCHORD
  c("ZFP974", "AC121871.1", "2010300C02RIK", "CADPS2", "TRAF3IP3", "GM10767", "GM28172", "GM26794", "M5C1000I18RIK", "INPPL1"),#DEFINITIVE ENDODERM
  c("LEFTY2", "MESP1", "MESP2", "CBLN1", "PMAIP1", "PRSS22", "TMEM229A", "ARL4D", "GM2694", "RIMBP2"),#NASCEND MESODERM
  c("TFAP2A", "GAP43", "BCAM", "RBBP7", "DLX5", "GATA3", "MED7", "PLATR22", "OLIG1", "DLX2"),#surface ectoderm
  c("NEUROG1", "PCSK2", "PDZPH1", "MIF", "MIP", "IRX3", "LY6H", "CRYBA1", "FOXI2", "SLC16A3"),#reustral neuroectoderm
  c("PTN", "HES3", "USP3", "AKR1B3", "DLGAP3", "TMEM74B", "GBX2", "CASP3", "MCM3", "PLPP4"),#caudal neuroectoderm
  c("TDO2", "ACTA2", "TNNT2", "TAGLN", "COLEC11", "AHNAK", "TMEM108", "GDF6", "SNAI2", "DAAM1")#mesenchyme
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

# Flatten the gene_groups list and convert to uppercase for consistency
all_genes <- toupper(unlist(gene_groups))

# Process upregulated data (no absolute value transformation)
up_data <- process_homer_data(homer_up_m1_annotated, all_genes) %>%
  mutate(Regulation = "Upregulated")

# Process downregulated data with absolute value for Peak.Score
down_data <- process_homer_data(homer_down_m1_annotated, all_genes, abs_peak_score = TRUE) %>%
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

# Get the genes that are actually present in our data
present_genes <- unique(peak_data$Gene.Name)

# Calculate positions for separating lines between groups
separator_positions <- cumsum(sapply(gene_groups, function(group) {
  sum(toupper(group) %in% present_genes)
}))
separator_positions <- separator_positions[-length(separator_positions)] + 0.5

# Ensure the gene order matches the original group order
peak_data$Gene.Name <- factor(peak_data$Gene.Name, 
                              levels = rev(all_genes[all_genes %in% present_genes]))

# Create the plot
p <- ggplot(peak_data, aes(x = Regulation, y = Gene.Name)) +
  geom_point(aes(size = Peak_Count, color = Mean_Peak_Score)) +
  
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

# Save the plot
ggsave("list2_M4.png", plot = p, width = 12, height = 20, dpi = 300)



# Calculate fold change and statistical significance
comparison_data <- peak_data %>%
  pivot_wider(names_from = Regulation, values_from = c(Mean_Peak_Score, Peak_Count)) %>%
  mutate(
    Fold_Change = log2(Mean_Peak_Score_Upregulated / Mean_Peak_Score_Downregulated),
    Total_Peaks = Peak_Count_Upregulated + Peak_Count_Downregulated
  ) %>%
  arrange(desc(abs(Fold_Change))) %>%
  mutate(Significant = abs(Fold_Change) > 0.001)  # You may want to adjust this threshold

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
