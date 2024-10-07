library(dplyr)
library(ggplot2)
library(gridExtra)
library(grid)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(grid)

# Define gene groups
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

process_homer_data <- function(data, genes_of_interest, regulation) {
  cols <- colnames(data)
  gene_col <- cols[grep("Gene.Name", cols, ignore.case = TRUE)]
  score_col <- cols[grep("Peak.Score", cols, ignore.case = TRUE)]
  
  if(length(gene_col) == 0 || length(score_col) == 0) {
    stop("Could not find appropriate columns for gene names or peak scores")
  }
  
  processed_data <- data %>%
    select(all_of(c(gene_col, score_col))) %>%
    rename(Gene.Name = !!gene_col, Peak.Score = !!score_col) %>%
    mutate(Gene.Name = toupper(Gene.Name)) %>%
    filter(Gene.Name %in% toupper(genes_of_interest)) %>%
    mutate(Regulation = regulation,
           Abs_Peak_Score = abs(Peak.Score))
  
  return(processed_data)
}

# Flatten the gene_groups list and convert to uppercase for consistency
all_genes_ordered <- toupper(unlist(gene_groups))

gene_to_group <- data.frame(
  Gene.Name = all_genes_ordered,
  Group = rep(names(gene_groups), sapply(gene_groups, length))
)

# Process upregulated and downregulated data
up_data <- process_homer_data(homer_up_m4_annotated, all_genes_ordered, "Upregulated")
down_data <- process_homer_data(homer_down_m4_annotated, all_genes_ordered, "Downregulated")

# Combine all processed data
all_data <- bind_rows(up_data, down_data)

# Summarize peak data
peak_data <- all_data %>%
  group_by(Gene.Name, Regulation) %>%
  summarize(
    Max_Abs_Peak_Score = max(Abs_Peak_Score),
    Peak_Count = n(),
    .groups = 'drop'
  ) %>%
  pivot_wider(
    names_from = Regulation, 
    values_from = c(Max_Abs_Peak_Score, Peak_Count),
    values_fill = list(Max_Abs_Peak_Score = 0, Peak_Count = 0)
  ) %>%
  mutate(
    Score_Difference = Max_Abs_Peak_Score_Upregulated - Max_Abs_Peak_Score_Downregulated,
    Regulation = case_when(
      Score_Difference > 0 ~ "Upregulated",
      Score_Difference < 0 ~ "Downregulated",
      TRUE ~ "No Change"
    )
  ) %>%
  right_join(gene_to_group, by = "Gene.Name") %>%
  arrange(match(Gene.Name, toupper(all_genes_ordered)))

# Ensure all genes are included and in the correct order
peak_data <- peak_data %>%
  mutate(Gene.Name = factor(Gene.Name, levels = rev(toupper(all_genes_ordered))))

# Calculate separator positions
group_sizes <- sapply(gene_groups, length)
separator_positions <- cumsum(group_sizes)[-length(group_sizes)] + 0.5

# Calculate label positions
group_labels <- data.frame(
  Group = names(gene_groups),
  y_position = cumsum(group_sizes) - group_sizes/2
)

# Create main plot
main_plot <- ggplot(peak_data, aes(x = Regulation, y = Gene.Name)) +
  geom_point(aes(size = abs(Score_Difference), color = Score_Difference)) +
  geom_hline(yintercept = separator_positions, color = "grey80", size = 0.5) +
  scale_size_continuous(range = c(1, 8), name = "Absolute Score Difference") +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", 
                        midpoint = 0, name = "Score Difference") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 8),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.title.y = element_blank()
  ) +
  labs(
    x = "", 
    y = "", 
    title = "ATAC-seq Peak Analysis: Upregulated vs Downregulated Genes"
  )

# Create label plot
label_plot <- ggplot(group_labels, aes(x = 1, y = y_position, label = Group)) +
  geom_text(hjust = 1, fontface = "bold", size = 4) +
  scale_y_continuous(limits = c(0, nrow(peak_data))) +
  theme_void() +
  theme(plot.margin = margin(5.5, 0, 5.5, 0))

# Combine plots
combined_plot <- grid.arrange(
  label_plot, main_plot, 
  ncol = 2, 
  widths = c(1, 4),
  respect = TRUE
)

# Display the combined plot
print(combined_plot)

# Save the combined plot
ggsave(
  filename = "E7_5_M4.png",
  plot = combined_plot,
  width = 49,
  height = 15,
  dpi = 400,
  units = "in",
  bg = "white",
  device = "png",
  limitsize = FALSE
)