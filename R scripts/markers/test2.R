library(dplyr)
library(ggplot2)

# E6.5: list1
gene_groups <- list(
  epiblast = c("SLC7A3", "STXBP6", "SLC35D3", "SNCA", "UCHL1", "ENPP3", "SYT11", "DNMT3A", "LYPD6B","GM42637"),
  PS = c("CDKN1A", "PDLIM4", "EPHA1", "AXIN2", "LGFBP3", "APELA", "NEFM", "IRF1", "WNT5B", "SALL3"),
  `nascent mesoderm` = c("DDX5", "UTP4", "EIF1A", "STX3", "MTX2", "TMX1", "EXOC5", "GNG10", "TNPO3", "MESP1"),
  `rostral neuroectoderm` = c("XKR4", "GM1992", "GM37381", "RP1", "SOX17", "GM37323", "MRPL15", "LYPLA1", "GM37988", "TCEA1"),
  `visceral endoderm` = c("CD59A", "GPC4", "AMOT", "SFRP1", "SMIM14", "CER1", "LBH", "GGCT", "FAM89A", "KRT8"),
  `EXE endoderm` = c("CYBA", "RBP4", "NPL", "APOA4", "H19", "TTR", "APOA1", "APOC2", "LGALS2", "APOE"),
  `exe ectoderm` = c("GJB3", "CLDN4", "WNT6", "BMP4", "GM9", "AKAP17B", "PRSS50", "DUOXA2", "INPP1", "AA467197"),
  `parietal endoderm` = c("SERPINH1", "CRYAB", "PTH1R", "FST", "TAX1BP3", "F3", "HS3ST1", "SPARC", "SAT1", "SRGN")
)


# E7.5



# Function to process homer annotated data with an option to take absolute value
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
all_genes_ordered <- unlist(gene_groups)

gene_to_group <- data.frame(
  Gene.Name = all_genes_ordered,
  Group = rep(names(gene_groups), sapply(gene_groups, length))
) %>%
  mutate(Gene.Name = toupper(Gene.Name))



# Process upregulated and downregulated data
up_data <- process_homer_data(homer_up_m3_annotated, all_genes, "Upregulated")
down_data <- process_homer_data(homer_down_m3_annotated, all_genes, "Downregulated")

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


group_labels <- peak_data %>%
  group_by(Group) %>%
  summarize(
    y_position = mean(which(levels(peak_data$Gene.Name) %in% Gene.Name)),
    .groups = 'drop'
  )



group_sizes <- sapply(gene_groups, length)
separator_positions <- cumsum(group_sizes)[-length(group_sizes)] + 0.5


peak_data$Gene.Name <- factor(peak_data$Gene.Name, levels = rev(peak_data$Gene.Name))



# Set gene order for plotting
peak_data$Gene.Name <- factor(peak_data$Gene.Name, levels = rev(peak_data$Gene.Name))

# Create the plot
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

# Créer le graphique des étiquettes
label_plot <- ggplot(group_labels, aes(x = 1, y = y_position, label = Group)) +
  geom_text(hjust = 1, fontface = "bold") +
  scale_y_continuous(limits = c(0, nrow(peak_data))) +
  theme_void() +
  theme(plot.margin = margin(5.5, 0, 5.5, 0))

# Combiner les graphiques
combined_plot <- grid.arrange(label_plot, main_plot, ncol = 2, widths = c(1, 4))

# Afficher le graphique combiné
print(combined_plot)

# Sauvegarder le graphique combiné
ggsave("ATAC_seq_peak_analysis_with_external_labels.png", combined_plot, width = 18, height = 20, dpi = 300)




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

# Créer le graphique des étiquettes
label_plot <- ggplot(group_labels, aes(x = 1, y = y_position, label = Group)) +
  geom_text(hjust = 1, fontface = "bold") +
  scale_y_continuous(limits = c(0, nrow(peak_data))) +
  theme_void() +
  theme(plot.margin = margin(5.5, 0, 5.5, 0))

# Combiner les graphiques
combined_plot <- grid.arrange(label_plot, main_plot, ncol = 2, widths = c(1, 4))

# Afficher le graphique combiné
print(combined_plot)

# Sauvegarder le graphique combiné
ggsave("ATAC_seq_peak_analysis_with_external_labels.png", combined_plot, width = 18, height = 20, dpi = 300)