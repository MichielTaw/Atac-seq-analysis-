contr <- makeContrasts(t24_BMP.WNT_06h_L - t06_BMP.WNT, levels = colnames(design))
#Define FDR threshold and logFC threshold
FDR.thresh <- 0.01
logFC.thresh <- 1



#________________METHOD 1:______________
results_m1 <- glmQLFTest(fit.peak.counts.tmm, contrast = contr)  
rowData(peak.counts.tmm) <- cbind(rowData(fit.peak.counts.tmm), results_m1$table) # combine GRanges rowdata with differential statistics
merged.peaks_m1 <- mergeWindows(rowRanges(peak.counts.tmm), tol=500L, max.width=5000L)
summary(width(merged.peaks_m1$region))
#use most significant window as statistical representation for p-value and FDR for merged windows
tab.best <- getBestTest(merged.peaks_m1$id, results_m1$table)

final.merged.peaks_m1 <- GRanges(cbind(as.data.frame(merged.peaks_m1$region), results_m1$table[tab.best$rep.test, -4], tab.best[,-c(7:8)]))
# sort by FDR
final.merged.peaks_m1 <- final.merged.peaks_m1[order(final.merged.peaks_m1@elementMetadata$FDR), ]
final.merged.peaks_m1

#Filter for significant differentially accessible peaks based on FDR and logFC
significant_peaks_m1 <- final.merged.peaks_m1[
  final.merged.peaks_m1@elementMetadata$FDR < FDR.thresh & 
    abs(final.merged.peaks_m1@elementMetadata$logFC) > logFC.thresh, 
]


# Sort by absolute logFC
significant_peaks_m1 <- significant_peaks_m1[order(-abs(significant_peaks_m1@elementMetadata$logFC)), ]
saveRDS(significant_peaks_m1, file='/collignon/Tawdarous/atac_seq_lab/results/csaw_result/new_csaw/comparaison_paper/t24_BMP.WNT_06h_L_VS_t06_BMP.WNT/significant_peaks_m1.rds')



library(ggplot2)

final.merged.peaks_m1$sig <- "n.s."
final.merged.peaks_m1$sig[
  final.merged.peaks_m1$FDR < FDR.thresh & 
    abs(final.merged.peaks_m1$logFC) > logFC.thresh
] <- "significant"

ggplot(data=data.frame(final.merged.peaks_m1),
       aes(x = logCPM, y = logFC, col = factor(sig, levels=c("n.s.", "significant")))) + 
  geom_point() + scale_color_manual(values = c("black", "red")) + 
  geom_smooth(inherit.aes=F, aes(x = logCPM, y = logFC), method = "loess") + # smoothed loess fit; can add span=0.5 to reduce computation load/time
  geom_hline(yintercept = 0) + labs(col = NULL)

df_ma1 <- as.data.frame(final.merged.peaks_m1)

# Create the MA plot with additional smooth lines
ggplot(data = df_ma1, aes(x = logCPM, y = logFC, col = factor(sig, levels = c("n.s.", "significant")))) + 
  geom_point() + 
  scale_color_manual(values = c("black", "red")) + 
  # Smoothed line for significant points only
  geom_smooth(data = subset(df_ma1, sig == "significant"), aes(x = logCPM, y = logFC), 
              method = "gam", col = "plum", se = FALSE) + 
  geom_hline(yintercept = 0) + 
  labs(col = NULL) +
  theme_bw(base_size = 15)


anno_m1 <- detailRanges(significant_peaks_m1,
                        txdb = TxDb.Mmusculus.UCSC.mm39.knownGene,
                        orgdb = org.Mm.eg.db, 
                        promoter = c(1500, 1500), 
                        dist = 10 * 1e3L)
anno_df_m1 <- as.data.frame(anno_m1)

# Extract relevant metadata from upregulated peaks
metadata_df_m1 <- as.data.frame(significant_peaks_m1)

# Combine the annotation with the metadata
combined_df_m1<- cbind(anno_df_m1, metadata_df_m1[, c("logCPM", "FDR", "direction", "logFC")])
combined_df_m1 <- combined_df_m1[order(combined_df_m1$logCPM),]
saveRDS(combined_df_m1, file='/collignon/Tawdarous/atac_seq_lab/results/csaw_result/new_csaw/comparaison_paper/t24_BMP.WNT_06h_L_VS_t06_BMP.WNT/combined_df_m1.rds')


# View the combined data
combined_df_m1 <- combined_df_m1[order(combined_df_m1$logFC),]
head(combined_df_m1)


#################


#________________METHOD 2:______________
results_m2 <- glmQLFTest(fit.peak.counts.loess, contrast = contr)
rowData(peak.counts.loess) <- cbind(rowData(peak.counts.loess), results_m2$table) # combine GRanges rowdata with differential statistics
merged.peaks_m2 <- mergeWindows(rowRanges(peak.counts.loess), tol=500L, max.width=5000L)
summary(width(merged.peaks_m2$region))
#use most significant window as statistical representation for p-value and FDR for merged windows
tab.best <- getBestTest(merged.peaks_m2$id, results_m2$table)

final.merged.peaks_m2 <- GRanges(cbind(as.data.frame(merged.peaks_m2$region), results_m2$table[tab.best$rep.test, -4], tab.best[,-c(7:8)]))
# sort by FDR
final.merged.peaks_m2 <- final.merged.peaks_m2[order(final.merged.peaks_m2@elementMetadata$FDR), ]
final.merged.peaks_m2

#Filter for significant differentially accessible peaks based on FDR and logFC
significant_peaks_m2 <- final.merged.peaks_m2[
  final.merged.peaks_m2@elementMetadata$FDR < FDR.thresh & 
    abs(final.merged.peaks_m2@elementMetadata$logFC) > logFC.thresh, 
]


# Sort by absolute logFC
significant_peaks_m2 <- significant_peaks_m2[order(-abs(significant_peaks_m2@elementMetadata$logFC)), ]
saveRDS(significant_peaks_m2, file='/collignon/Tawdarous/atac_seq_lab/results/csaw_result/new_csaw/comparaison_paper/t24_BMP.WNT_06h_L_VS_t06_BMP.WNT/significant_peaks_m2.rds')


final.merged.peaks_m2$sig <- "n.s."
final.merged.peaks_m2$sig[
  final.merged.peaks_m2$FDR < FDR.thresh & 
    abs(final.merged.peaks_m2$logFC) > logFC.thresh
] <- "significant"

ggplot(data=data.frame(final.merged.peaks_m2),
       aes(x = logCPM, y = logFC, col = factor(sig, levels=c("n.s.", "significant")))) + 
  geom_point() + scale_color_manual(values = c("black", "red")) + 
  geom_smooth(inherit.aes=F, aes(x = logCPM, y = logFC), method = "loess") + # smoothed loess fit; can add span=0.5 to reduce computation load/time
  geom_hline(yintercept = 0) + labs(col = NULL)

# Convert to data frame for ggplot
df_ma2 <- as.data.frame(final.merged.peaks_m2)

# Create the MA plot with additional smooth lines
ggplot(data = df_ma2, aes(x = logCPM, y = logFC, col = factor(sig, levels = c("n.s.", "significant")))) + 
  geom_point() + 
  scale_color_manual(values = c("black", "red")) + 
  # Smoothed line for significant points only
  geom_smooth(data = subset(df_ma2, sig == "significant"), aes(x = logCPM, y = logFC), 
              method = "loess", col = "plum", se = FALSE) + 
  geom_hline(yintercept = 0) + 
  labs(col = NULL) +
  theme_bw(base_size = 15)


library(MASS)
ggplot(data = df_ma2, aes(x = logCPM, y = logFC, col = factor(sig, levels = c("n.s.", "significant")))) + 
  geom_point() + 
  scale_color_manual(values = c("black", "red")) + 
  geom_smooth(data = subset(df_ma1, sig == "significant"), aes(x = logCPM, y = logFC), 
              method = "rlm", col = "plum", se = FALSE) + 
  geom_hline(yintercept = 0) + 
  labs(col = NULL) +
  theme_bw(base_size = 15)

ggplot(data = df_ma2, aes(x = log(logCPM + 0.1), y = logFC, col = factor(sig, levels = c("n.s.", "significant")))) + 
  geom_point() + 
  scale_color_manual(values = c("black", "red")) + 
  geom_smooth(data = subset(df_ma1, sig == "significant"), aes(x = log(logCPM + 1), y = logFC), 
              method = "loess", col = "plum", se = FALSE) + 
  geom_hline(yintercept = 0) + 
  labs(col = NULL) +
  theme_bw(base_size = 15)

df_ma2$logCPM_clipped <- pmin(df_ma2$logCPM, 1)
ggplot(data = df_ma2, aes(x = logCPM_clipped, y = logFC, col = factor(sig, levels = c("n.s.", "significant")))) + 
  geom_point() + 
  scale_color_manual(values = c("black", "red")) + 
  geom_smooth(data = subset(df_ma2, sig == "significant"), aes(x = logCPM_clipped, y = logFC), 
              method = "loess", col = "plum", se = FALSE) + 
  geom_hline(yintercept = 0) + 
  labs(col = NULL) +
  theme_bw(base_size = 15)

anno_m2 <- detailRanges(significant_peaks_m2,
                        txdb = TxDb.Mmusculus.UCSC.mm39.knownGene,
                        orgdb = org.Mm.eg.db, 
                        promoter = c(3000, 3000), 
                        dist = 10 * 1e3L)
anno_df_m2 <- as.data.frame(anno_m2)

# Extract relevant metadata from upregulated peaks
metadata_df_m2 <- as.data.frame(significant_peaks_m2)

# Combine the annotation with the metadata
combined_df_m2<- cbind(anno_df_m2, metadata_df_m2[, c("logCPM", "FDR", "direction", "logFC")])
combined_df_m2 <- combined_df_m2[order(combined_df_m2$logCPM),]
saveRDS(combined_df_m2, file='/collignon/Tawdarous/atac_seq_lab/results/csaw_result/new_csaw/comparaison_paper/t24_BMP.WNT_06h_L_VS_t06_BMP.WNT/combined_df_m2.rds')


# View the combined data
combined_df_m2 <- combined_df_m2[order(combined_df_m2$logFC),]
head(combined_df_m2)




#________________METHOD 3:______________
results_m3 <- glmQLFTest(fit.counts.local.tmm, contrast = contr)

rowData(counts.local.tmm) <- cbind(rowData(counts.local.tmm), results_m3$table) # combine GRanges rowdata with differential statistics
merged.peaks_m3 <- mergeWindows(rowRanges(counts.local.tmm), tol=500L, max.width=5000L)
summary(width(merged.peaks_m3$region))
#use most significant window as statistical representation for p-value and FDR for merged windows
tab.best <- getBestTest(merged.peaks_m3$id, results_m3$table)

final.merged.peaks_m3 <- GRanges(cbind(as.data.frame(merged.peaks_m3$region), results_m3$table[tab.best$rep.test, -4], tab.best[,-c(7:8)]))
# sort by FDR
final.merged.peaks_m3 <- final.merged.peaks_m3[order(final.merged.peaks_m3@elementMetadata$FDR), ]
final.merged.peaks_m3

#Filter for significant differentially accessible peaks based on FDR and logFC
significant_peaks_m3 <- final.merged.peaks_m3[
  final.merged.peaks_m3@elementMetadata$FDR < FDR.thresh & 
    abs(final.merged.peaks_m3@elementMetadata$logFC) > logFC.thresh, 
]


# Sort by absolute logFC
significant_peaks_m3 <- significant_peaks_m3[order(-abs(significant_peaks_m3@elementMetadata$logFC)), ]
saveRDS(significant_peaks_m3, file='/collignon/Tawdarous/atac_seq_lab/results/csaw_result/new_csaw/comparaison_paper/t24_BMP.WNT_06h_L_VS_t06_BMP.WNT/significant_peaks_m3.rds')

library(ggplot2)
final.merged.peaks_m3$sig <- "n.s."
final.merged.peaks_m3$sig[
  final.merged.peaks_m3$FDR < FDR.thresh & 
    abs(final.merged.peaks_m3$logFC) > logFC.thresh
] <- "significant"
ggplot(data=data.frame(final.merged.peaks_m3),
       aes(x = logCPM, y = logFC, col = factor(sig, levels=c("n.s.", "significant")))) + 
  geom_point() + scale_color_manual(values = c("black", "red")) + 
  geom_smooth(inherit.aes=F, aes(x = logCPM, y = logFC), method = "loess") + # smoothed loess fit; can add span=0.5 to reduce computation load/time
  geom_hline(yintercept = 0) + labs(col = NULL)

# Convert to data frame for ggplot
df_ma3 <- as.data.frame(final.merged.peaks_m3)

# Create the MA plot with additional smooth lines
ggplot(data = df_ma3, aes(x = logCPM, y = logFC, col = factor(sig, levels = c("n.s.", "significant")))) + 
  geom_point() + 
  scale_color_manual(values = c("black", "red")) + 
  # Smoothed line for significant points only
  geom_smooth(data = subset(df_ma3, sig == "significant"), aes(x = logCPM, y = logFC), 
              method = "loess", col = "plum", se = FALSE) + 
  geom_hline(yintercept = 0) + 
  labs(col = NULL) +
  theme_bw(base_size = 15)




anno_m3 <- detailRanges(significant_peaks_m3,
                        txdb = TxDb.Mmusculus.UCSC.mm39.knownGene,
                        orgdb = org.Mm.eg.db, 
                        promoter = c(3000, 3000), 
                        dist = 10 * 1e3L)
anno_df_m3 <- as.data.frame(anno_m3)

# Extract relevant metadata from upregulated peaks
metadata_df_m3 <- as.data.frame(significant_peaks_m3)

# Combine the annotation with the metadata
combined_df_m3<- cbind(anno_df_m3, metadata_df_m3[, c("logCPM", "FDR", "direction", "logFC")])
combined_df_m3 <- combined_df_m3[order(combined_df_m2$logCPM),]
saveRDS(combined_df_m3, file='/collignon/Tawdarous/atac_seq_lab/results/csaw_result/new_csaw/comparaison_paper/t24_BMP.WNT_06h_L_VS_t06_BMP.WNT/combined_df_m3.rds')


# View the combined data
combined_df_m3 <- combined_df_m3[order(combined_df_m3$logFC),]
head(combined_df_m3)



#________________METHOD 4 :______________

results_m4 <- glmQLFTest(fit.counts.local.loess, contrast = contr)

rowData(counts.local.loess) <- cbind(rowData(counts.local.loess), results_m4$table) # combine GRanges rowdata with differential statistics
merged.peaks_m4 <- mergeWindows(rowRanges(counts.local.loess), tol=500L, max.width=5000L)
summary(width(merged.peaks_m4$region))
#use most significant window as statistical representation for p-value and FDR for merged windows
tab.best <- getBestTest(merged.peaks_m4$id, results_m4$table)

final.merged.peaks_m4 <- GRanges(cbind(as.data.frame(merged.peaks_m4$region), results_m4$table[tab.best$rep.test, -4], tab.best[,-c(7:8)]))
# sort by FDR
final.merged.peaks_m4 <- final.merged.peaks_m4[order(final.merged.peaks_m4@elementMetadata$FDR), ]
final.merged.peaks_m4

#Filter for significant differentially accessible peaks based on FDR and logFC
significant_peaks_m4 <- final.merged.peaks_m4[
  final.merged.peaks_m4@elementMetadata$FDR < FDR.thresh & 
    abs(final.merged.peaks_m4@elementMetadata$logFC) > logFC.thresh, 
]


# Sort by absolute logFC
significant_peaks_m4 <- significant_peaks_m4[order(-abs(significant_peaks_m4@elementMetadata$logFC)), ]
saveRDS(significant_peaks_m4, file='/collignon/Tawdarous/atac_seq_lab/results/csaw_result/new_csaw/comparaison_paper/t24_BMP.WNT_06h_L_VS_t06_BMP.WNT/significant_peaks_m4.rds')


final.merged.peaks_m4$sig <- "n.s."
final.merged.peaks_m4$sig[
  final.merged.peaks_m4$FDR < FDR.thresh & 
    abs(final.merged.peaks_m4$logFC) > logFC.thresh
] <- "significant"

# Convert to data frame for ggplot
df_ma4 <- as.data.frame(final.merged.peaks_m4)

# Create the MA plot with additional smooth lines
ggplot(data = df_ma4, aes(x = logCPM, y = logFC, col = factor(sig, levels = c("n.s.", "significant")))) + 
  geom_point() + 
  scale_color_manual(values = c("black", "red")) + 
  # Smoothed line for significant points only
  geom_smooth(data = subset(df_ma4, sig == "significant"), aes(x = logCPM, y = logFC), 
              method = "loess", col = "plum", se = FALSE) + 
  geom_hline(yintercept = 0) + 
  labs(col = NULL) +
  theme_bw(base_size = 15)


anno_m4 <- detailRanges(significant_peaks_m4,
                        txdb = TxDb.Mmusculus.UCSC.mm39.knownGene,
                        orgdb = org.Mm.eg.db, 
                        promoter = c(3000, 3000), 
                        dist = 10 * 1e3L)
anno_df_m4 <- as.data.frame(anno_m4)

# Extract relevant metadata from upregulated peaks
metadata_df_m4 <- as.data.frame(significant_peaks_m4)

# Combine the annotation with the metadata
combined_df_m4<- cbind(anno_df_m4, metadata_df_m4[, c("logCPM", "FDR", "direction", "logFC")])
combined_df_m4 <- combined_df_m4[order(combined_df_m4$logCPM),]
saveRDS(combined_df_m4, file='/collignon/Tawdarous/atac_seq_lab/results/csaw_result/new_csaw/comparaison_paper/t24_BMP.WNT_06h_L_VS_t06_BMP.WNT/combined_df_m4.rds')


# View the combined data
combined_df_m4 <- combined_df_m4[order(combined_df_m4$logFC),]
head(combined_df_m4)


#------------IMPORT-------
setwd("/collignon/Tawdarous/atac_seq_lab/results/csaw_result/new_csaw/comparaison_paper/t24_BMP.WNT_06h_L_VS_t06_BMP.WNT")


significant_peaks_m1 <- readRDS("/collignon/Tawdarous/atac_seq_lab/results/csaw_result/new_csaw/comparaison_paper/t24_BMP.WNT_06h_L_VS_t06_BMP.WNT/significant_peaks_m1.rds")
significant_peaks_m1$name <- paste0("peak_", 1:length(significant_peaks_m1))


significant_peaks_m2 <- readRDS("/collignon/Tawdarous/atac_seq_lab/results/csaw_result/new_csaw/comparaison_paper/t24_BMP.WNT_06h_L_VS_t06_BMP.WNT/significant_peaks_m2.rds")
significant_peaks_m2$name <- paste0("peak_", 1:length(significant_peaks_m2))


significant_peaks_m3 <- readRDS("/collignon/Tawdarous/atac_seq_lab/results/csaw_result/new_csaw/comparaison_paper/t24_BMP.WNT_06h_L_VS_t06_BMP.WNT/significant_peaks_m3.rds")
significant_peaks_m3$name <- paste0("peak_", 1:length(significant_peaks_m3))


significant_peaks_m4 <- readRDS("/collignon/Tawdarous/atac_seq_lab/results/csaw_result/new_csaw/comparaison_paper/t24_BMP.WNT_06h_L_VS_t06_BMP.WNT/significant_peaks_m4.rds")
significant_peaks_m4$name <- paste0("peak_", 1:length(significant_peaks_m4))

#------------
up_m1 <- significant_peaks_m1[significant_peaks_m1$direction=="up"]
down_m1 <- significant_peaks_m1[significant_peaks_m1$direction=="down"]

up_m2 <- significant_peaks_m2[significant_peaks_m2$direction=="up"]
down_m2 <- significant_peaks_m2[significant_peaks_m2$direction=="down"]

up_m3 <- significant_peaks_m3[significant_peaks_m3$direction=="up"]
down_m3 <- significant_peaks_m3[significant_peaks_m3$direction=="down"]

up_m4 <- significant_peaks_m4[significant_peaks_m4$direction=="up"]
down_m4 <- significant_peaks_m4[significant_peaks_m4$direction=="down"]

#___________
create_homer_df <- function(significant_peaks) {
  
  # Create homer_df
  homer_df <- data.frame(
    chr = seqnames(significant_peaks),
    start = start(significant_peaks) - 1,  # BED is 0-based, GRanges is 1-based
    end = end(significant_peaks),
    name = significant_peaks$name,
    score = significant_peaks$logFC,  # Using logFC as score, can be modified
    strand = strand(significant_peaks)
  )
  
  # Modify the strand column: replace "*" with "."
  homer_df$strand <- ifelse(homer_df$strand == "*", ".", as.character(homer_df$strand))
  
  # Add additional columns (logCPM, F_statistic, PValue, FDR)
  homer_df_extended <- homer_df %>%
    mutate(
      logCPM = significant_peaks$logCPM,
      F_statistic = significant_peaks$F,
      PValue = significant_peaks$PValue,
      FDR = significant_peaks$FDR
    )
  
  # Return the extended homer_df
  return(homer_df_extended)
}



homer_up_m1 <- create_homer_df(up_m1)
homer_down_m1 <- create_homer_df(down_m1)

homer_up_m2 <- create_homer_df(up_m2)
homer_down_m2 <- create_homer_df(down_m2)


homer_up_m3 <- create_homer_df(up_m3)
homer_down_m3 <- create_homer_df(down_m3)


homer_up_m4 <- create_homer_df(up_m4)
homer_down_m4 <- create_homer_df(down_m4)


setwd("/collignon/Tawdarous/atac_seq_lab/results/csaw_result/new_csaw/comparaison_paper/t24_BMP.WNT_06h_L_VS_t06_BMP.WNT/homer")

write.table(homer_up_m1, "homer_up_m1.txt", 
            sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(homer_down_m1, "homer_down_m1.txt", 
            sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

write.table(homer_up_m2, "homer_up_m2.txt", 
            sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(homer_down_m2, "homer_down_m2.txt", 
            sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)


write.table(homer_up_m3, "homer_up_m3.txt", 
            sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(homer_down_m3, "homer_down_m3.txt", 
            sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

write.table(homer_up_m4, "homer_up_m4.txt", 
            sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(homer_down_m4, "homer_down_m4.txt", 
            sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)



# Write the HOMER-specific peak file


