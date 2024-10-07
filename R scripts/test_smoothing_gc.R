library(qsmooth)
library(SummarizedExperiment)
library(bodymapRat)
library(quantro)



bm_dat <- bodymapRat()
# select brain and liver samples, stage 21 weeks, and only bio reps
keepColumns = (colData(bm_dat)$organ %in% c("Brain", "Liver")) &
  (colData(bm_dat)$stage == 21) & (colData(bm_dat)$techRep == 1)
keepRows = rowMeans(assay(bm_dat)) > 10 # Filter out low counts
bm_dat_e1 <- bm_dat[keepRows,keepColumns]
bm_dat_e1


pd2 <- colData(bm_dat_e1)
counts2 <- assay(bm_dat_e1)[!grepl("^ERCC", 
                                   rownames( assay(bm_dat_e1))), ]
pd2$group <- paste(pd2$organ, pd2$sex, sep="_")



par(mfrow=c(2,2))
#par(mfrow=c(1,1))

pd1 <- colData(processed_peak)
counts1 <- assay(processed_peak)[!grepl("^ERCC", 
                                   rownames( assay(processed_peak))), ]

matboxplot(log2(counts1+1), groupFactor = factor(pd1$CellType),
           main = "Raw data", xaxt="n", 
           ylab = "Expression (log2 scale)")
axis(1, at=seq_len(length(as.character(pd1$CellType))),
     labels=FALSE)
text(seq_len(length(pd1$CellType)), par("usr")[3] -2, 
     labels = pd1$CellType, srt = 90, pos = 1, xpd = TRUE)

matdensity(log2(counts1+1), groupFactor = pd1$CellType, 
           main = "Raw data", ylab= "density",
           xlab = "Expression (log2 scale)")
legend('topright', levels(factor(pd1$CellType)), 
       col = 1:2, lty = 1)

qs_norm_e1 <- qsmooth(object = counts1, group_factor = pd1$CellType)
qs_norm_e1 
