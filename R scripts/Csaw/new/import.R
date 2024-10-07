library(dplyr)

#------------------------IMPORT METHODS-----------------
design <- readRDS("/collignon/Tawdarous/atac_seq_lab/results/csaw_result/new_csaw/paper_csaw/design.rds")


#method1:
y.peak.counts.tmm <- readRDS("/collignon/Tawdarous/atac_seq_lab/results/csaw_result/new_csaw/paper_csaw/y.peak.counts.tmm.rds")
fit.peak.counts.tmm <- readRDS("/collignon/Tawdarous/atac_seq_lab/results/csaw_result/new_csaw/paper_csaw/fit.peak.counts.tmm.rds")



#method2:
y.peak.counts.loess <- readRDS("/collignon/Tawdarous/atac_seq_lab/results/csaw_result/new_csaw/paper_csaw/y.peak.counts.loess.rds")
fit.peak.counts.loess <- readRDS("/collignon/Tawdarous/atac_seq_lab/results/csaw_result/new_csaw/paper_csaw/fit.peak.counts.loess.rds")


#method3:
y.counts.local.tmm <- readRDS("/collignon/Tawdarous/atac_seq_lab/results/csaw_result/new_csaw/paper_csaw/y.counts.local.tmm.rds")
fit.counts.local.tmm <- readRDS("/collignon/Tawdarous/atac_seq_lab/results/csaw_result/new_csaw/paper_csaw/fit.counts.local.tmm.rds")


#method 4:

y.counts.local.loess <- readRDS("/collignon/Tawdarous/atac_seq_lab/results/csaw_result/new_csaw/paper_csaw/y.counts.local.loess.rds")
fit.counts.local.loess <- readRDS("/collignon/Tawdarous/atac_seq_lab/results/csaw_result/new_csaw/paper_csaw/fit.counts.local.loess.rds")

