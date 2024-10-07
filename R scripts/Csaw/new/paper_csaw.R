#-------------LIBS
library(GenomicRanges)
library(csaw)
library(edgeR)



#------------------------PATHS--------
bam_files_path="/collignon/Tawdarous/atac_seq_lab/results/peak_calling/macs3/peaks_annotate"
replicat_peaks_path="/collignon/Tawdarous/atac_seq_lab/results/peak_calling/macs3/replicated_peaks/chrom_mapped"

#------------------------PATHS--------
blacklist <-read.table('/collignon/Tawdarous/atac_seq_lab/scripts/files/mm10-blacklist.v2.Liftover.mm39.bed.txt', sep="\t")[,1:3]
colnames(blacklist) <- c("chrom", "start", "end")
blacklist <- GRanges(blacklist)

#______BROADPEAK______
peak_annotated_path="/collignon/Tawdarous/atac_seq_lab/results/peak_calling/macs3/peaks_annotate"
#_____________________________
t06_BMP_WNT_rep1_path <- paste0( peak_annotated_path ,"/t06_BMP+WNT_rep1_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted_chrom_mapped.broadPeak")
t06_BMP_WNT_rep1 <- read.table(t06_BMP_WNT_rep1_path, sep="\t")[,1:3]
colnames(t06_BMP_WNT_rep1) <- c("chrom", "start", "end")
t06_BMP_WNT_rep1 <- t06_BMP_WNT_rep1[!is.na(t06_BMP_WNT_rep1[, 1]) & t06_BMP_WNT_rep1[, 1] != "", ]
t06_BMP_WNT_rep1 <- GRanges(t06_BMP_WNT_rep1)
#_____________________________

t06_BMP_WNT_rep2_path <- paste0( peak_annotated_path ,"/t06_BMP+WNT_rep2_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted_chrom_mapped.broadPeak")
t06_BMP_WNT_rep2 <- read.table(t06_BMP_WNT_rep2_path, sep="\t")[,1:3]
colnames(t06_BMP_WNT_rep2) <- c("chrom", "start", "end")
t06_BMP_WNT_rep2 <- t06_BMP_WNT_rep2[!is.na(t06_BMP_WNT_rep2[, 1]) & t06_BMP_WNT_rep2[, 1] != "", ]
t06_BMP_WNT_rep2 <- GRanges(t06_BMP_WNT_rep2)

#_____________________________

t06_BMP_WNT_rep3_path <- paste0( peak_annotated_path ,"/t06_BMP+WNT_rep3_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted_chrom_mapped.broadPeak")
t06_BMP_WNT_rep3 <- read.table(t06_BMP_WNT_rep3_path, sep="\t")[,1:3]
colnames(t06_BMP_WNT_rep3) <- c("chrom", "start", "end")
t06_BMP_WNT_rep3 <- t06_BMP_WNT_rep3[!is.na(t06_BMP_WNT_rep3[, 1]) & t06_BMP_WNT_rep3[, 1] != "", ]
t06_BMP_WNT_rep3 <- GRanges(t06_BMP_WNT_rep3)

#_____________________________
t06_WNT_rep1_path <- paste0( peak_annotated_path ,"/t06_WNT_rep1_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted_chrom_mapped.broadPeak")
t06_WNT_rep1 <- read.table(t06_WNT_rep1_path, sep="\t")[,1:3]
colnames(t06_WNT_rep1) <- c("chrom", "start", "end")
t06_WNT_rep1 <- t06_WNT_rep1[!is.na(t06_WNT_rep1[, 1]) & t06_WNT_rep1[, 1] != "", ]
t06_WNT_rep1 <- GRanges(t06_WNT_rep1)

#_____________________________
t06_WNT_rep2_path <- paste0( peak_annotated_path ,"/t06_WNT_rep2_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted_chrom_mapped.broadPeak")
t06_WNT_rep2 <- read.table(t06_WNT_rep2_path, sep="\t")[,1:3]
colnames(t06_WNT_rep2) <- c("chrom", "start", "end")
t06_WNT_rep2 <- t06_WNT_rep2[!is.na(t06_WNT_rep2[, 1]) & t06_WNT_rep2[, 1] != "", ]
t06_WNT_rep2 <- GRanges(t06_WNT_rep2)

#_____________________________
t06_WNT_rep3_path <- paste0( peak_annotated_path ,"/t06_WNT_rep3_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted_chrom_mapped.broadPeak")
t06_WNT_rep3 <- read.table(t06_WNT_rep3_path, sep="\t")[,1:3]
colnames(t06_WNT_rep3) <- c("chrom", "start", "end")
t06_WNT_rep3 <- t06_WNT_rep3[!is.na(t06_WNT_rep3[, 1]) & t06_WNT_rep3[, 1] != "", ]
t06_WNT_rep3 <- GRanges(t06_WNT_rep3)


#_____________________________
t24_BMP_WNT_06h_LI_rep1_path <- paste0( peak_annotated_path ,"/t24_BMP+WNT-06h_LI_rep1_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted_chrom_mapped.broadPeak")
t24_BMP_WNT_06h_LI_rep1 <- read.table(t24_BMP_WNT_06h_LI_rep1_path, sep="\t")[,1:3]
colnames(t24_BMP_WNT_06h_LI_rep1) <- c("chrom", "start", "end")
t24_BMP_WNT_06h_LI_rep1 <- t24_BMP_WNT_06h_LI_rep1[!is.na(t24_BMP_WNT_06h_LI_rep1[, 1]) & t24_BMP_WNT_06h_LI_rep1[, 1] != "", ]
t24_BMP_WNT_06h_LI_rep1 <- GRanges(t24_BMP_WNT_06h_LI_rep1)

#_____________________________
t24_BMP_WNT_06h_LI_rep2_path <- paste0( peak_annotated_path ,"/t24_BMP+WNT-06h_LI_rep2_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted_chrom_mapped.broadPeak")
t24_BMP_WNT_06h_LI_rep2 <- read.table(t24_BMP_WNT_06h_LI_rep2_path, sep="\t")[,1:3]
colnames(t24_BMP_WNT_06h_LI_rep2) <- c("chrom", "start", "end")
t24_BMP_WNT_06h_LI_rep2 <- t24_BMP_WNT_06h_LI_rep2[!is.na(t24_BMP_WNT_06h_LI_rep2[, 1]) & t24_BMP_WNT_06h_LI_rep2[, 1] != "", ]
t24_BMP_WNT_06h_LI_rep2 <- GRanges(t24_BMP_WNT_06h_LI_rep2)

#_____________________________
t24_BMP_WNT_06h_LI_rep3_path <- paste0( peak_annotated_path ,"/t24_BMP+WNT-06h_LI_rep3_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted_chrom_mapped.broadPeak")
t24_BMP_WNT_06h_LI_rep3 <- read.table(t24_BMP_WNT_06h_LI_rep3_path, sep="\t")[,1:3]
colnames(t24_BMP_WNT_06h_LI_rep3) <- c("chrom", "start", "end")
t24_BMP_WNT_06h_LI_rep3 <- t24_BMP_WNT_06h_LI_rep3[!is.na(t24_BMP_WNT_06h_LI_rep3[, 1]) & t24_BMP_WNT_06h_LI_rep3[, 1] != "", ]
t24_BMP_WNT_06h_LI_rep3 <- GRanges(t24_BMP_WNT_06h_LI_rep3)


#_____________________________
t24_BMP_WNT_06h_L_rep1_path <- paste0( peak_annotated_path ,"/t24_BMP+WNT-06h_L_rep1_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted_chrom_mapped.broadPeak")
t24_BMP_WNT_06h_L_rep1 <- read.table(t24_BMP_WNT_06h_L_rep1_path, sep="\t")[,1:3]
colnames(t24_BMP_WNT_06h_L_rep1) <- c("chrom", "start", "end")
t24_BMP_WNT_06h_L_rep1 <- t24_BMP_WNT_06h_L_rep1[!is.na(t24_BMP_WNT_06h_L_rep1[, 1]) & t24_BMP_WNT_06h_L_rep1[, 1] != "", ]
t24_BMP_WNT_06h_L_rep1 <- GRanges(t24_BMP_WNT_06h_L_rep1)

#_____________________________
t24_BMP_WNT_06h_L_rep2_path <- paste0( peak_annotated_path ,"/t24_BMP+WNT-06h_L_rep2_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted_chrom_mapped.broadPeak")
t24_BMP_WNT_06h_L_rep2 <- read.table(t24_BMP_WNT_06h_L_rep2_path, sep="\t")[,1:3]
colnames(t24_BMP_WNT_06h_L_rep2) <- c("chrom", "start", "end")
t24_BMP_WNT_06h_L_rep2 <- t24_BMP_WNT_06h_L_rep2[!is.na(t24_BMP_WNT_06h_L_rep2[, 1]) & t24_BMP_WNT_06h_L_rep2[, 1] != "", ]
t24_BMP_WNT_06h_L_rep2 <- GRanges(t24_BMP_WNT_06h_L_rep2)

#_____________________________
t24_BMP_WNT_06h_L_rep3_path <- paste0( peak_annotated_path ,"/t24_BMP+WNT-06h_L_rep3_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted_chrom_mapped.broadPeak")
t24_BMP_WNT_06h_L_rep3 <- read.table(t24_BMP_WNT_06h_L_rep3_path, sep="\t")[,1:3]
colnames(t24_BMP_WNT_06h_L_rep3) <- c("chrom", "start", "end")
t24_BMP_WNT_06h_L_rep3 <- t24_BMP_WNT_06h_L_rep3[!is.na(t24_BMP_WNT_06h_L_rep3[, 1]) & t24_BMP_WNT_06h_L_rep3[, 1] != "", ]
t24_BMP_WNT_06h_L_rep3 <- GRanges(t24_BMP_WNT_06h_L_rep3)


#_____________________________
t24_BMP_WNT_06h_LX1_rep1_path <- paste0( peak_annotated_path ,"/t24_BMP+WNT-06h_LX1_rep1_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted_chrom_mapped.broadPeak")
t24_BMP_WNT_06h_LX1_rep1 <- read.table(t24_BMP_WNT_06h_LX1_rep1_path, sep="\t")[,1:3]
colnames(t24_BMP_WNT_06h_LX1_rep1) <- c("chrom", "start", "end")
t24_BMP_WNT_06h_LX1_rep1 <- t24_BMP_WNT_06h_LX1_rep1[!is.na(t24_BMP_WNT_06h_LX1_rep1[, 1]) & t24_BMP_WNT_06h_LX1_rep1[, 1] != "", ]
t24_BMP_WNT_06h_LX1_rep1 <- GRanges(t24_BMP_WNT_06h_LX1_rep1)

#_____________________________
t24_BMP_WNT_06h_LX1_rep2_path <- paste0( peak_annotated_path ,"/t24_BMP+WNT-06h_LX1_rep2_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted_chrom_mapped.broadPeak")
t24_BMP_WNT_06h_LX1_rep2 <- read.table(t24_BMP_WNT_06h_LX1_rep2_path, sep="\t")[,1:3]
colnames(t24_BMP_WNT_06h_LX1_rep2) <- c("chrom", "start", "end")
t24_BMP_WNT_06h_LX1_rep2 <- t24_BMP_WNT_06h_LX1_rep2[!is.na(t24_BMP_WNT_06h_LX1_rep2[, 1]) & t24_BMP_WNT_06h_LX1_rep2[, 1] != "", ]
t24_BMP_WNT_06h_LX1_rep2 <- GRanges(t24_BMP_WNT_06h_LX1_rep2)

#_____________________________
t24_BMP_WNT_06h_LX1_rep3_path <- paste0( peak_annotated_path ,"/t24_BMP+WNT-06h_LX1_rep3_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted_chrom_mapped.broadPeak")
t24_BMP_WNT_06h_LX1_rep3 <- read.table(t24_BMP_WNT_06h_LX1_rep3_path, sep="\t")[,1:3]
colnames(t24_BMP_WNT_06h_LX1_rep3) <- c("chrom", "start", "end")
t24_BMP_WNT_06h_LX1_rep3 <- t24_BMP_WNT_06h_LX1_rep3[!is.na(t24_BMP_WNT_06h_LX1_rep3[, 1]) & t24_BMP_WNT_06h_LX1_rep3[, 1] != "", ]
t24_BMP_WNT_06h_LX1_rep3 <- GRanges(t24_BMP_WNT_06h_LX1_rep3)


#_____________________________
t24_WNT_06h_L_rep1_path <- paste0( peak_annotated_path ,"/t24_WNT-06h_L_rep1_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted_chrom_mapped.broadPeak")
t24_WNT_06h_L_rep1 <- read.table(t24_WNT_06h_L_rep1_path, sep="\t")[,1:3]
colnames(t24_WNT_06h_L_rep1) <- c("chrom", "start", "end")
t24_WNT_06h_L_rep1 <- t24_WNT_06h_L_rep1[!is.na(t24_WNT_06h_L_rep1[, 1]) & t24_WNT_06h_L_rep1[, 1] != "", ]
t24_WNT_06h_L_rep1 <- GRanges(t24_WNT_06h_L_rep1)

#_____________________________
t24_WNT_06h_L_rep2_path <- paste0( peak_annotated_path ,"/t24_WNT-06h_L_rep2_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted_chrom_mapped.broadPeak")
t24_WNT_06h_L_rep2 <- read.table(t24_WNT_06h_L_rep2_path, sep="\t")[,1:3]
colnames(t24_WNT_06h_L_rep2) <- c("chrom", "start", "end")
t24_WNT_06h_L_rep2 <- t24_WNT_06h_L_rep2[!is.na(t24_WNT_06h_L_rep2[, 1]) & t24_WNT_06h_L_rep2[, 1] != "", ]
t24_WNT_06h_L_rep2 <- GRanges(t24_WNT_06h_L_rep2)

#_____________________________
t24_WNT_06h_L_rep3_path <- paste0( peak_annotated_path ,"/t24_WNT-06h_L_rep3_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted_chrom_mapped.broadPeak")
t24_WNT_06h_L_rep3 <- read.table(t24_WNT_06h_L_rep3_path, sep="\t")[,1:3]
colnames(t24_WNT_06h_L_rep3) <- c("chrom", "start", "end")
t24_WNT_06h_L_rep3 <- t24_WNT_06h_L_rep3[!is.na(t24_WNT_06h_L_rep3[, 1]) & t24_WNT_06h_L_rep3[, 1] != "", ]
t24_WNT_06h_L_rep3 <- GRanges(t24_WNT_06h_L_rep3)

#____________________consensus peak set____________________
# Intersection for t06_BMP_WNT_consensus
temp_intersect_1 <- intersect(t06_BMP_WNT_rep1, t06_BMP_WNT_rep2, ignore.strand = TRUE)
t06_BMP_WNT_consensus <- intersect(temp_intersect_1, t06_BMP_WNT_rep3, ignore.strand = TRUE)

# Union for t06_WNT_consensus
temp_union_1 <- intersect(t06_WNT_rep1, t06_WNT_rep2, ignore.strand = TRUE)
t06_WNT_consensus <- intersect(temp_union_1, t06_WNT_rep3, ignore.strand = TRUE)

# Union for t24_BMP_WNT_06h_LI_consensus
temp_union_2 <- intersect(t24_BMP_WNT_06h_LI_rep1, t24_BMP_WNT_06h_LI_rep2, ignore.strand = TRUE)
t24_BMP_WNT_06h_LI_consensus <- intersect(temp_union_2, t24_BMP_WNT_06h_LI_rep3, ignore.strand = TRUE)

# Union for t24_BMP_WNT_06h_L_consensus
temp_union_3 <- intersect(t24_BMP_WNT_06h_L_rep1, t24_BMP_WNT_06h_L_rep2, ignore.strand = TRUE)
t24_BMP_WNT_06h_L_consensus <- intersect(temp_union_3, t24_BMP_WNT_06h_L_rep3, ignore.strand = TRUE)

# Union for t24_BMP_WNT_06h_LX1_consensus
temp_union_4 <- intersect(t24_BMP_WNT_06h_LX1_rep1, t24_BMP_WNT_06h_LX1_rep2, ignore.strand = TRUE)
t24_BMP_WNT_06h_LX1_consensus <- intersect(temp_union_4, t24_BMP_WNT_06h_LX1_rep3, ignore.strand = TRUE)

# Union for t24_WNT_06h_consensus
temp_union_5 <- intersect(t24_WNT_06h_L_rep1, t24_WNT_06h_L_rep2, ignore.strand = TRUE)
t24_WNT_06h_consensus <- intersect(temp_union_5, t24_WNT_06h_L_rep3, ignore.strand = TRUE)

# Sequential union of all consensus peaks
temp_union_1 <- union(t06_BMP_WNT_consensus, t06_WNT_consensus, ignore.strand = TRUE)
temp_union_2 <- union(temp_union_1, t24_BMP_WNT_06h_LI_consensus, ignore.strand = TRUE)
temp_union_3 <- union(temp_union_2, t24_BMP_WNT_06h_L_consensus, ignore.strand = TRUE)
temp_union_4 <- union(temp_union_3, t24_BMP_WNT_06h_LX1_consensus, ignore.strand = TRUE)
all.peaks <- union(temp_union_4, t24_WNT_06h_consensus, ignore.strand = TRUE)



#______BAMFILES____

#______BAMFILES____
dir.bam= "/collignon/Tawdarous/atac_seq_lab/results/post_alignment/shift_chr_renamed"


t06_BMP_WNT.rep1.bam=file.path(dir.bam,"t06_BMP+WNT_rep1_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted_chr_renamed.sorted.bam")
t06_BMP_WNT.rep1.bam.sizes <- getPESizes("t06_BMP_WNT.rep1.bam")
hist(t06_BMP_WNT.rep1.bam.sizes$sizes) 

t06_BMP_WNT.rep2.bam=file.path(dir.bam,"t06_BMP+WNT_rep2_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted_chr_renamed.bam")
t06_BMP_WNT.rep2.bam.sizes <- getPESizes(t06_BMP_WNT.rep2.bam)
hist(t06_BMP_WNT.rep2.bam.sizes$sizes) 

t06_BMP_WNT.rep3.bam=file.path(dir.bam,"t06_BMP+WNT_rep3_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted_chr_renamed.bam")
t06_BMP_WNT.rep3.bam.sizes <- getPESizes(t06_BMP_WNT.rep3.bam)
hist(t06_BMP_WNT.rep3.bam.sizes$sizes) 

#-----------
t06_WNT.rep1.bam=file.path(dir.bam,"t06_WNT_rep1_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted_chr_renamed.bam")
t06_WNT.rep1.bam.sizes <- getPESizes(t06_WNT.rep1.bam)
hist(t06_WNT.rep1.bam.sizes$sizes) 


t06_WNT.rep2.bam=file.path(dir.bam,"t06_WNT_rep2_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted_chr_renamed.bam")
t06_WNT.rep2.bam.sizes <- getPESizes(t06_WNT.rep2.bam)
hist(t06_WNT.rep2.bam.sizes$sizes) 


t06_WNT.rep3.bam=file.path(dir.bam,"t06_WNT_rep3_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted_chr_renamed.bam")
t06_WNT.rep3.bam.sizes <- getPESizes(t06_WNT.rep3.bam)
hist(t06_WNT.rep3.bam.sizes$sizes) 

#-----------
t24_BMP_WNT_06h_LI.rep1.bam=file.path(dir.bam,"t24_BMP+WNT-06h_LI_rep1_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted_chr_renamed.bam")
t24_BMP_WNT_06h_LI.rep1.bam.sizes <- getPESizes(t24_BMP_WNT_06h_LI.rep1.bam)
hist(t24_BMP_WNT_06h_LI.rep1.bam.sizes$sizes)


t24_BMP_WNT_06h_LI.rep2.bam=file.path(dir.bam,"t24_BMP+WNT-06h_LI_rep2_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted_chr_renamed.bam")
t24_BMP_WNT_06h_LI.rep2.bam.sizes <- getPESizes(t24_BMP_WNT_06h_LI.rep2.bam)
hist(t24_BMP_WNT_06h_LI.rep2.bam.sizes$sizes)


t24_BMP_WNT_06h_LI.rep3.bam=file.path(dir.bam,"t24_BMP+WNT-06h_LI_rep3_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted_chr_renamed.bam")
t24_BMP_WNT_06h_LI.rep3.bam.sizes <- getPESizes(t24_BMP_WNT_06h_LI.rep3.bam)
hist(t24_BMP_WNT_06h_LI.rep3.bam.sizes$sizes) 


#--------------
t24_BMP_WNT_06h_L.rep1.bam=file.path(dir.bam,"t24_BMP+WNT-06h_L_rep1_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted_chr_renamed.bam")
t24_BMP_WNT_06h_L.rep1.bam.sizes <- getPESizes(t24_BMP_WNT_06h_L.rep1.bam)
hist(t24_BMP_WNT_06h_L.rep1.bam.sizes$sizes) 

t24_BMP_WNT_06h_L.rep2.bam=file.path(dir.bam,"t24_BMP+WNT-06h_L_rep2_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted_chr_renamed.bam")
t24_BMP_WNT_06h_L.rep2.bam.sizes <- getPESizes(t24_BMP_WNT_06h_L.rep2.bam)
hist(t24_BMP_WNT_06h_L.rep2.bam.sizes$sizes) 

t24_BMP_WNT_06h_L.rep3.bam=file.path(dir.bam,"t24_BMP+WNT-06h_L_rep3_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted_chr_renamed.bam")
t24_BMP_WNT_06h_L.rep3.bam.sizes <- getPESizes(t24_BMP_WNT_06h_L.rep3.bam)
hist(t24_BMP_WNT_06h_L.rep3.bam.sizes$sizes) 


#--------------
t24_BMP_WNT_06h_LX1.rep1.bam=file.path(dir.bam,"t24_BMP+WNT-06h_LX1_rep1_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted_chr_renamed.bam")
t24_BMP_WNT_06h_LX1.rep1.bam.sizes <- getPESizes(t24_BMP_WNT_06h_LX1.rep1.bam)
hist(t24_BMP_WNT_06h_LX1.rep1.bam.sizes$sizes)


t24_BMP_WNT_06h_LX1.rep2.bam=file.path(dir.bam,"t24_BMP+WNT-06h_LX1_rep2_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted_chr_renamed.bam")
t24_BMP_WNT_06h_LX1.rep2.bam.sizes <- getPESizes(t24_BMP_WNT_06h_LX1.rep2.bam)
hist(t24_BMP_WNT_06h_LX1.rep2.bam.sizes$sizes) 


t24_BMP_WNT_06h_LX1.rep3.bam=file.path(dir.bam,"t24_BMP+WNT-06h_LX1_rep3_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted_chr_renamed.bam")
t24_BMP_WNT_06h_LX1.rep3.bam.sizes <- getPESizes(t24_BMP_WNT_06h_LX1.rep3.bam)
hist(t24_BMP_WNT_06h_LX1.rep3.bam.sizes$sizes) 
#-----------------
t24_WNT_06h_L_rep1.bam=file.path(dir.bam,"t24_WNT-06h_L_rep1_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted_chr_renamed.bam")
t24_WNT_06h_L_rep1.bam.sizes <- getPESizes(t24_WNT_06h_L_rep1.bam)
hist(t24_WNT_06h_L_rep1.bam.sizes$sizes)


t24_WNT_06h_L_rep2.bam=file.path(dir.bam,"t24_WNT-06h_L_rep2_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted_chr_renamed.bam")
t24_WNT_06h_L_rep2.bam.sizes <- getPESizes(t24_WNT_06h_L_rep2.bam)
hist(t24_WNT_06h_L_rep2.bam.sizes$sizes) 


t24_WNT_06h_L_rep3.bam=file.path(dir.bam,"t24_WNT-06h_L_rep3_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted_chr_renamed.bam")
t24_WNT_06h_L_rep3.bam.sizes <- getPESizes(t24_WNT_06h_L_rep3.bam)
hist(t24_WNT_06h_L_rep3.bam.sizes$sizes) 


pe.bams<- c(t06_BMP_WNT.rep1.bam,t06_BMP_WNT.rep2.bam,t06_BMP_WNT.rep3.bam,
            t06_WNT.rep1.bam,t06_WNT.rep2.bam,t06_WNT.rep3.bam,
            t24_BMP_WNT_06h_LI.rep1.bam,t24_BMP_WNT_06h_LI.rep2.bam,t24_BMP_WNT_06h_LI.rep3.bam,
            t24_BMP_WNT_06h_L.rep1.bam,t24_BMP_WNT_06h_L.rep2.bam,t24_BMP_WNT_06h_L.rep3.bam,
            t24_BMP_WNT_06h_LX1.rep1.bam,t24_BMP_WNT_06h_LX1.rep2.bam,t24_BMP_WNT_06h_LX1.rep3.bam,
            t24_WNT_06h_L_rep1.bam,t24_WNT_06h_L_rep2.bam,t24_WNT_06h_L_rep3.bam)



#----------------------------PARAMETERS
# define read parameters
standard.chr <- paste0("chr", c(1:19, "X", "Y")) # only use standard chromosomes
param <- readParam(max.frag=1000, pe="both", discard=blacklist, restrict=standard.chr)
#------------------------------------------------COUNT READS IN PEAKS------------
peak.counts <- regionCounts(pe.bams, all.peaks, param=param)
rownames(peak.counts) <- paste0(seqnames(all.peaks), ":", start(all.peaks), "-", end(all.peaks))
colnames(peak.counts) <- basename(colData(peak.counts)$bam.files)
#saveRDS(peak.counts, file='/collignon/Tawdarous/atac_seq_lab/results/csaw_result/new_csaw/paper_csaw/peak_counts.rds')
peak.abundances <- aveLogCPM(asDGEList(peak.counts)) 
peak.counts.filt <- peak.counts[peak.abundances > -3, ] # only use peaks logCPM > -3
#saveRDS(peak.counts.filt, file='/collignon/Tawdarous/atac_seq_lab/results/csaw_result/new_csaw/paper_csaw/peak.counts.filt.rds')
#------------------------------------------------COUNT READS IN WINDOW------------------------
counts <- windowCounts(pe.bams, width=250, param=param) # set width as desired from the fragment length distribution analyses
saveRDS(counts, file='/collignon/Tawdarous/atac_seq_lab/results/csaw_result/new_csaw/paper_csaw/counts_250W.rds')

#__________filter window by local enrichment
# -----------------local background estimator: 2kb neighborhood
neighbor <- suppressWarnings(resize(rowRanges(counts), width=2000, fix="center")) # change width parameter as desired
wider <- regionCounts(pe.bams, regions=neighbor, param=param) # count reads in neighborhoods
saveRDS(wider, file='/collignon/Tawdarous/atac_seq_lab/results/csaw_result/new_csaw/paper_csaw/wider_2K.rds')
# filter.stat <- filterWindows(counts, wider, type="local") # the filterWindows() function is deprecated and has been replaced by filterWindowsLocal(). This is an archived step.
filter.stat <- filterWindowsLocal(counts, wider)
counts.local.filt <- counts[filter.stat$filter > log2(3),] # threshold of 3-fold increase in enrichment over 2kb neighborhood abundance; change as desired
saveRDS(counts.local.filt, file='/collignon/Tawdarous/atac_seq_lab/results/csaw_result/new_csaw/paper_csaw/counts.local.filt_2K.rds')

#-------------------count BAM background bins (for TMM normalization)-------------
binned <- windowCounts(pe.bams, bin=TRUE, width=10000, param=param)
saveRDS(binned, file='/collignon/Tawdarous/atac_seq_lab/results/csaw_result/new_csaw/paper_csaw/binned.rds')


#------------------------------------------NORMALIZATION

# method 1: MACS2 peaks only, TMM normalization based on binned counts
peak.counts.tmm <- peak.counts.filt
peak.counts.tmm <- normFactors(binned, se.out=peak.counts.tmm)
peak.counts.tmm
saveRDS(peak.counts.tmm, file='/collignon/Tawdarous/atac_seq_lab/results/csaw_result/new_csaw/paper_csaw/peak.counts.tmm.rds')

y.peak.counts.tmm <- asDGEList(peak.counts.tmm,
                               genes = as.data.frame(rowData(peak.counts.tmm)),
                               group = colData(peak.counts.tmm)$CellType)
colnames(y.peak.counts.tmm$counts) <- colnames(peak.counts.tmm)
rownames(y.peak.counts.tmm$samples) <- colnames(peak.counts.tmm)
#y.peak.counts.tmm$samples$group
design <- model.matrix(~0+group, data=y.peak.counts.tmm$samples)
colnames(design) <- sub('group','',colnames(design))
design
y.peak.counts.tmm <- estimateDisp(y.peak.counts.tmm, design)
saveRDS(y.peak.counts.tmm, file='/collignon/Tawdarous/atac_seq_lab/results/csaw_result/new_csaw/paper_csaw/y.peak.counts.tmm.rds')
fit.peak.counts.tmm <- glmQLFit(y.peak.counts.tmm, design, robust=TRUE)
saveRDS(fit.peak.counts.tmm, file='/collignon/Tawdarous/atac_seq_lab/results/csaw_result/new_csaw/paper_csaw/fit.peak.counts.tmm.rds')


colnames(design) <- gsub("-", "_", colnames(design))
colnames(design) <- make.names(colnames(design), unique = TRUE)


#__________________________t06_BMP.WNT VS t24_BMP.WNT_06h_L _______________________
contr <- makeContrasts(t06_BMP.WNT - t24_BMP.WNT_06h_L, levels = colnames(design))
results <- glmQLFTest(fit.peak.counts.tmm, contrast = contr)
rowData(peak.counts.tmm) <- cbind(rowData(peak.counts.tmm), results$table) # combine GRanges rowdata with differential statistics

# merge nearby windows
# up to "tol" distance apart: 500 bp in this case; max merged window width: 5000 bp
merged.peaks <- mergeWindows(rowRanges(peak.counts.tmm), tol=500L, max.width=5000L)
# summary(width(merged.peaks$region))
# should merge some peaks; change as desired

#use most significant window as statistical representation for p-value and FDR for merged windows
tab.best <- getBestTest(merged.peaks$id, results$table)
final.merged.peaks <- GRanges(cbind(as.data.frame(merged.peaks$region), results$table[tab.best$rep.test, -4], tab.best[,-c(7:8)]))

# sort by FDR
final.merged.peaks <- final.merged.peaks[order(final.merged.peaks@elementMetadata$FDR), ]
final.merged.peaks


FDR.thresh <- 0.05 # set as desired
final.merged.peaks.sig <- final.merged.peaks[final.merged.peaks@elementMetadata$FDR < FDR.thresh, ]
final.merged.peaks.sig # significant differentially-accessible windows

library(ggplot2)

final.merged.peaks$sig <- "n.s."
final.merged.peaks$sig[final.merged.peaks$FDR < FDR.thresh] <- "significant"

ggplot(data=data.frame(final.merged.peaks),
       aes(x = logCPM, y = logFC, col = factor(sig, levels=c("n.s.", "significant")))) + 
  geom_point() + scale_color_manual(values = c("black", "red")) + 
  geom_smooth(inherit.aes=F, aes(x = logCPM, y = logFC), method = "loess") + # smoothed loess fit; can add span=0.5 to reduce computation load/time
  geom_hline(yintercept = 0) + labs(col = NULL)

# method 2: MACS2 peaks only, csaw loess-normalization
peak.counts.loess <- peak.counts.filt
peak.counts.loess <- normOffsets(peak.counts.loess, se.out=TRUE) # type="loess" is now default
colData(peak.counts.loess)$CellType <- cell_types
peak.counts.loess
saveRDS(peak.counts.loess, file='/collignon/Tawdarous/atac_seq_lab/results/csaw_result/new_csaw/paper_csaw/peak.counts.loess.rds')



y.peak.counts.loess<- asDGEList(peak.counts.loess,
                               genes = as.data.frame(rowData(peak.counts.loess)),
                               group = colData(peak.counts.loess)$CellType)
colnames(y.peak.counts.loess$counts) <- colnames(peak.counts.loess)
rownames(y.peak.counts.loess$samples) <- colnames(peak.counts.loess)
#y.peak.counts.loess$samples$group
design <- model.matrix(~0+group, data=y.peak.counts.loess$samples)
colnames(design) <- sub('group','',colnames(design))
design
y.peak.counts.loess <- estimateDisp(y.peak.counts.loess, design)
saveRDS(y.peak.counts.loess, file='/collignon/Tawdarous/atac_seq_lab/results/csaw_result/new_csaw/paper_csaw/y.peak.counts.loess.rds')
fit.peak.counts.loess <- glmQLFit(y.peak.counts.loess, design, robust=TRUE)
saveRDS(fit.peak.counts.loess, file='/collignon/Tawdarous/atac_seq_lab/results/csaw_result/new_csaw/paper_csaw/fit.peak.counts.loess.rds')


colnames(design) <- gsub("-", "_", colnames(design))
colnames(design) <- make.names(colnames(design), unique = TRUE)


#__________________________t06_BMP.WNT VS t24_BMP.WNT_06h_L _______________________
contr <- makeContrasts(t06_BMP.WNT - t24_BMP.WNT_06h_L, levels = colnames(design))
results <- glmQLFTest(fit.peak.counts.loess, contrast = contr)
rowData(peak.counts.loess) <- cbind(rowData(peak.counts.loess), results$table) # combine GRanges rowdata with differential statistics

# merge nearby windows
# up to "tol" distance apart: 500 bp in this case; max merged window width: 5000 bp
merged.peaks <- mergeWindows(rowRanges(peak.counts.tmm), tol=500L, max.width=5000L)
# summary(width(merged.peaks$region))
# should merge some peaks; change as desired

#use most significant window as statistical representation for p-value and FDR for merged windows
tab.best <- getBestTest(merged.peaks$id, results$table)
final.merged.peaks <- GRanges(cbind(as.data.frame(merged.peaks$region), results$table[tab.best$rep.test, -4], tab.best[,-c(7:8)]))

# sort by FDR
final.merged.peaks <- final.merged.peaks[order(final.merged.peaks@elementMetadata$FDR), ]
final.merged.peaks


FDR.thresh <- 0.05 # set as desired
final.merged.peaks.sig <- final.merged.peaks[final.merged.peaks@elementMetadata$FDR < FDR.thresh, ]
final.merged.peaks.sig # significant differentially-accessible windows

library(ggplot2)

final.merged.peaks$sig <- "n.s."
final.merged.peaks$sig[final.merged.peaks$FDR < FDR.thresh] <- "significant"

ggplot(data=data.frame(final.merged.peaks),
       aes(x = logCPM, y = logFC, col = factor(sig, levels=c("n.s.", "significant")))) + 
  geom_point() + scale_color_manual(values = c("black", "red")) + 
  geom_smooth(inherit.aes=F, aes(x = logCPM, y = logFC), method = "loess") + # smoothed loess fit; can add span=0.5 to reduce computation load/time
  geom_hline(yintercept = 0) + labs(col = NULL)

upregulated <- final.merged.peaks[final.merged.peaks$direction == "up"]
upregulated <- upregulated[order(upregulated$logCPM),]

anno <- detailRanges(upregulated,
                     txdb = TxDb.Mmusculus.UCSC.mm39.knownGene,
                     orgdb = org.Mm.eg.db, 
                     promoter = c(1500, 1500), 
                     dist = 10 * 1e3L)
anno_df <- as.data.frame(anno)

# Extract relevant metadata from upregulated peaks
metadata_df <- as.data.frame(upregulated)

# Combine the annotation with the metadata
combined_df <- cbind(anno_df, metadata_df[, c("logCPM", "FDR", "direction")])

# View the combined data
combined_df <- combined_df[order(combined_df$FDR),]
head(combined_df)








# from vignette: "For type="loess", a numeric matrix of the same dimensions as counts, containing the log-based offsets for use in GLM fitting."

# method 3: csaw de novo peaks by local enrichment, TMM normalization based on binned counts
counts.local.tmm <- counts.local.filt
counts.local.tmm <- normFactors(binned, se.out=counts.local.tmm)
colData(counts.local.tmm)$CellType <- cell_types

saveRDS(counts.local.tmm, file='/collignon/Tawdarous/atac_seq_lab/results/csaw_result/new_csaw/paper_csaw/counts.local.tmm.rds')


y.counts.local.tmm<- asDGEList(counts.local.tmm,
                                genes = as.data.frame(rowData(counts.local.tmm)),
                                group = colData(counts.local.tmm)$CellType)
colnames(y.counts.local.tmm$counts) <- colnames(counts.local.tmm)
rownames(y.counts.local.tmm$samples) <- colnames(counts.local.tmm)
#y.counts.local.tmm$samples$group
design <- model.matrix(~0+group, data=y.counts.local.tmm$samples)
colnames(design) <- sub('group','',colnames(design))
design
y.counts.local.tmm <- estimateDisp(y.counts.local.tmm, design)
saveRDS(y.counts.local.tmm, file='/collignon/Tawdarous/atac_seq_lab/results/csaw_result/new_csaw/paper_csaw/y.counts.local.tmm.rds')
fit.counts.local.tmm <- glmQLFit(y.counts.local.tmm, design, robust=TRUE)
saveRDS(fit.counts.local.tmm, file='/collignon/Tawdarous/atac_seq_lab/results/csaw_result/new_csaw/paper_csaw/fit.counts.local.tmm.rds')


colnames(design) <- gsub("-", "_", colnames(design))
colnames(design) <- make.names(colnames(design), unique = TRUE)


#__________________________t06_BMP.WNT VS t24_BMP.WNT_06h_L _______________________
contr <- makeContrasts(t06_BMP.WNT - t24_BMP.WNT_06h_L, levels = colnames(design))
results <- glmQLFTest(fit.counts.local.tmm, contrast = contr)
rowData(counts.local.tmm) <- cbind(rowData(counts.local.tmm), results$table) # combine GRanges rowdata with differential statistics

# merge nearby windows
# up to "tol" distance apart: 500 bp in this case; max merged window width: 5000 bp
merged.peaks <- mergeWindows(rowRanges(counts.local.tmm), tol=500L, max.width=5000L)
# summary(width(merged.peaks$region))
# should merge some peaks; change as desired

#use most significant window as statistical representation for p-value and FDR for merged windows
tab.best <- getBestTest(merged.peaks$id, results$table)
final.merged.peaks <- GRanges(cbind(as.data.frame(merged.peaks$region), results$table[tab.best$rep.test, -4], tab.best[,-c(7:8)]))

# sort by FDR
final.merged.peaks <- final.merged.peaks[order(final.merged.peaks@elementMetadata$FDR), ]
final.merged.peaks


FDR.thresh <- 0.05 # set as desired
final.merged.peaks.sig <- final.merged.peaks[final.merged.peaks@elementMetadata$FDR < FDR.thresh, ]
z # significant differentially-accessible windows

library(ggplot2)

final.merged.peaks$sig <- "n.s."
final.merged.peaks$sig[final.merged.peaks$FDR < FDR.thresh] <- "significant"

ggplot(data=data.frame(final.merged.peaks),
       aes(x = logCPM, y = logFC, col = factor(sig, levels=c("n.s.", "significant")))) + 
  geom_point() + scale_color_manual(values = c("black", "red")) + 
  geom_smooth(inherit.aes=F, aes(x = logCPM, y = logFC), method = "loess") + # smoothed loess fit; can add span=0.5 to reduce computation load/time
  geom_hline(yintercept = 0) + labs(col = NULL)

upregulated <- final.merged.peaks[final.merged.peaks$direction == "up"]
upregulated <- upregulated[order(upregulated$logCPM),]

anno <- detailRanges(upregulated,
                     txdb = TxDb.Mmusculus.UCSC.mm39.knownGene,
                     orgdb = org.Mm.eg.db, 
                     promoter = c(1500, 1500), 
                     dist = 10 * 1e3L)
anno_df <- as.data.frame(anno)

# Extract relevant metadata from upregulated peaks
metadata_df <- as.data.frame(upregulated)

# Combine the annotation with the metadata
combined_df <- cbind(anno_df, metadata_df[, c("logCPM", "FDR", "direction")])

# View the combined data
combined_df <- combined_df[order(combined_df$logCPM),]
head(combined_df)

# method 4: csaw de novo peaks by local enrichment, csaw loess-normalization
counts.local.loess <- counts.local.filt
counts.local.loess <- normOffsets(counts.local.loess, se.out=TRUE) # type="loess" is now default
colData(counts.local.loess)$CellType <- cell_types

saveRDS(counts.local.loess, file='/collignon/Tawdarous/atac_seq_lab/results/csaw_result/new_csaw/paper_csaw/counts.local.loess.rds')


y.counts.local.loess<- asDGEList(counts.local.loess,
                               genes = as.data.frame(rowData(counts.local.loess)),
                               group = colData(counts.local.loess)$CellType)
colnames(y.counts.local.loess$counts) <- colnames(y.counts.local.loess)
rownames(y.counts.local.loess$samples) <- colnames(y.counts.local.loess)
#y.counts.local.loess$samples$group
design <- model.matrix(~0+group, data=y.counts.local.loess$samples)
colnames(design) <- sub('group','',colnames(design))
design
y.counts.local.loess <- estimateDisp(y.counts.local.loess, design)
saveRDS(y.counts.local.loess, file='/collignon/Tawdarous/atac_seq_lab/results/csaw_result/new_csaw/paper_csaw/y.counts.local.loess.rds')
fit.counts.local.loess <- glmQLFit(y.counts.local.loess, design, robust=TRUE)
saveRDS(fit.counts.local.loess, file='/collignon/Tawdarous/atac_seq_lab/results/csaw_result/new_csaw/paper_csaw/fit.counts.local.loess.rds')


colnames(design) <- gsub("-", "_", colnames(design))
colnames(design) <- make.names(colnames(design), unique = TRUE)


#__________________________t06_BMP.WNT VS t24_BMP.WNT_06h_L _______________________
contr <- makeContrasts(t06_BMP.WNT - t24_BMP.WNT_06h_L, levels = colnames(design))
results <- glmQLFTest(fit.counts.local.loess, contrast = contr)
rowData(counts.local.loess) <- cbind(rowData(counts.local.loess), results$table) # combine GRanges rowdata with differential statistics

# merge nearby windows
# up to "tol" distance apart: 500 bp in this case; max merged window width: 5000 bp
merged.peaks <- mergeWindows(rowRanges(counts.local.loess), tol=500L, max.width=5000L)
# summary(width(merged.peaks$region))
# should merge some peaks; change as desired

#use most significant window as statistical representation for p-value and FDR for merged windows
tab.best <- getBestTest(merged.peaks$id, results$table)
final.merged.peaks <- GRanges(cbind(as.data.frame(merged.peaks$region), results$table[tab.best$rep.test, -4], tab.best[,-c(7:8)]))

# sort by FDR
final.merged.peaks <- final.merged.peaks[order(final.merged.peaks@elementMetadata$FDR), ]
final.merged.peaks


FDR.thresh <- 0.05 # set as desired
final.merged.peaks.sig <- final.merged.peaks[final.merged.peaks@elementMetadata$FDR < FDR.thresh, ]
z # significant differentially-accessible windows

library(ggplot2)

final.merged.peaks$sig <- "n.s."
final.merged.peaks$sig[final.merged.peaks$FDR < FDR.thresh] <- "significant"

ggplot(data=data.frame(final.merged.peaks),
       aes(x = logCPM, y = logFC, col = factor(sig, levels=c("n.s.", "significant")))) + 
  geom_point() + scale_color_manual(values = c("black", "red")) + 
  geom_smooth(inherit.aes=F, aes(x = logCPM, y = logFC), method = "loess") + # smoothed loess fit; can add span=0.5 to reduce computation load/time
  geom_hline(yintercept = 0) + labs(col = NULL)

upregulated <- final.merged.peaks[final.merged.peaks$direction == "up"]
upregulated <- upregulated[order(upregulated$logCPM),]

anno <- detailRanges(upregulated,
                     txdb = TxDb.Mmusculus.UCSC.mm39.knownGene,
                     orgdb = org.Mm.eg.db, 
                     promoter = c(1500, 1500), 
                     dist = 10 * 1e3L)
anno_df <- as.data.frame(anno)

# Extract relevant metadata from upregulated peaks
metadata_df <- as.data.frame(upregulated)

# Combine the annotation with the metadata
combined_df <- cbind(anno_df, metadata_df[, c("logCPM", "FDR", "direction")])

# View the combined data
combined_df <- combined_df[order(combined_df$logCPM),]
head(combined_df)



