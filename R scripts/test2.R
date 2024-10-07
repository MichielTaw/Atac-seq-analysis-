# Create a mapping dictionary
chr_map <- list(
  NC_000067.7 = "chr1",
  NC_000068.8 = "chr2",
  NC_000069.7 = "chr3",
  NC_000070.7 = "chr4",
  NC_000071.7 = "chr5",
  NC_000072.7 = "chr6",
  NC_000073.7 = "chr7",
  NC_000074.7 = "chr8",
  NC_000075.7 = "chr9",
  NC_000076.7 = "chr10",
  NC_000077.7 = "chr11",
  NC_000078.7 = "chr12",
  NC_000079.7 = "chr13",
  NC_000080.7 = "chr14",
  NC_000081.7 = "chr15",
  NC_000082.7 = "chr16",
  NC_000083.7 = "chr17",
  NC_000084.7 = "chr18",
  NC_000085.7 = "chr19",
  NC_000086.8 = "chrX",
  NC_000087.8 = "chrY"
)

# Define a vector of accession numbers for mouse chromosomes 1-19, X, and Y
accession_numbers <- names(chr_map)  # Pulls all keys from the chr_map list

# Translate to chromosome names using the map
standard_chr <- sapply(accession_numbers, function(x) chr_map[[x]])

# Print or use the translated chromosome names
print(standard_chr)


-----------------------------------------------------------------------
library(GenomicRanges)
library(csaw)


dir.bam= "/collignon/Tawdarous/atac_seq_lab/results/peak_calling/macs3"
dr.csaw= "/collignon/Tawdarous/atac_seq_lab/results/tests/csaw"


t06_BMP_WNT.rep1=file.path(dir.bam,"t06_BMP+WNT_rep1_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted.bam/t06_BMP+WNT_rep1_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted.bam_peaks.broadPeak")
t06_BMP_WNT.rep2=file.path(dir.bam,"t06_BMP+WNT_rep2_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted.bam/t06_BMP+WNT_rep2_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted.bam_peaks.broadPeak")
t06_WNT.rep1=file.path(dir.bam,"t06_WNT_rep1_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted.bam/t06_WNT_rep1_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted.bam_peaks.broadPeak")
t06_WNT.rep2=file.path(dir.bam,"t06_WNT_rep2_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted.bam/t06_WNT_rep2_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted.bam_peaks.broadPeak")
treat_overlap_peaks.filt.broadPeak=file.path(dr.csaw,"t06_BMP+WNT_overlapPeaks.broadPeak")
control_overlap_peaks.filt.broadPeak=file.path(dr.csaw, "t06_WNT_overlapPeaks.broadPeak")


treat1.peaks <- read.table(t06_BMP_WNT.rep1, sep="\t")[,1:3]
treat2.peaks <- read.table(t06_BMP_WNT.rep2, sep="\t")[,1:3]
control1.peaks <- read.table(t06_WNT.rep1, sep="\t")[,1:3]
control2.peaks <- read.table(t06_WNT.rep2, sep="\t")[,1:3]



colnames(treat1.peaks) <- c("chrom", "start", "end")
colnames(treat2.peaks) <- c("chrom", "start", "end")
colnames(control1.peaks) <- c("chrom", "start", "end")
colnames(control2.peaks) <- c("chrom", "start", "end")

#------------------------------------------------------------------------------

# read naive overlap broadPeak files
treat.overlap.peaks <- read.table(treat_overlap_peaks.filt.broadPeak, sep="\t")[,1:3]
control.overlap.peaks <- read.table(control_overlap_peaks.filt.broadPeak, sep="\t")[,1:3]
colnames(treat.overlap.peaks) <- c("chrom", "start", "end")
colnames(control.overlap.peaks) <- c("chrom", "start", "end")






# convert to GRanges objects
treat1.peaks <- GRanges(treat1.peaks)
treat2.peaks <- GRanges(treat2.peaks)
#treatn.peaks <- GRanges(treatn.peaks)
treat.overlap.peaks <- GRanges(treat.overlap.peaks)
control1.peaks <- GRanges(control1.peaks)
control2.peaks <- GRanges(control2.peaks)
#controln.peaks <- GRanges(controln.peaks)
control.overlap.peaks <- GRanges(control.overlap.peaks)



# define consensus peakset

# one method: union of all replicate peak sets for both conditions
treat.peaks <- union(treat1.peaks, treat2.peaks)
control.peaks <- union(control1.peaks, control2.peaks)
all.peaks <- union(treat.peaks, control.peaks)

# another method: intersect between biological replicates; union between both experimental conditions
treat.peaks <- intersect(treat1.peaks, treat2.peaks)
control.peaks <- intersect(control1.peaks, control2.peaks)
all.peaks <- union(treat.peaks, control.peaks)

blacklist_bed=file.path("/collignon/Tawdarous/data/Tosic19/atac_seq_PRJNA369204/raw_data/Blacklist/lists/mm10-blacklist.v2.bed.gz")  # Path to the blacklist BED file
blacklist <- read.table(blacklist_bed, sep="\t")
colnames(blacklist) <- c("chrom", "start", "end")
blacklist <- GRanges(blacklist)

#specify paired-end BAMs
dir.encode="/collignon/Tawdarous/atac_seq_lab/results/post_alignment/encode"
control1.sorted.noDups.filt.noMT.bam=file.path(dir.encode,"t06_WNT_rep1_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.bam")
control2.sorted.noDups.filt.noMT.bam=file.path(dir.encode,"t06_WNT_rep2_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.bam")
treat1.sorted.noDups.filt.noMT.bam=file.path(dir.encode,"t06_BMP+WNT_rep1_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.bam")
treat2.sorted.noDups.filt.noMT.bam=file.path(dir.encode,"t06_BMP+WNT_rep2_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.bam")
pe.bams <- c(control1.sorted.noDups.filt.noMT.bam, control2.sorted.noDups.filt.noMT.bam,
             treat1.sorted.noDups.filt.noMT.bam, treat2.sorted.noDups.filt.noMT.bam)


# define read parameters


#standard.chr <- c("NC_000067.7", "NC_000068.8", "NC_000069.7") # Use the exact names as found
param <- readParam(max.frag=1000, pe="both")

# count reads in windows specified by MACS2                                      
peak.counts <- regionCounts(pe.bams, all.peaks, param=param)

library("edgeR")
peak.abundances <- aveLogCPM(asDGEList(peak.counts)) 
peak.counts.filt <- peak.counts[peak.abundances > -3, ] # only use peaks logCPM > -3
# few or no peaks should be removed; modify as desired





control1.pe.sizes <- getPESizes(control1.sorted.noDups.filt.noMT.bam)
control2.pe.sizes <- getPESizes(control2.sorted.noDups.filt.noMT.bam)
treat1.pe.sizes <- getPESizes(treat1.sorted.noDups.filt.noMT.bam)
treat2.pe.sizes <- getPESizes(treat2.sorted.noDups.filt.noMT.bam)
gc()
# plot
hist(treat1.pe.sizes$sizes) # repeat for all replicates and conditions
hist(treat2.pe.sizes$sizes)
hist(control1.pe.sizes$sizes)
hist(control2.pe.sizes$sizes)

--------------
# for analysis with csaw de novo enriched query windows, select a window size that is greater than the majority of fragments
  
  ##############################
# count BAM reads in, e.g. 300 bp windows
counts <- windowCounts(pe.bams, width=300, param=param) # set width as desired from the fragment length distribution analyses

# filter uninteresting features (windows) by local enrichment
# local background estimator: 2kb neighborhood
neighbor <- suppressWarnings(resize(rowRanges(counts), width=2000, fix="center")) # change width parameter as desired
wider <- regionCounts(pe.bams, regions=neighbor, param=param) # count reads in neighborhoods
# filter.stat <- filterWindows(counts, wider, type="local") # the filterWindows() function is deprecated and has been replaced by filterWindowsLocal(). This is an archived step.
filter.stat <- filterWindowsLocal(counts, wider)
counts.local.filt <- counts[filter.stat$filter > log2(3),] # threshold of 3-fold increase in enrichment over 2kb neighborhood abundance; change as desired

###############################
# count BAM background bins (for TMM normalization)
binned <- windowCounts(pe.bams, bin=TRUE, width=10000, param=param)

##########################################
# NORMALIZATION

# method 1: MACS2 peaks only, TMM normalization based on binned counts
peak.counts.tmm <- peak.counts.filt
peak.counts.tmm <- normFactors(binned, se.out=peak.counts.tmm)

# method 2: MACS2 peaks only, csaw loess-normalization
peak.counts.loess <- peak.counts.filt
peak.counts.loess <- normOffsets(peak.counts.loess, se.out=TRUE) # type="loess" is now default
# from vignette: "For type="loess", a numeric matrix of the same dimensions as counts, containing the log-based offsets for use in GLM fitting."

# method 3: csaw de novo peaks by local enrichment, TMM normalization based on binned counts
counts.local.tmm <- counts.local.filt
counts.local.tmm <- normFactors(binned, se.out=counts.local.tmm)

# method 4: csaw de novo peaks by local enrichment, csaw loess-normalization
counts.local.loess <- counts.local.filt
counts.local.loess <- normOffsets(counts.local.loess, se.out=TRUE) # type="loess" is now default
# from vignette: "For type="loess", a numeric matrix of the same dimensions as counts, containing the log-based offsets for use in GLM fitting."

#########################################
# DIFFERENTIAL ACCESSIBILITY ANALYSIS

# set working windows for the desired analysis
working.windows <- peak.counts.tmm # MACS2 peaks only, standard TMM normalization based on binned counts
# working.windows <- peak.counts.loess # MACS2 peaks only, for trended biases
# working.windows <- counts.local.tmm # csaw de novo peaks by local enrichment, standard TMM normalization based on binned counts
# working.windows <- counts.local.loess # csaw de novo peaks by local enrichment, for trended biases
# SEE THE CSAW MANUAL FOR MORE INFO ON NORMALIZATION METHODS
###########

# setup design matrix
# see edgeR manual for more information
y <- asDGEList(working.windows)
colnames(y$counts) <- c("control1", "control2", "treat1", "treat2")
rownames(y$samples) <- c("control1", "control2", "treat1", "treat2")
y$samples$group <- c("control", "control", "treat", "treat")
design <- model.matrix(~0+group, data=y$samples)
colnames(design) <- c("control", "treat") # CONFIRM THAT THESE COLUMNS CORRECTLY ALIGN!!
# design
# IMPORTANT: the user should manually confirm that the design matrix is correctly labeled according to sample metadata!

# stabilize dispersion estimates with empirical bayes
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust=TRUE)

# testing for differentially-accessible windows
results <- glmQLFTest(fit, contrast=makeContrasts(treat-control, levels=design))
# head(results$table)
rowData(working.windows) <- cbind(rowData(working.windows), results$table) # combine GRanges rowdata with differential statistics
# working.windows@rowRanges

# merge nearby windows
# up to "tol" distance apart: 500 bp in this case; max merged window width: 5000 bp
merged.peaks <- mergeWindows(rowRanges(working.windows), tol=500L, max.width=5000L)
# summary(width(merged.peaks$region))
# should merge some peaks; change as desired

# use most significant window as statistical representation for p-value and FDR for merged windows
tab.best <- getBestTest(merged.peaks$id, results$table)

# concatenating all relevant statistical data for final merged windows (no redundant columns)
final.merged.peaks <- GRanges(cbind(as.data.frame(merged.peaks$region), results$table[tab.best$rep.test, -4], tab.best[,-c(7:8)]))

# sort by FDR
final.merged.peaks <- final.merged.peaks[order(final.merged.peaks@elementMetadata$FDR), ]
final.merged.peaks
# filter by FDR threshold
FDR.thresh <- 0.05 # set as desired
final.merged.peaks.sig <- final.merged.peaks[final.merged.peaks@elementMetadata$FDR < FDR.thresh, ]
final.merged.peaks.sig # significant differentially-accessible windows

write.table(final.merged.peaks, "treat_vs_control_csaw_DA-windows_all.txt", sep="\t", quote=F, col.names=T, row.names=F)
write.table(final.merged.peaks.sig, "treat_vs_control_csaw_DA-windows_significant.txt", sep="\t", quote=F, col.names=T, row.names=F)



###########################################

# Generate MA plot
library(ggplot2)

final.merged.peaks$sig <- "n.s."
final.merged.peaks$sig[final.merged.peaks$FDR < FDR.thresh] <- "significant"

ggplot(data=data.frame(final.merged.peaks),
       aes(x = logCPM, y = logFC, col = factor(sig, levels=c("n.s.", "significant")))) + 
  geom_point() + scale_color_manual(values = c("black", "red")) + 
  geom_smooth(inherit.aes=F, aes(x = logCPM, y = logFC), method = "loess") + # smoothed loess fit; can add span=0.5 to reduce computation load/time
  geom_hline(yintercept = 0) + labs(col = NULL)