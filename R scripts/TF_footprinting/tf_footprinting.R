#-----------------------------Library---------------------------------------:
library(ATACseqQC)
library(MotifDb)
library(BSgenome.Mmusculus.UCSC.mm39)
library(clusterProfiler)
library(enrichR)
library(BiFET)
library(ATACseqTFEA)
library(ChIPpeakAnno)
genome <- BSgenome.Mmusculus.UCSC.mm39



#-----------------------------PATH
dir.bam= "/collignon/Tawdarous/atac_seq_lab/results/post_alignment/shift_chr_renamed"

t06_BMP_WNT.rep1=file.path(dir.bam,"t06_BMP+WNT_rep1_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted_chr_renamed.sorted.bam")
t06_BMP_WNT.rep2=file.path(dir.bam,"t06_BMP+WNT_rep2_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted_chr_renamed.bam")
t06_BMP_WNT.rep3=file.path(dir.bam,"t06_BMP+WNT_rep3_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted_chr_renamed.bam")

t06_WNT.rep1=file.path(dir.bam,"t06_WNT_rep1_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted_chr_renamed.bam")
t06_WNT.rep2=file.path(dir.bam,"t06_WNT_rep2_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted_chr_renamed.bam")
t06_WNT.rep3=file.path(dir.bam,"t06_WNT_rep3_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted_chr_renamed.bam")

t24_BMP_WNT_06h_LI.rep1=file.path(dir.bam,"t24_BMP+WNT-06h_LI_rep1_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted_chr_renamed.bam")
t24_BMP_WNT_06h_LI.rep2=file.path(dir.bam,"t24_BMP+WNT-06h_LI_rep2_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted_chr_renamed.bam")
t24_BMP_WNT_06h_LI.rep3=file.path(dir.bam,"t24_BMP+WNT-06h_LI_rep3_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted_chr_renamed.bam")

t24_BMP_WNT_06h_L.rep1=file.path(dir.bam,"t24_BMP+WNT-06h_L_rep1_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted_chr_renamed.bam")
t24_BMP_WNT_06h_L.rep2=file.path(dir.bam,"t24_BMP+WNT-06h_L_rep2_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted_chr_renamed.bam")
t24_BMP_WNT_06h_L.rep3=file.path(dir.bam,"t24_BMP+WNT-06h_L_rep3_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted_chr_renamed.bam")

t24_BMP_WNT_06h_LX1.rep1=file.path(dir.bam,"t24_BMP+WNT-06h_LX1_rep1_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted_chr_renamed.bam")
t24_BMP_WNT_06h_LX1.rep2=file.path(dir.bam,"t24_BMP+WNT-06h_LX1_rep2_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted_chr_renamed.bam")
t24_BMP_WNT_06h_LX1.rep3=file.path(dir.bam,"t24_BMP+WNT-06h_LX1_rep3_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted_chr_renamed.bam")

t24_WNT_06h_L_rep1=file.path(dir.bam,"t24_WNT-06h_L_rep1_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted_chr_renamed.bam")
t24_WNT_06h_L_rep2=file.path(dir.bam,"t24_WNT-06h_L_rep2_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted_chr_renamed.bam")
t24_WNT_06h_L_rep3=file.path(dir.bam,"t24_WNT-06h_L_rep3_sorted.rmChrM_RG.marked.sorted.blacklist-filtered.shifted_chr_renamed.bam")


#------------------CTFC FOOTPRINTING--------------------------:
CTCF <- query(MotifDb, c("CTCF"))
CTCF <- as.list(CTCF)
print(CTCF[[1]], digits=2)

seqlev <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", 
            "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chrX", "chrY")



#--------FoxA2------------
FoxA2 <- query(MotifDb, c("FoxA2"))
FoxA2 <- as.list(FoxA2)
print(FoxA2[[1]], digits=2)
#__________t06_BMP_WNT
FoxA2_t06_BMP_WNT.rep1<- factorFootprints(t06_BMP_WNT.rep1, pfm=FoxA2[[1]],
                                           genome=genome,
                                           min.score="90%", seqlev=seqlev,
                                           upstream=100, downstream=100)




featureAlignedHeatmap(FoxA2_t06_BMP_WNT.rep1$signal, 
                      feature.gr=reCenterPeaks(FoxA2_t06_BMP_WNT.rep1$bindingSites,
                                               width=200+width(FoxA2_t06_BMP_WNT.rep1$bindingSites[1])), 
                      annoMcols="score",
                      sortBy="score",
                      n.tile=ncol(FoxA2_t06_BMP_WNT.rep1$signal[[1]]))

FoxA2_t06_BMP_WNT.rep1$spearman.correlation
FoxA2_t06_BMP_WNT.rep1$Profile.segmentation
#________t06_WNT

FoxA2_t06_WNT.rep1<- factorFootprints(t06_WNT.rep1, pfm=FoxA2[[1]],
                                      genome=genome,
                                      min.score="90%", seqlev=seqlev,
                                      upstream=100, downstream=100)
FoxA2_t06_WNT.rep1$spearman.correlation
FoxA2_t06_WNT.rep1$Profile.segmentation

#t24_WNT_06h_L_rep1
FoxA2_t24_WNT_06h_L_rep1<- factorFootprints(t24_WNT_06h_L_rep1, pfm=FoxA2[[1]],
                                      genome=genome,
                                      min.score="90%", seqlev=seqlev,
                                      upstream=100, downstream=100)
FoxA2_t24_WNT_06h_L_rep1$spearman.correlation
FoxA2_t24_WNT_06h_L_rep1$Profile.segmentation




