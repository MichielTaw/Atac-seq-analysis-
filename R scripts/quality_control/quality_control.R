library(ChIPseeker)
library(ATACseqQC)



#-------------------------PATHS-----------------

output.dir <- "/collignon/Tawdarous/atac_seq_lab/results/bamQC_output"

if (!dir.exists(output.dir)) {
  dir.create(output.dir, recursive = TRUE)
}

dir.bam= "/collignon/Tawdarous/atac_seq_lab/results/post_alignment/encode_chr_renamed"


t06_BMP_WNT.rep1=file.path(dir.bam,"t06_BMP+WNT_rep1_sorted.rmChrM_RG.marked.sorted.blacklist-filtered_chr_renamed.bam")
t06_BMP_WNT.rep2=file.path(dir.bam,"t06_BMP+WNT_rep2_sorted.rmChrM_RG.marked.sorted.blacklist-filtered_chr_renamed.bam")
t06_BMP_WNT.rep3=file.path(dir.bam,"t06_BMP+WNT_rep3_sorted.rmChrM_RG.marked.sorted.blacklist-filtered_chr_renamed.bam")

bamQC(t06_BMP_WNT.rep2, outPath = NULL)
bamQC(t06_BMP_WNT.rep3, outPath = NULL)



t06_WNT.rep1=file.path(dir.bam,"t06_WNT_rep1_sorted.rmChrM_RG.marked.sorted.blacklist-filtered_chr_renamed.bam")
t06_WNT.rep2=file.path(dir.bam,"t06_WNT_rep2_sorted.rmChrM_RG.marked.sorted.blacklist-filtered_chr_renamed.bam")
t06_WNT.rep3=file.path(dir.bam,"t06_WNT_rep3_sorted.rmChrM_RG.marked.sorted.blacklist-filtered_chr_renamed.bam")

t24_BMP_WNT_06h_LI.rep1=file.path(dir.bam,"t24_BMP+WNT-06h_LI_rep1_sorted.rmChrM_RG.marked.sorted.blacklist-filtered_chr_renamed.bam")
t24_BMP_WNT_06h_LI.rep2=file.path(dir.bam,"t24_BMP+WNT-06h_LI_rep2_sorted.rmChrM_RG.marked.sorted.blacklist-filtered_chr_renamed.bam")
t24_BMP_WNT_06h_LI.rep3=file.path(dir.bam,"t24_BMP+WNT-06h_LI_rep3_sorted.rmChrM_RG.marked.sorted.blacklist-filtered_chr_renamed.bam")

t24_BMP_WNT_06h_L.rep1=file.path(dir.bam,"t24_BMP+WNT-06h_L_rep1_sorted.rmChrM_RG.marked.sorted.blacklist-filtered_chr_renamed.bam")
t24_BMP_WNT_06h_L.rep2=file.path(dir.bam,"t24_BMP+WNT-06h_L_rep2_sorted.rmChrM_RG.marked.sorted.blacklist-filtered_chr_renamed.bam")
t24_BMP_WNT_06h_L.rep3=file.path(dir.bam,"t24_BMP+WNT-06h_L_rep3_sorted.rmChrM_RG.marked.sorted.blacklist-filtered_chr_renamed.bam")

t24_BMP_WNT_06h_LX1.rep1=file.path(dir.bam,"t24_BMP+WNT-06h_LX1_rep1_sorted.rmChrM_RG.marked.sorted.blacklist-filtered_chr_renamed.bam")
t24_BMP_WNT_06h_LX1.rep2=file.path(dir.bam,"t24_BMP+WNT-06h_LX1_rep2_sorted.rmChrM_RG.marked.sorted.blacklist-filtered_chr_renamed.bam")
t24_BMP_WNT_06h_LX1.rep3=file.path(dir.bam,"t24_BMP+WNT-06h_LX1_rep3_sorted.rmChrM_RG.marked.sorted.blacklist-filtered_chr_renamed.bam")

t24_WNT_06h_L_rep1=file.path(dir.bam,"t24_WNT-06h_L_rep1_sorted.rmChrM_RG.marked.sorted.blacklist-filtered_chr_renamed.bam")
t24_WNT_06h_L_rep2=file.path(dir.bam,"t24_WNT-06h_L_rep2_sorted.rmChrM_RG.marked.sorted.blacklist-filtered_chr_renamed.bam")
t24_WNT_06h_L_rep3=file.path(dir.bam,"t24_WNT-06h_L_rep3_sorted.rmChrM_RG.marked.sorted.blacklist-filtered_chr_renamed.bam")


#----------------fragment size distribution---------

bamfile.labels <- gsub(".bam", "", basename(t06_BMP_WNT.rep3))
fragSize <- fragSizeDist(t06_BMP_WNT.rep3, bamfile.labels)

#-----------Promoter/Transcript body (PT) score)---------------



