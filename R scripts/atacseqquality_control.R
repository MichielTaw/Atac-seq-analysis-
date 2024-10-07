

# Define the refseq to chromosome mapping
refseq_to_chr <- c(
  "NC_000067.7" = "chr1",
  "NC_000068.8" = "chr2",
  "NC_000069.7" = "chr3",
  "NC_000070.7" = "chr4",
  "NC_000071.7" = "chr5",
  "NC_000072.7" = "chr6",
  "NC_000073.7" = "chr7",
  "NC_000074.7" = "chr8",
  "NC_000075.7" = "chr9",
  "NC_000076.7" = "chr10",
  "NC_000077.7" = "chr11",
  "NC_000078.7" = "chr12",
  "NC_000079.7" = "chr13",
  "NC_000080.7" = "chr14",
  "NC_000081.7" = "chr15",
  "NC_000082.7" = "chr16",
  "NC_000083.7" = "chr17",
  "NC_000084.7" = "chr18",
  "NC_000085.7" = "chr19",
  "NC_000086.8" = "chrX",
  "NC_000087.8" = "chrY",
  "NC_005089.1" = "chrM"
)
library(ATACseqQC)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(GenomicAlignments)

path <- "/collignon/Tawdarous/atac_seq_lab/results/post_alignment/encode"
bamFiles <- dir(path, "*.bam$", full.name=TRUE)

# Use NC_000067.7 (RefSeq for chr1) instead of "chr1"
gals <- lapply(bamFiles, function(bamfile){
  readBamFile(bamFile=bamfile, tag=character(0),
              which=GRanges("NC_000067.7", IRanges(1, 1e6)),
              asMates=FALSE)
})

# Modify the transcripts to use RefSeq accessions
txs <- transcripts(TxDb.Mmusculus.UCSC.mm10.knownGene)
seqlevels(txs) <- refseq_to_chr[seqlevels(txs)]

# Plot correlation using NC_000067.7 instead of "chr1"
plotCorrelation(GAlignmentsList(gals), txs, seqlev="NC_000067.7")


_________________________

# Load the BAM files
path <- "/collignon/Tawdarous/atac_seq_lab/results/post_alignment/encode"
bamFiles <- dir(path, "*.bam$", full.name=TRUE)

# Use NC_000067.7 (RefSeq for chr1) instead of "chr1"
gals <- lapply(bamFiles, function(bamfile){
  readBamFile(bamFile=bamfile, tag=character(0),
              which=GRanges("NC_000067.7", IRanges(1, 1e6)),
              asMates=FALSE)
})

# Modify the transcripts to use RefSeq accessions
txs <- transcripts(TxDb.Mmusculus.UCSC.mm10.knownGene)
seqlevelsStyle(txs) <- "NCBI"  # Ensure the seqlevels style matches the RefSeq format

# Map chromosome names to RefSeq accessions in txs
refseq_levels <- names(refseq_to_chr)
txs <- keepSeqlevels(txs, refseq_levels, pruning.mode="coarse")

# Plot correlation using NC_000067.7 instead of "chr1"
plotCorrelation(GAlignmentsList(gals), txs, seqlev="NC_000067.7")


# Verify the path to the BAM files
path <- "/collignon/Tawdarous/atac_seq_lab/results/post_alignment/encode"
bamFiles <- dir(path, pattern = "\\.bam$", full.names = TRUE)

# Check if bamFiles contains valid file paths
print(bamFiles)


library(Rsamtools)

# Load one BAM file to inspect the seqlevels
bamFile <- bamFiles[1]
bamHeader <- scanBamHeader(bamFile)

# Extract seqlevels from the BAM header
bamSeqlevels <- names(bamHeader[[1]]$targets)
print(bamSeqlevels)
