# install BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# List of packages to be installed
packages <- c(
  "ATACseqQC",
  "ChIPpeakAnno",
  "GenomicAlignments",
  "BSgenome.Mmusculus.UCSC.mm10",
  "TxDb.Mmusculus.UCSC.mm10.knownGene",
  "BSgenome.Mmusculus.UCSC.mm39",
  "TxDb.Mmusculus.UCSC.mm39.knownGene",
  "MotifDb",
  "motifStack"
)
BiocManager::install(packages)

lapply(packages, library, character.only = TRUE)


bamFile = "mtawdarous@rstudio2:/collignon/Tawdarous/atac_seq_lab/results/alignment/rmChrM__AND_duplicate_removed/t24_WNT-06h_L_rep2_sorted.rmChrM_RG.marked.bam"
bamQC(bamFile, outPath = NULL)


print
Sys.setenv(CFLAGS="-I/home/mtawdarous/softwares/gsl/include")
Sys.setenv(LDFLAGS="-L/home/mtawdarous/softwares/gsl/lib")


export PATH="/home/mtawdarous/softwares/gsl/bin:$PATH"
export LD_LIBRARY_PATH="/home/mtawdarous/softwares/gsl/lib:$LD_LIBRARY_PATH"
export C_INCLUDE_PATH="/home/mtawdarous/softwares/gsl/include:$C_INCLUDE_PATH"


PKG_CPPFLAGS=-I/home/mtawdarous/softwares/gsl/include
PKG_LIBS=-L/home/mtawdarous/softwares/gsl/lib -lgsl -lgslcblas




Sys.setenv(PATH = paste(Sys.getenv("PATH"), "/home/linuxbrew/.linuxbrew/opt/gsl/bin", sep=":"))
Sys.setenv(LD_LIBRARY_PATH = paste(Sys.getenv("LD_LIBRARY_PATH"), "/home/linuxbrew/.linuxbrew/opt/gsl/lib", sep=":"))
Sys.setenv(LD_LIBRARY_PATH = paste(Sys.getenv("LD_LIBRARY_PATH"), "/home/linuxbrew/.linuxbrew/opt/gsl/lib", sep=":"))
LD_LIBRARY_PATH=/usr/local/lib
export LD_LIBRARY_PATH


BiocManager::install("ATACseqQC")




library(readr)

# Read the delimited file with a specified separator
peak <- read_delim(file = "/collignon/Tawdarous/atac_seq_lab/results/peak_analysis/peak_intersections.txt", delim = "\t")

# Display the data
print(peak)



# Load necessary libraries
library(ggplot2)
library(readr)
library(tidyr)
library(dplyr)

# Read the data

# Print the data
print(peak)

# Melt the data for plotting
peak_long <- peak %>%
  pivot_longer(cols = c(Genrich, MACS3, HMMRATAC, Genrich_MACS3, Genrich_HMMRATAC, MACS3_HMMRATAC, All_Intersect),
               names_to = "Category",
               values_to = "Count")

# Plot the peak counts for each software
ggplot(peak_long %>% filter(Category %in% c("Genrich", "MACS3", "HMMRATAC")), aes(x = Replicate, y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Peak Counts by Software", x = "Replicate", y = "Count")

# Plot the intersection counts between pairs of software
ggplot(peak_long %>% filter(Category %in% c("Genrich_MACS3", "Genrich_HMMRATAC", "MACS3_HMMRATAC")), aes(x = Replicate, y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Intersection Counts by Software Pair", x = "Replicate", y = "Count")

# Plot the intersection count of all three software
ggplot(peak_long %>% filter(Category == "All_Intersect"), aes(x = Replicate, y = Count)) +
  geom_bar(stat = "identity", fill = "blue") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Intersection Count of All Three Software", x = "Replicate", y = "Count")

# Plot the intersection percentage
ggplot(peak, aes(x = Replicate, y = `Intersection%`)) +
  geom_bar(stat = "identity", fill = "green") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Intersection Percentage", x = "Replicate", y = "Intersection Percentage (%)")


peak2 <- read_delim(file = "/collignon/Tawdarous/atac_seq_lab/results/peak_analysis/peak_intersections2.txt", delim = "\t")

peak3 <- read_delim(file = "/collignon/Tawdarous/atac_seq_lab/results/peak_analysis/peak_intersections3.txt", delim = "\t")

