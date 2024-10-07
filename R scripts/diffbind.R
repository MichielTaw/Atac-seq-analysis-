if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DiffBind")


library(DiffBind)
library(tidyverse)
library(rtracklayer)
