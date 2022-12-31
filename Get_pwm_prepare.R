if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install('Biostrings')
BiocManager::install('TFBSTools')
BiocManager::install('JASPAR2020')