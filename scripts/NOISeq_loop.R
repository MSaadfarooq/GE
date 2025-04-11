# NOISeq loop
library(NOISeq)

load("~/rna_seq/EDA/Outputs/Ens_info_NOISeq.RData") # load ensembl info object

# EXTRACT NOISEQ QC ANNOTATION INFO ----
length <- unique(info[c("gene","length")])
rownames(length) <- length[,1]
length$gene <- NULL

gc <- unique(info[c("gene","percent.gc")])
rownames(gc) <- gc[,1]
gc$gene <- NULL

biotype <- unique(info[c("gene","biotype")])
rownames(biotype) <- biotype[,1]
biotype$gene <- NULL

coords <- unique(info[c("gene","chr","start","end")])
rownames(coords) <- coords[,1]
coords$gene <- NULL

# Read the metadata file
meta <- read.csv("~/rna_seq/Meta_CB.csv")

# Get a list of all gene_count files (assuming they are named like '*_gene_counts.csv')
file_list <- list.files(path = "~/rna_seq", pattern = "*\\_gene_counts.csv$", full.names = TRUE)

# Loop over each gene count file
for (file in file_list) {
  # Extract the tissue name from the file name (without path and extension)
  tissue_name <- gsub("*\\_gene_counts.csv$", "\\1", basename(file))
  
  # Read the current gene count file
  gene_counts <- read.csv(file, row.names = 1) # check files for rownames
  
  # CREATE NOISEQ QC DATASET 
  raw.noiseq <- readData(data=gene_counts, length=length, gc=gc,
                         biotype=biotype, chromosome=coords, factors=meta)
  
  biodetection <- dat(raw.noiseq, k=0, type="biodetection", factor="bw")
  pdf(paste("Biotype_Distribution_", tissue_name, ".pdf", sep = "")) # width=10, height=5
  explo.plot(biodetection, samples=1)
  explo.plot(biodetection, plottype="persample", toplot="protein_coding")
  dev.off()
  
  # COUNT DISTRIBUTION global:
  counts <- dat(raw.noiseq, factor="bw", type="countsbio")
  pdf(paste("Count_Distribution_", tissue_name, ".pdf", sep = ""))
  explo.plot(counts, toplot="global", plottype="boxplot")
  explo.plot(counts, toplot = "protein_coding", plottype = "boxplot")
  explo.plot(counts, samples=1, toplot="global", plottype="boxplot")
  explo.plot(counts, samples=2, toplot="global", plottype="boxplot")
  dev.off()
  
  # SATURATION:
  saturation <- dat(raw.noiseq,type="saturation")
  pdf(paste("DetectedFeatures_vs_SeqDepth_", tissue_name, ".pdf", sep = ""))
  explo.plot(saturation, toplot="global", samples=1:nrow(meta))
  explo.plot(saturation, toplot="protein_coding", samples=1:nrow(meta))
  dev.off()
  
  # PROTEIN CODING COUNT DISTRIBUTION:
  counts.samp <- dat(raw.noiseq, factor=NULL, type="countsbio")
  pdf(paste("ProteinCoding_CountDistribution_", tissue_name, ".pdf", sep = ""),width=50, height=5)
  explo.plot(counts.samp, toplot="protein_coding", samples=NULL, plottype="boxplot")
  dev.off()
  
  # PERCENTAGE LOW COUNT FEATURES PER SAMPLE:
  pdf(paste("ProteinCoding_PercentLowCount_", tissue_name, ".pdf", sep = ""),width=50, height=5)
  explo.plot(counts.samp, toplot="protein_coding", samples=NULL, plottype="barplot")
  dev.off()
  
  # CHECK FOR LENGTH BIAS (significant p-value + R2 > 70%):
  len <- dat(raw.noiseq, factor="bw", type="lengthbias")
  pdf(paste("ProteinCoding_LenghBias_", tissue_name, ".pdf", sep = ""), width=7, height=5)
  explo.plot(len, samples=NULL, toplot="protein_coding")
  dev.off()
  
  # CHECK FOR GC BIAS (significant p-value + R2 > 70%):
  raw.gc <- dat(raw.noiseq, factor="bw", type="GCbias")
  pdf(paste("ProteinCoding_GCBias_", tissue_name, ".pdf", sep = ""))
  explo.plot(raw.gc, samples=NULL, toplot="protein_coding")
  dev.off()
  
}
  
  