# load packages 
library(tidyverse)
library(NOISeq)
library(RColorBrewer)
library(biomaRt)

# load data [load info BM object from PC]
gene_counts <- read.csv("Outputs/gene_counts.csv", row.names = 1)
#rownames(gene_counts) <- gene_counts$X
#gene_counts$X <- NULL
meta <- read.csv("~/rna_seq/RNA-seq/data/Metadata.csv")

meta <- meta[match(colnames(gene_counts), meta$cordBlood_ID),]

# get ensemble annotation
ensembl <- useEnsembl(biomart="ensembl",version=112,dataset="hsapiens_gene_ensembl")

# Pull out information of interest
chrs <- c(1:22,"X","Y")
info <- getBM(attributes=c("external_gene_name", "ensembl_gene_id", "chromosome_name", "start_position", "end_position", "percentage_gene_gc_content", "gene_biotype"), mart=ensembl)

info$length <- info$end_position - info$start_position
# rename columns
info <- dplyr::rename(info, c("gene.symbol"="external_gene_name","gene"="ensembl_gene_id","chr"="chromosome_name","start"="start_position","end"="end_position","percent.gc"="percentage_gene_gc_content","biotype"="gene_biotype"))

info <- subset(info, chr %in% chrs)
# save all above
save(info,file="Outputs/Ens_info_NOISeq.RData")
load("~/rna_seq/EDA/Outputs/Ens_info_NOISeq.RData") # alternatively

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

# CREATE NOISEQ QC DATASET 
raw.noiseq <- readData(data=gene_counts, length=length, gc=gc, biotype=biotype, chromosome=coords, factors=meta)

# BIOTYPE DISTRIBUTION:
# https://www.rdocumentation.org/packages/NOISeq/versions/2.16.0/topics/Biodetection

biodetection <- dat(raw.noiseq, k=0, type="biodetection", factor="bw")

pdf("EDA_PL/BiotypeDistribution.pdf", width=10, height=5)
explo.plot(biodetection, samples=1)
explo.plot(biodetection, plottype="persample", toplot="protein_coding")
dev.off()

# BIOTYPE DISTRIBUTION / sample of each group
biodetection.comp <- dat(raw.noiseq, k=0, type="biodetection", factor=NULL)
explo.plot(biodetection.comp, samples = c(1, 4), toplot="protein_coding", plottype="comparison") 

# COUNT DISTRIBUTION global:
counts <- dat(raw.noiseq, factor="bw", type="countsbio")

pdf("EDA_PL/CountDistribution.pdf", width=10, height=5)
explo.plot(counts, toplot="global", plottype="boxplot")
explo.plot(counts, toplot = "protein_coding", plottype = "boxplot")
explo.plot(counts, samples=1, toplot="global", plottype="boxplot")
explo.plot(counts, samples=2, toplot="global", plottype="boxplot")
dev.off()

# SATURATION:
saturation <- dat(raw.noiseq,type="saturation")

pdf("EDA_PL/DetectedFeatures_vs_SeqDepth.pdf", width=10, height=10)
explo.plot(saturation, toplot="global", samples=1:nrow(meta))
explo.plot(saturation, toplot="protein_coding", samples=1:nrow(meta))
dev.off()

# PROTEIN CODING COUNT DISTRIBUTION:
counts.samp <- dat(raw.noiseq, factor=NULL, type="countsbio")

pdf("EDA_PL/ProteinCoding_CountDistribution.pdf", width=50, height=5)
explo.plot(counts.samp, toplot="protein_coding", samples=NULL, plottype="boxplot")
dev.off()

# PERCENTAGE LOW COUNT FEATURES PER SAMPLE:
pdf("EDA_PL/ProteinCoding_PercentLowCount.pdf", width=50, height=5)
explo.plot(counts.samp, toplot="protein_coding", samples=NULL, plottype="barplot")
dev.off()

# CHECK FOR LENGTH BIAS (significant p-value + R2 > 70%):
# https://www.rdocumentation.org/packages/NOISeq/versions/2.16.0/topics/lengthbias
len <- dat(raw.noiseq, factor="bw", type="lengthbias")

pdf("EDA_PL/ProteinCoding_LenghBias.pdf", width=7, height=5)
explo.plot(len, samples=NULL, toplot="protein_coding")
dev.off()

# CHECK FOR GC BIAS (significant p-value + R2 > 70%):
raw.gc <- dat(raw.noiseq, factor="bw", type="GCbias")

pdf("EDA_PL/ProteinCoding_GCBias.pdf", width=7, height=5)
explo.plot(raw.gc, samples=NULL, toplot="protein_coding")
dev.off()
