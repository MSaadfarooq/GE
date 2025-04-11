# Load libraries ----
suppressPackageStartupMessages({
  library(DESeq2)
  library(vsn)
  library(tidyverse)
  library(apeglm)
  library(ComplexHeatmap)
  library(RColorBrewer)
  library(ggvenn)
  library(BiocParallel)
})
## Simple project ID
project <- "EDA_CB"

# Re-create DESeq2 dataset if the design formula has changed after QC
dds_CB <- DESeqDataSet(CB_gse, design = ~ sex + bw) # variable of interest at the end of formula
dds_CB
levels(dds_CB$bw) # make sure the control level is the first (ref) level, else relevel
dds_CB$bw <- factor(dds_CB$bw, levels = c("normal","low"))

## Run DESeq ----
dds_CB <- DESeq(dds_CB, parallel = TRUE, BPPARAM = MulticoreParam(8))

# Plot dispersions
plotDispEsts(dds_CB, main="Dispersion plot")

# approach-to-count-outliers
pdf("EDA_CB/CooksCutoff_NoBatchCorrection.pdf", width=40, height=7)
boxplot(log10(assays(dds_CB)[["cooks"]]), range=0, las=2)
dev.off()

# dds object transformation
rld_CB <- vst(dds_CB, blind=F) # check vst
saveRDS(rld_CB,file =paste0(project,"_rld.RDS"))
saveRDS(dds_CB,file =paste0(project,"_dds.RDS"))

# Extract DESeq2 results
CB_result <- results(dds_CB, alpha=0.05)
CB_result
resultsNames(dds_CB)
summary(CB_result, alpha=0.05)

head(CB_result[order(CB_result$pvalue), ]) # order by lowest p-value
sum(CB_result$padj < 0.05, na.rm=TRUE)

# Output normalized counts to save as a file
CB_normalized_counts <- counts(dds_CB, normalized=TRUE)
write.csv(CB_normalized_counts, 'CordBlood_normalized_counts.csv', row.names = TRUE)

# distribution of p-values
hist(CB_result$pvalue[CB_result$baseMean > 1], col="tan", main = "Non-adjusted p-value distribution")
# distribution of adjusted p-values
hist(CB_result$padj[CB_result$baseMean > 1], col="lightblue", main = "Adjusted p-value distribution")

# filtering threshold that has been used to filter low count genes
metadata(CB_result)$filterThreshold #genes with basemean < [output] have been filtered

as_tibble(metadata(CB_result)$filterNumRej) %>%
  ggplot(aes(x = theta, y = numRej)) +
  geom_point() +
  geom_vline(xintercept = 0.232, # must be set with basemean threshold
             color = 'red')

## Gene-level filtering ----
# Filter genes by zero expression
CB_result[which(CB_result$baseMean == 0),] %>% 
  data.frame() %>% 
  View()

# Filter genes that have an extreme outlier - summary(CB_result)
CB_result[which(is.na(CB_result$pvalue) & 
                  is.na(CB_result$padj) &
                  CB_result$baseMean > 0),] %>% 
  data.frame() %>% View()

# explore the relationship between LFC, significant DEGs and the genes mean count
plotMA(CB_result)

# Independent Filtering and log-fold shrinkage
CB_result_Lfc <- lfcShrink(dds_CB, coef = "bw_normal_vs_low", type="apeglm", parallel = T, BPPARAM = MulticoreParam(8))
plotMA(CB_result_Lfc)

head(as.data.frame(CB_result))
head(as.data.frame(rowRanges(CB_gse)))

# select the gene to plot by rowname or by numeric index
plotCounts(dds_CB, gene=which.min(CB_result$padj), intgroup=c("sex", 'bw'))

## Extract DE genes ----
# Set thresholds
padj.cutoff <- 0.05
lfc.cutoff <- 1
CB_DElist <- CB_result %>%
  data.frame() %>%
  tibble::rownames_to_column(var="gene") %>% 
  mutate(threshold=case_when(padj<=padj.cutoff & abs(log2FoldChange)>=lfc.cutoff~"Significant",T~"Not Significant")) %>%
  arrange(padj) %>% as_tibble()

# Subset the tibble to keep only significant genes
sigDE_CB <- CB_DElist %>%
  dplyr::filter(threshold=="Significant") # *look for 

# or to get df of sig p-adj
sigDE_CB.df <- as.data.frame(subset(CB_result, padj < 0.05))
sig_DE_genes <- rownames(sigDE_CB.df)
length(sig_DE_genes)

# add annotations like the gene name and biotype ----
library(AnnotationHub)
ah <- AnnotationHub() # create a link to the annotation hub
query(ah,c("EnsDb","112","Homo")) # Search annotationHub for Ensembl annotation

human_ens <- ah[["AH116860"]] #Load the the ensembldb object

# Extract a table of gene annotations or load("~/rna_seq/EDA/Outputs/ENS_anno_genes.Rdata")
anno.genes <- genes(human_ens,return.type = "data.frame") %>% 
  dplyr::select(gene_id,gene_name,gene_biotype,entrezid) %>% 
  as_tibble()
anno.genes

# resolve entrez-Ensembl ID conflict by selecting the first in each list
anno.genes$entrezid <- map(anno.genes$entrezid,1) %>% unlist()
anno.genes

# merge our result tables [sig or all DE] with the annotations to add extra columns
sigDE_anno <- sigDE_CB %>%
  left_join(anno.genes, by = c("gene" = "gene_id")) %>%
  mutate(gene_biotype = as.factor(gene_biotype)) # biotype as a factor

write.csv(sigDE_anno, file = "CB_sigDEGs.csv")
readr::rds(sigDE_anno, 'CB_sigDE_anno.RDS')

##===========================##
# Visualization of DE results
## ==========================##

## Volcano plot ----
CB_DElist %>%
  dplyr::filter(!is.na(padj)) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj),
             color = threshold)) +
  scale_colour_manual(values = c("gray", "red")) +
  geom_point(size = 0.5) +
  geom_hline(yintercept = -log10(0.05)) +
  geom_vline(xintercept = 1) +
  geom_vline(xintercept = -1) + # lines optional
  theme_bw()

## Barplot of gene_biotype
ggplot(sigDE_anno,aes("Significant genes",fill=gene_biotype)) +
  geom_bar() + 
  theme_bw() + 
  xlab("") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

## Heatmap of results ----
# Transform counts
rld_CB <- rlog(dds, blind=F)

# Get top 300 DE genes
HM_genes <- CB_result[order(CB_result$padj), ] |>
  head(300) |>
  rownames()
heatmapData <- assay(rld_CB)[HM_genes, ]

# Scale counts for visualization
heatmapData <- t(scale(t(heatmapData)))

# Add annotation & colors
heatmapColAnnot <- data.frame(colData(rld_CB)[, c("bw", "sex")])
heatmapColAnnot <- HeatmapAnnotation(df = heatmapColAnnot)
colors <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(255)

# Plot as Heatmap
ComplexHeatmap::Heatmap(heatmapData,
                        col = colors,
                        top_annotation = heatmapColAnnot,
                        cluster_rows = TRUE, show_row_names = FALSE,
                        cluster_columns = T, column_dend_reorder = T)

# Top 10 differentially expressed genes (by padj values) ----
# Order results by padj values
top10_sigOE <- CB_DElist %>% 
  arrange(padj) %>% 	# Arrange rows by padj
  #pull(gene) %>% 		# Extract character vector of ordered genes
  head(10)		# Extract the first 10 genes

meta<-coldata %>% rename(sample=names)
  
as_tibble(counts(dds[top10_sigOE$gene, ]),
          rownames = 'gene') %>%
  pivot_longer(names_to = "sample", values_to = "counts", -gene) %>%  
  left_join(as_tibble(coldata, rownames = "sample")) %>% 
  ggplot(aes(x = sample, y = counts, fill = bw)) + #geom_point()
  geom_bar(stat = 'identity', color = "gray30") +
  facet_wrap( ~ gene, scales = "free", ncol = 3) +
  theme(axis.text.x = element_text(size = 7, angle = 90),
        axis.title.x = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7))

# Venn diagram of the results ----
# prepare the data with a column for Gene
PL_result_vn <- as.data.frame(PL_result) %>%
     rownames_to_column("Gene")
CB_result_vn <- as.data.frame(CB_result) %>%
     rownames_to_column("Gene")

# function to get genes
getGenes <- function(shrTab, direction) {
  sign <- ifelse(direction == "up", 1, -1)
  shrTab %>%
    filter(padj < 0.05) %>%
    filter(sign * log2FoldChange > 0) %>%
    pull("Gene")
}
vennList <- list(Upregulated_CB = getGenes(CB_result_vn, "up"),
                 Downregulated_CB = getGenes(CB_result_vn, "down"),
                 Upregulated_PL = getGenes(PL_result_vn, "up"),
                 Downregulated_PL = getGenes(PL_result_vn, "down"))
# draw venn
ggvenn(vennList, set_name_size = 3)

# common significant DE genes in both groups
genes_PL_CB <- inner_join(sigDE_PL, sigDE_CB, by = "gene")

sigDE_annoPL_CB <- genes_PL_CB %>%
  inner_join(anno.genes, by = c("gene" = "gene_id"))
write.csv(sigDE_annoPL_CB, "DE_sig_PL-CB.csv")

# plot their Venn Diagram
sigVen <- list(CordBlood=sigDE_CB$gene, Placenta=sigDE_PL$gene)
ggvenn(sigVen, fill_color = c("lightgreen", "purple3"), auto_scale = T, stroke_color = FALSE )
