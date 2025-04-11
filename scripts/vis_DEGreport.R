# additional visualizations
library(pheatmap)
library(DEGreport)
library(ggrepel)
library(EnhancedVolcano)

# load anno.genes
load("~/rna_seq/EDA/Outputs/ENS_anno_genes.Rdata")
#meta <- meta |> dplyr::filter(cordBlood_ID %in% colnames(normalized_counts))

# get normalized_counts_anno
normalized_counts_anno <- counts(dds, normalized=T) %>% 
  data.frame() %>%
  rownames_to_column(var="gene") %>%
  as_tibble() %>%
  left_join(anno.genes, by=c("gene" = "gene_id"))

### Extract normalized expression for significant genes from the OE and control samples (2:4 and 7:9)
norm_OEsig <- normalized_counts_anno %>% 
  dplyr::filter(gene %in% sigDE$gene)  

### Set a color palette
heat_colors <- brewer.pal(6, "YlOrRd")

### Run pheatmap using the metadata data frame for the annotation
pheatmap(norm_OEsig[2:12], 
         color = heat_colors, 
         cluster_rows = T, 
         show_rownames = F,
         annotation = heatmapColAnnot, 
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         height = 20)

# Size factor QC
degCheckFactors(normalized_counts)

# Mean-Variance QC plots
degQC(normalized_counts, meta[["bw"]], pvalue = res_LFC[["pvalue"]])

# Covariates effect on count data
resCov <- degCovariates(log2(counts(dds)+0.5), colData(dds))

# Plot top genes coloring by group
degPlot(dds = dds, res = res_LFC, n = 6, xs = "bw", group = "sex")

# plotting genes in a wide format
degPlotWide(dds, rownames(dds)[1:6], group="bw")

# volcano plot
degVolcano(data.frame(resLFC_tb[,c("log2FoldChange","padj")]), # table - 2 columns
           plot_text = data.frame(resLFC_tb[1:10,c("log2FoldChange","padj","gene")])) # table to add names

resreport <- degResults(dds = dds, name = "test", org = NULL,
                        do_go = FALSE, group = "bw", xs = "bw",
                        path_results = NULL)

cor <- degCorCov(colData(dds))

filter_count <- degFilter(counts(dds),
                          meta, "bw",
                          min=1, minreads = 50)
cat("gene in final count matrix", nrow(filter_count))

# plot top 20 genes
top20_sigOE_genes <- resLFC_tb %>% 
  arrange(padj) %>% 	# Order results by padj values
  pull(gene) %>% 		#Extract character vector of ordered genes
  head(n=20)		#Extract the first 20 genes

## normalized counts for top 20 significant genes
top20_sigOE_norm <- normalized_counts_anno %>%
  dplyr::filter(gene %in% top20_sigOE_genes)

# Gathering the columns to have normalized counts to a single column
top20_sigOE_long <- top20_sigOE_norm %>% pivot_longer(colnames(top20_sigOE_norm)[2:12], names_to = "cordBlood_ID", values_to = "normalized_counts")

top20_sigOE_long <- inner_join(meta, top20_sigOE_long)

# Save plotcounts to a data frame object
d <- plotCounts(dds, gene="ENSG00000155363", intgroup="bw", returnData=TRUE)

# What is the data output of plotCounts()?
d %>% View()

# Plot the normalized counts, using the samplenames (rownames(d) as labels)
ggplot(d, aes(x = bw, y = count, color = bw)) + 
  geom_point(position=position_jitter(w = 0.1,h = 0)) +
  geom_text_repel(aes(label = rownames(d))) + 
  theme_bw() +
  ggtitle("Birthweight gene of interest") +
  theme(plot.title = element_text(hjust = 0.5))

# Labelled volcano plot
EnhancedVolcano(res_LFC,
                lab = rownames(res_LFC),
                title = 'CB GE', subtitle = '',
                x = 'log2FoldChange',
                y = 'pvalue', # or padj
                axisLabSize = 10,
                labSize = 3.0,
                pointSize = 1.0,
                legendPosition = "right",
                legendLabSize = 12,
                caption = '', 
                col = c("grey30", "lightblue", "royalblue", "red3"))

# heatmap
raw.dist <- dist(t(assay(vsd)))
raw.distMatrix <- as.matrix(raw.dist)
colz <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255) ## Set a colour pallette in shades of blue
pheatmap(raw.distMatrix,
         clustering_distance_rows=raw.dist,
         clustering_distance_cols=raw.dist,
         col=colz)

# Instead of using GeneRatio, use the fraction of DE genes in the gene sets which are like a “scaled” value for all gene sets
n_11 = resGOTable$Count
n_01 = as.numeric(gsub("/.*$", "", resGOTable$BgRatio))
resGOTable$DE_Ratio = n_11/n_01
ggplot(resGOTable[1:10, ], 
       aes(x = zScore, y = factor(Description, levels = rev(Description)), 
           col = DE_Ratio, size = Count)) +
  geom_point() +
  ylab("")