# set environment
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggpubr)
  library(grid)
  library(gridExtra)
  library(data.table)
  library(cowplot)
  library(gplots)
  library(dendextend)
  library(reshape2)
  library(igraph)
  library(STRINGdb)
})
options(repr.plot.width =20, repr.plot.height = 17, repr.plot.resolution=300)

# load data 
meta <- read_csv("~/treat/meta.csv")
sigDEGs <- read_csv("~/treat/EDA_DE/sigDEGs.csv")

normalized_counts <- read_csv("~/treat/EDA_DE/Outputs/normalized_counts.csv")
normalized_counts <- column_to_rownames(normalized_counts, var = '...1')
all(colnames(normalized_counts) %in% meta$sample)

# select genes for cluster
selec <- as.list(sigDEGs$gene)
total_selec <- list()
total_selec <- append(total_selec, selec)
total_selec <- c(unique(total_selec))
total_selec <- t(as.data.frame(total_selec))
DEgenes <- normalized_counts[total_selec[,1],]

# scale the data
DEgenes <- as.matrix(DEgenes)
scaledata <- t(scale(t(DEgenes)))

# Hierarchical tree of the samples
hc <- hclust(as.dist(1-cor(scaledata, method="spearman")), method="complete")

sampleTree = as.dendrogram(hc, method="average")
plot(sampleTree,
     main = "Sample Clustering",
     ylab = "Height")

# Hierarchical tree of the genes
hr <- hclust(as.dist(1-cor(t(scaledata), method="pearson")), method="complete")

geneTree = as.dendrogram(hr, method="average")
plot(geneTree,
     leaflab = "none",             
     main = "Gene Clustering",
     ylab = "Height")
# construct the heatmap
heatmap.2(DEgenes,
          Rowv=as.dendrogram(hr), 
          Colv=as.dendrogram(hc),
          col=redgreen(100),
          scale="row",
          margins = c(7, 7),
          cexCol = 0.7,
          labRow = F,
          main = "Heatmap of DEGs",
          trace = "none")
# extract discrete clusters of genes as a means to identify co-expression modules
hclusth1.5 = cutree(hr, h=1.5) #cut tree at height of 1.5
hclusth1.0 = cutree(hr, h=1.0) #cut tree at height of 1.0
hclusth0.5 = cutree(hr, h=0.5) #cut tree at height of 0.5

# number of clusters
length(unique(hclusth0.5))
length(unique(hclusth1.5))
head(hclusth0.5)

#plot the tree
plot(geneTree,
     leaflab = "none",
     main = "Gene Clustering",
     ylab = "Height")

#add the three cluster vectors
the_bars <- cbind(hclusth0.5, hclusth1.0, hclusth1.5)

colored_bars(the_bars, geneTree, sort_by_labels_order = T, y_shift=-0.1,
             rowLabels = c("h=0.5","h=1.0","h=1.5"),cex.rowLabels=0.7)
# add lines showing the cut heights
abline(h=1.5, lty = 2, col="grey")
abline(h=1.0, lty = 2, col="grey")
abline(h=0.5, lty = 2, col="grey")

# alternatively it is also possible to set the number of cluster you want.
hclustk3 = cutree(hr, k=3)

plot(geneTree,
     leaflab = "none",
     main = "Gene Clustering",
     ylab = "Height")
colored_bars(hclustk3, geneTree, sort_by_labels_order = T, y_shift=-0.1, rowLabels = c("k=5"),cex.rowLabels=0.7)

# visualize the bar indicating the clusters in combination with a heatmap
clustColBar <- rainbow(length(unique(hclustk3)), start=0.1, end=0.9)
clustColBar <- clustColBar[as.vector(hclustk3)]

heatmap.2(DEgenes,
          Rowv=as.dendrogram(hr), 
          Colv=as.dendrogram(hc),
          col=redgreen(100),
          scale="row",
          margins = c(7, 7),
          cexCol = 0.7,
          labRow = F,
          main = "Heatmap.2",
          trace = "none",
          RowSideColors=clustColBar,
          key = FALSE)

# Vizualizing the cores of the clusters ----
clust.core = function(i, dat, clusters) {
  ind = (clusters == i)
  colMeans(dat[ind,])
}
clusters <- hclustk3
cores <- sapply(unique(clusters), clust.core, scaledata, clusters)

moltenCores <- melt(cores) ## get the data in a "long" format
colnames(moltenCores) <- c('sample','cluster','value')

ggplot(moltenCores, aes(x=sample,y=value, group=cluster, colour=as.factor(cluster))) + 
  geom_point() + geom_line() +
  xlab("Sample") + ylab("Expression") +
  labs(title= "Cluster Expression of the samples",color = "Cluster") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('../treat/EDA_DE/Outputs/co-exp_clusters_DEGs.png')

# Isolate a cluster of interrest.
clust2 <- t(scaledata[hclustk3==2,])

#get the data frame into long format for plotting
clust2Molten <- melt(clust2, id.vars = "Sample")
colnames(clust2Molten) <- c('sample','gene','value')

#Subset the cores molten dataframe so we can plot core1
core <- moltenCores[moltenCores$cluster==2,]

ggplot(clust2Molten, aes(x=sample,y=value)) + 
  geom_line(color="grey", aes(color="grey", group=gene)) +
  geom_line(data=core, aes(sample,value, group=cluster), color="blue",inherit.aes=FALSE) + #adds the core 
  geom_point(data=core, aes(sample,value, group=cluster), color="blue",inherit.aes=FALSE) +
  xlab("Samples") + ylab("Expression") +
  labs(title= paste0("Cluster 2 consisting ", ncol(clust2), " genes"),color = "Score") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('../treat/EDA_DE/Outputs/cluster2_Genes.png')

# How many clusters to choose ----
wss <- (nrow(scaledata)-1)*sum(apply(scaledata,2,var))
for (i in 2:20) wss[i] <- sum(kmeans(scaledata,
                                     centers=i)$withinss)
plot(1:20, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

library(cluster)
sil <- rep(0, 20)
#repeat k-means for 1:20 and extract silhouette:
for(i in 2:20){
  k1to20 <- kmeans(scaledata, centers = i, nstart = 25, iter.max = 20)
  ss <- silhouette(k1to20$cluster, dist(scaledata))
  sil[i] <- mean(ss[, 3])
}
# Plot the  average silhouette width
plot(1:20, sil, type = "b", pch = 19, xlab = "Number of clusters k", ylab="Average silhouette width")
abline(v = which.max(sil), lty = 2)

cat("Average silhouette width optimal number of clusters:", which.max(sil), "\n")

# subset DEGs for string network ----
clust2_DEGs <- sigDEGs %>% # filter cluster2 genes
  filter(gene %in% colnames(clust2)) %>% select(gene_name, log2FoldChange, padj)
clust2_DEGs <- as.data.frame(clust2_DEGs)

# creating an instance of the STRINGdb reference class
string_db <- STRINGdb$new(version="12.0", species=9606,
                          score_threshold=400, network_type="full",
                          link_data='combined_only', input_directory="")

# map the gene names to the STRING database identifiers - adds an additional column
sigDEGs_map <- string_db$map(clust2_DEGs, "gene_name", removeUnmappedRows = TRUE)

# extract the most significant genes & produce the STRING network
hits <- sigDEGs_map$STRING_id # or subset[1:200]
string_db$plot_network(hits)

# get the STRING identifier of proteins of interest
p53 = string_db$mp( "tp53" )
atm = string_db$mp( "atm" )

# retrieve the interactions that connect input proteins
string_db$get_interactions(c(p53, atm)) #string_db$get_interactions(sigDEGs_map$STRING_id)
string_db$plot_network(c(p53, atm))

# color in green the genes: down-regulated, red genes: up-regulated
sigDEGs_map_UPdown <- string_db$add_diff_exp_color(subset(sigDEGs_map,log10(padj) >= -log10(0.05) | abs(log2FoldChange) >= 0.5), logFcColStr ="log2FoldChange")
head(sigDEGs_map_UPdown)
table(sigDEGs_map_UPdown$color)
#sigDEGs_map_UPdown <- string_db$add_diff_exp_color(screen = sigDEGs_map, logFcColStr = 'log2FoldChange')                                                   

payload_id <- string_db$post_payload(sigDEGs_map_UPdown$STRING_id,
                                     colors = sigDEGs_map_UPdown$color)
# display a STRING network png with the "halo"
string_db$plot_network(hits, payload_id=payload_id)

# Functional enrichment ----
enrichmentGO <- string_db$get_enrichment(hits, category = "Process") # GO for top 200

#KEGG enrichment analysis
enrichmentKEGG <- string_db$get_enrichment(sigDEGs_map$STRING_id, category = "KEGG")

# bar plot terms
ggplot(subset(enrichmentGO, category=='Process'),
       aes(number_of_genes, description, fill = p_value))+
  geom_bar(stat = 'identity')+
  geom_text(aes(label = number_of_genes, hjust = -0.5)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10))

go_immune <- subset(enrichmentGO, description == 'Regulation of immune response')

# get subnetwork from GO analysis
string_id <- str_split(go_immune$inputGenes, pattern=',')
gene_sym <- str_split(go_immune$preferredNames, pattern=',')
go_immune <- data.frame(nodes=string_id, gene_sym=gene_sym)
colnames(go_immune) <- c('nodes','gene_symbol')
head(go_immune)
# plot this network
string_db$plot_network(go_immune$nodes)

# get interactions of this network
sub_edges <- string_db$get_interactions(go_immune$nodes)
sub_edges <- left_join(sub_edges, go_immune, by = c('from'='nodes'))
sub_edges <- left_join(sub_edges, go_immune, by = c('to'='nodes'))
sub_edges
sub_edges <- unique(sub_edges)
sub_edges <- sub_edges |> select(gene_symbol.x, gene_symbol.y) |>
  rename(from=gene_symbol.x, to=gene_symbol.y)

# igraph to plot ----
g <- graph_from_data_frame(sub_edges, directed = FALSE, vertices = go_immune$gene_symbol)
vcount(g) # how many proteins

sort(degree(g), decreasing = T) # top proteins with degree - most connected nodes
set.seed(122)
plot(g)

# select top connected proteins, hub gene will be 'gold'
go_immune$color <- ifelse(grepl('CCR2|FCGR3B|FCGR1A|MNDA',go_immune$gene_symbol),'gold','grey')
go_immune

V(g)$color <- go_immune$color
set.seed(122)
plot(g)
# can add size
go_immune$size <- ifelse(grepl('CCR2|FCGR3B|FCGR1A|MNDA',go_immune$gene_symbol),20,15)
go_immune

V(g)$size <- go_immune$size
set.seed(122)
plot(g, edge.color='orange2', layout=layout_on_sphere)
plot(g, edge.color='orange2', layout=layout_with_kk)

# clustering
ceb <- cluster_edge_betweenness(g)
plot_dendrogram(ceb, mode = 'hclust')

plot(ceb, g)

# know terms are assigned to your set of proteins (and not necessary enriched)
sigDEG_annot <- string_db$get_annotations(hits)

# get some clusters by igraph
clustersList <- string_db$get_clusters(sigDEGs_map$STRING_id[1:500])

for(i in seq(1:4)){   # plot first 4 clusters, replace i for any other number
string_db$plot_network(clustersList[[i]])
}
string_db$plot_network(clustersList[[4]])
