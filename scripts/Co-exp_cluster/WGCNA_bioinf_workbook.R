library(tidyverse)
library(WGCNA)
library(DESeq2)

# import data
CB_gene_counts <- read_csv("EDA_CB/CB_gene_counts.csv")

names(CB_gene_counts)[1] = "GeneId"
names(CB_gene_counts)

# clean and tidy the dataset
col_sel = names(CB_gene_counts)[-1]     # Get all but first column name
mdata <- CB_gene_counts %>%
  tidyr::pivot_longer(
    .,                        # The dot is the the input data, magrittr tutorial
    col = all_of(col_sel)
  ) %>%  
  mutate(
    group = gsub("-.*","", name) %>% gsub("[.].*","", .)   # Get the shorter treatment names
  )

# ==== Plot groups (Sample Groups vs RNA Seq Counts) to identify outliers
(
  p <- mdata %>%
    ggplot(., aes(x = name, y = value)) +             # x = treatment, y = RNA Seq count
    geom_violin() +                                   # violin plot, show distribution
    geom_point(alpha = 0.2) +                         # scatter plot
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90)          # Rotate treatment text
    ) +
    labs(x = "Treatment Groups", y = "RNA Seq Counts") +
    facet_grid(cols = vars(group), drop = TRUE, scales = "free_x")      # Facet by hour
)
# Prepare DESeq input, which is expecting a matrix of integers
de_input = as.matrix(CB_gene_counts[,-1])
row.names(de_input) = CB_gene_counts$GeneId
de_input[1:5,1:10]

dds <- DESeqDataSetFromMatrix(round(de_input),
                              Meta_CB,
                              design = ~bw)
dds <- DESeq(dds)

vsd <- vst(dds)

wpn_vsd <- getVarianceStabilizedData(dds)
rv_wpn <- rowVars(wpn_vsd)
summary(rv_wpn)

q75_wpn <- quantile( rowVars(wpn_vsd), .75)  # <= original
q95_wpn <- quantile( rowVars(wpn_vsd), .95)  # <= changed to 95 quantile to reduce dataset
expr_normalized <- wpn_vsd[ rv_wpn > q95_wpn, ]
dim(expr_normalized)

expr_normalized_df <- data.frame(expr_normalized) %>%
  mutate(
    Gene_id = row.names(expr_normalized)
  ) %>%
  pivot_longer(-Gene_id)

expr_normalized_df %>% ggplot(., aes(x = name, y = value)) +
  geom_violin() +
  geom_point() +
  theme_bw() +
  theme(
    axis.text.x = element_text( angle = 90)
  ) +
  ylim(0, NA) +
  labs(
    title = "Normalized and 95 quantile Expression",
    x = "treatment",
    y = "normalized expression"
  )
# transpose the data
input_mat = t(expr_normalized)

input_mat[1:5,1:10]

allowWGCNAThreads(nThreads = 32)


gsg <-goodSamplesGenes(expression.data)

summary(gsg)
gsg$allOK


# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to = 20, by = 2))

# Call the network topology analysis function
spt = pickSoftThreshold(
  input_mat,             # <= Input data
  #blockSize = 30,
  powerVector = powers,
  verbose = 5
)

par(mfrow = c(1,2));
cex1 = 0.9;

plot(spt$fitIndices[, 1],
     -sign(spt$fitIndices[, 3]) * spt$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence")
)
text(spt$fitIndices[, 1],
     -sign(spt$fitIndices[, 3]) * spt$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red"
)
abline(h = 0.90, col = "red")
plot(spt$fitIndices[, 1],
     spt$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)
text(spt$fitIndices[, 1],
     spt$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red")

picked_power = 9
temp_cor <- cor       
cor <- WGCNA::cor         # Force it to use WGCNA cor function (fix a namespace conflict issue)

netwk <- blockwiseModules(input_mat,                # <= input here
                          
                          # == Adjacency Function ==
                          power = picked_power,                # <= power here
                          networkType = "signed",
                          
                          # == Tree and Block Options ==
                          deepSplit = 2,
                          pamRespectsDendro = F,
                          # detectCutHeight = 0.75,
                          minModuleSize = 30,
                          maxBlockSize = 8000,
                          
                          # == Module Adjustments ==
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,
                          
                          # == TOM == Archive the run results in TOM file (saves time)
                          saveTOMs = T,
                          saveTOMFileBase = "ER",
                          
                          # == Output Options
                          numericLabels = T,
                          verbose = 3)
cor <- temp_cor     # Return cor function to original namespace

write_rds(netwk, 'PL_wcgna_res.RDS')

# Convert labels to colors for plotting
mergedColors = labels2colors(netwk$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(
  netwk$dendrograms[[1]],
  mergedColors[netwk$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )

table(netwk$colors)

module_df <- data.frame(
  gene_id = names(netwk$colors),
  colors = labels2colors(netwk$colors)
)

module_df[1:5,]

write_delim(module_df,
            file = "gene_modules.txt",
            delim = "\t")

# Get Module Eigengenes per cluster
MEs0 <- moduleEigengenes(input_mat, mergedColors)$eigengenes

# Reorder modules so similar modules are next to each other
MEs0 <- orderMEs(MEs0)
module_order = names(MEs0) %>% gsub("ME","", .)

# Add sample names
MEs0$sampleID = row.names(MEs0)

# tidy & plot data
mME = MEs0 %>%
  pivot_longer(-sampleID) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)
  )

mME %>% ggplot(., aes(x=sampleID, y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-trait Relationships", y = "Modules", fill="corr")

# pick out a few modules of interest here
modules_of_interest = c("green", "pink", "tan")

# Pull out list of genes in that module
submod = module_df %>%
  subset(colors %in% modules_of_interest)

row.names(module_df) = module_df$gene_id

# Get normalized expression for those genes
expr_normalized[1:5,1:10]

subexpr = expr_normalized[submod$gene_id,]

submod_df = data.frame(subexpr) %>%
  mutate(
    gene_id = row.names(.)
  ) %>%
  pivot_longer(-gene_id) %>%
  mutate(
    module = module_df[gene_id,]$colors
  )

submod_df %>% ggplot(., aes(x=name, y=value, group=gene_id)) +
  geom_line(aes(color = module),
            alpha = 0.2) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90)
  ) +
  facet_grid(rows = vars(module)) +
  labs(x = "treatment",
       y = "normalized expression")

# network file can be generated for Cytoscape or as an edge/vertices file
genes_of_interest = module_df %>%
  subset(colors %in% modules_of_interest)

expr_of_interest = expr_normalized[genes_of_interest$gene_id,]
expr_of_interest[1:5,1:5]

TOM = TOMsimilarityFromExpr(t(expr_of_interest),
                            power = picked_power)
edge_list = data.frame(TOM) %>%
  mutate(
    gene1 = row.names(.)
  ) %>%
  pivot_longer(-gene1) %>%
  dplyr::rename(gene2 = name, correlation = value) %>%
  unique() %>%
  subset(!(gene1==gene2)) %>%
  mutate(
    module1 = module_df[gene1,]$colors,
    module2 = module_df[gene2,]$colors
  )
# Export Network file to be read into Cytoscape, VisANT, etc
write_delim(edge_list,
            file = "edgelist.tsv",
            delim = "\t")
head(edge_list)
# GWENA workflow ----
library(GWENA)

CB_gene_counts <- read.csv("EDA_CB/CB_gene_counts.csv")

names(CB_gene_counts)[1] = "GeneId"
de_input = as.matrix(CB_gene_counts[,-1])
row.names(de_input) = CB_gene_counts$GeneId
# The DESeq2 functions need the values to be converted to integers
df <- round(de_input) %>%
  # The next steps require a data frame and round() returns a matrix
  as.data.frame() %>%
  # Only keep rows that have total counts above the cutoff
  dplyr::filter(rowSums(.) >= 50)

# Create a `DESeqDataSet` object & transform it
dds <- DESeqDataSetFromMatrix(
  countData = df, # Our prepped data frame with counts
  colData = Meta_CB, # Data frame with annotation for our samples
  design = ~1 # Here we are not specifying a model
)
dds_norm <- vst(dds)

# Retrieve the normalized count from the `DESeqDataSet` & Transpose it
normalized_counts <- assay(dds_norm) %>%
  t()
is_data_expr(normalized_counts)

# Network building
net <- build_net(normalized_counts, cor_func = "spearman", # Spearman cor-less sensitive to outliers
                 n_threads = 12)
net$metadata$power

# Fit of the power law to data ($R^2$) :
fit_power_table <- net$metadata$fit_power_table
fit_power_table[fit_power_table$Power == net$metadata$power, "SFT.R.sq"]

# detect module using hierarchical clustering
modules <- detect_modules(normalized_counts, 
                          net$network, 
                          detailled_result = TRUE,
                          merge_threshold = 0.25)
# Number of modules before merging :
length(unique(modules$modules_premerge))

# Number of modules after merging: 
length(unique(modules$modules))


layout_mod_merge <- plot_modules_merge(
  modules_premerge = modules$modules_premerge, 
  modules_merged = modules$modules)

# Resulting modules contain more genes whose repartition can be seen
ggplot(data.frame(modules$modules %>% stack), aes(x = ind)) + stat_count() +
  ylab("Number of genes") + xlab("Module")

# separate the positive (+ facet) and negative (- facet) correlations profile
plot_expression_profiles(normalized_counts, modules$modules)

# Functional enrichment
enrichment <- bio_enrich(modules$modules)
plot_enrichment(enrichment)

# Phenotype enrichment [data.frame]
allTraits <- Meta_CB[, c("cordBlood_ID","newborn_weight_kg","newborn_muac_cm", "father_age","mother_age")] #pulling out only continuous traits  
allTraits <- as.data.frame(allTraits)

phenotype_association <- associate_phenotype(
  modules$modules_eigengenes, 
  allTraits)

plot_modules_phenotype(phenotype_association)

module_example <- modules$modules$`2`
graph <- build_graph_from_sq_mat(net$network[module_example, module_example])

layout_mod_2 <- plot_module(graph, upper_weight_th = 0.999995, 
                            vertex.label.cex = 0, 
                            node_scaling_max = 7, 
                            legend_cex = 1)

net_mod_2 <- net$network[modules$modules$`2`, modules$modules$`2`] 
sub_clusters <- get_sub_clusters(net_mod_2)

layout_mod_2_sub_clust <- plot_module(graph, upper_weight_th = 0.999995,
                                      groups = sub_clusters,
                                      vertex.label.cex = 0, 
                                      node_scaling_max = 7, 
                                      legend_cex = 1)

# Expression by condition with data.frame
samples_by_cond <- lapply(Meta_CB$bw %>% unique, function(cond){
  df <- Meta_CB %>% 
    dplyr::filter(bw == "low") %>%
    dplyr::select(newborn_muac_cm, father_age)
  apply(df, 1, paste, collapse = "_")
}) %>% setNames(Meta_CB$bw %>% unique)

expr_by_cond <- lapply(samples_by_cond %>% names, function(low){
  samples <- samples_by_cond[[low]]
  normalized_counts[which(rownames(normalized_counts) %in% samples),]
}) %>% setNames(samples_by_cond %>% names)

# Network building and modules detection by condition
net_by_cond <- lapply(expr_by_cond, build_net, cor_func = "spearman", 
                      n_threads = 16, keep_matrices = "both")

mod_by_cond <- mapply(detect_modules, expr_by_cond, 
                      lapply(net_by_cond, `[[`, "network"), 
                      MoreArgs = list(detailled_result = TRUE), 
                      SIMPLIFY = FALSE)


comparison <- compare_conditions(expr_by_cond, 
                                 lapply(net_by_cond, `[[`, "adja_mat"), 
                                 lapply(net_by_cond, `[[`, "cor_mat"),  
                                 lapply(mod_by_cond, `[[`, "modules"), 
                                 pvalue_th = 0.05)