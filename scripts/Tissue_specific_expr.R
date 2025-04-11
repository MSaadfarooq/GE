# https://scienceparkstudygroup.github.io/master-advanced-forensics/index.html
library(tidyverse)
library(GGally)
library("ggrepel")
library(cluster)
library(pheatmap)

# load and set data ----
CB_gene_counts <- read_csv("CB/CB_gene_counts.csv")
PL_gene_counts <- read_csv("Outputs/PL_gene_counts.csv")

CB_gene_counts <- column_to_rownames(CB_gene_counts, var = "...1")
gene_counts <- column_to_rownames(gene_counts, var = "...1")
PL_gene_counts <- column_to_rownames(PL_gene_counts, var='geneID')

CB_gene_counts <- CB_gene_counts[1:9000,]
PL_gene_counts <- PL_gene_counts[1:9000,]
gene_counts <- gene_counts[1:9000,]

CB_gene_counts$CB0011 <- NULL
CB_gene_counts$CB0017 <- NULL

PL_gene_counts$PL0011 <- NULL
PL_gene_counts$PL0017 <- NULL

all.equal(names(gene_counts), names(CB_gene_counts))

df <- cbind(gene_counts, CB_gene_counts)
df <- cbind(df, PL_gene_counts)

# convert to long format
df_tidy <- df %>% 
  pivot_longer(cols = everything(), 
               names_to = "tissue", 
               values_to = "expression")
# summary stats
df_tidy %>% 
  summarise(max = max(expression), 
            min = min(expression), 
            average = mean(expression), 
            median = median(expression)
  )
# One tissue or sample ploting
df_tidy %>%
  filter(tissue == "FS0030") %>%
  filter(expression != 0) %>% # filter unexpressed
  ggplot(data = ., aes(x = tissue, y = log10(expression+1))) + # log10 tranformed
  geom_boxplot()

# boxplot for multiple samples of same tissues
df_tidy %>%
  filter(expression != 0) %>% 
  filter(grepl(pattern = "FS*", x = tissue)) %>% 
  ggplot(data = ., aes(x = tissue, y = log10(expression + 1), fill = tissue)) + 
  geom_boxplot(alpha = 0.1) + 
  theme(axis.text.x = element_text(angle = 90)) +
  guides(fill=FALSE)

# multiple tissues of same ID
df_tidy %>%
  filter(expression != 0) %>% 
  filter(grepl(pattern = "*80", x = tissue)) %>% # select only two tissues? can rename samples
  ggplot(data = ., aes(x = tissue, y = log10(expression + 1), fill = tissue)) + 
  geom_boxplot(alpha = 0.1) + 
  theme(axis.text.x = element_text(angle = 90)) +
  guides(fill=FALSE)


# relationships between the different tissues
df_tidy %>%
  filter(expression != 0) %>% 
  #filter(grepl(pattern = "CB*", x = tissue)) %>% 
  pivot_wider(names_from = tissue, values_from = expression) %>% 
  #dplyr::select(- gene_id) %>%  # as ggpairs only accept numerical values
  ggpairs(title = "Scatterplot matrix of tissues", upper = "blank") +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1))

df |> select(contains('FS')) |> ggpairs()

# PCA using custom function ----
df_pca <- t(df)

mypca <- function(x, center = TRUE, scale = TRUE){  
  # Samples should be in rows
  # Variables in the columns
  
  # remove columns containing only 0 values
  # not informative + cause svd() error
  x_without_zero_columns <- x[,colSums(x != 0) != 0] 
  
  # perform SVD
  SVD <- svd(scale(x_without_zero_columns, center = center, scale = scale))
  
  # create scores data frame
  scores <- as.data.frame(SVD$u %*% diag(SVD$d))
  rownames(scores) <- rownames(x)
  colnames(scores) <- paste0("PC", c(1:dim(scores)[2]))
  
  # create loadings data frams
  loadings <- data.frame(SVD$v)
  colnames(loadings) <- paste0("PC", c(1:dim(loadings)[2]))
  row.names(loadings) <- colnames(x_without_zero_columns)
  
  # create data frame for explained variances
  explained_var <- as.data.frame(round((SVD$d^2) / sum(SVD$d^2)*100, digits = 1))
  rownames(explained_var) <- paste0("PC", c(1:dim(loadings)[2]))
  colnames(explained_var) <- "exp_var"
  
  # return result
  return (list("scores" = scores, "loadings" = loadings, "explained_var" = explained_var))
}
pca <- mypca(x = df_pca, 
             center = TRUE, 
             scale = TRUE)
# Create a dataframe with all PC components (n = number of tissues)
exp_var_df <- data.frame(PC = seq(1:nrow(pca$explained_var)), exp_var = pca$explained_var)
exp_var_df <- data.frame(PC = 1:10, exp_var = pca$explained_var[1:10,]) # for top 10 PCs

# make the complete screeplot
ggplot(exp_var_df, aes(x = PC, y = exp_var, label = PC)) +
  ylab('explained variance (%)') + 
  ggtitle('explained variance per component (all PCs') + 
  geom_bar(stat = "identity") +
  labs(x = "Principal Component numbers") + 
  scale_y_continuous(limits = c(0, 50)) +
  scale_x_continuous(breaks = seq(from = 1, to = 10, by = 1))

# want to explain 80% of the variance, how many PCs would you need?
exp_var_df = mutate(exp_var_df, cum_var = cumsum(exp_var))

ggplot(exp_var_df, aes(x = PC, y = cum_var)) +
  geom_point() + 
  geom_line(group = 1) + 
  labs(x = "Principal Component", y = "Cumulative Explained Variance (%)") +
  scale_x_continuous(breaks = 1:10) +
  geom_hline(yintercept = 80, color = "red")

# this gives you all the PCs for which the cumulative explained variance is above 80%
exp_var_df[which(exp_var_df$cum_var > 80),]
# adding min() before this expression gives you the first row where it happens. 
exp_var_df[min(which(exp_var_df$cum_var > 80)),]

# Score plot: create dataframe of gene scores in new dimensional space created by the computed PCs
scores <- pca$scores
scores[1:5,1:5]

# extract explained variance to add to the axis labels
explained_var = pca$explained_var$exp_var

# useful for labelling the dots in the plot 
tissue_names = row.names(scores)

ggplot(scores, aes(x = PC1, y = PC2, label = tissue_names)) + 
  geom_point(aes(colour = factor(tissue_names))) + 
  xlab(paste0('PC1 (',explained_var[1],'%)')) + 
  ylab(paste0('PC2 (',explained_var[2],'%)')) + 
  ggtitle('PCA score plot') + 
  geom_text_repel() +
  guides(colour = FALSE) # remove the color legend (not informative)

loadings <- pca$loadings %>%
  rownames_to_column("gene_id") 

loadings[1:5,1:5]

top10genes_PC1_PC2 <- 
  loadings %>% 
  pivot_longer(cols = - gene_id, names_to = "PC", values_to = "loadings") %>% 
  filter(PC == "PC1" | PC == "PC2") %>%                                         
  group_by(PC) %>% 
  arrange(desc(abs(loadings))) %>%
  slice(1:10)     # add back gene symbol

loadings4plot <- inner_join(loadings, top10genes_PC1_PC2, by = "gene_id") %>% 
  select(gene_id, PC1, PC3)

ggplot(loadings4plot) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC3), 
               arrow = arrow(length = unit(0.1, "in")), colour = "brown") +
  geom_text_repel(data = loadings4plot, aes(x = PC1, y = PC3, label = gene_id), size = 2) + 
  labs(x = "PC1", y = "PC3")

# Hierarchical Clustering of the tissues
# conversion to matrix for distance calculation later
mat_expr <- df %>% as.matrix

min(mat_expr)
max(mat_expr)

# Scaling- this function scales columns rather than rows
mat_expr_scaled <- mat_expr %>% 
  t() %>%  # transpose the matrix before scaling
  scale(center = TRUE, scale = TRUE) %>% 
  t() %>%  # return it to its original format
  na.omit()

mean(mat_expr_scaled[1,]) # mean for ENSG00000000003  (close to 0)
sd(mat_expr_scaled[1,])   # standard deviation for ENSG00000000003 (unit variance of 1)

# compute a distance matrix between every tissue in a pairwise manner
distance_tissues <- dist(
  t(mat_expr_scaled),     # notice the t() to calculate the distances between tissues, not between genes
  method = "euclidean")
as.matrix(distance_tissues)[1:5,1:3]

# AGlomerative NESting method from the cluster package as it produces smaller clusters
# The AGNES clustering method coupled to Ward's cluster dissimilarity estimation method
hcl_tissue_ward <- cluster::agnes(x = distance_tissues,method = "ward")

# create a dendrogram from this hierarchical cluster object 
plot(hcl_tissue_ward, 
     hang  = -1) # this ensures that labels end up at the same level

# To select a threshold for “highly expressed genes” meaningfully, compute the percentiles of the TPM
df_expr_tidy <- df %>%
  rownames_to_column('gene_id') %>%
  pivot_longer(- gene_id, names_to = "tissue", values_to = "tpm")

# getting the percentiles
df_expr_tidy %>% 
  group_by(gene_id) %>% 
  summarise(median_tpm = median(tpm)) %>% 
  with(., round(x = quantile(x = median_tpm, probs = c(
    seq(from = 0,
        to = 0.9,
        by = 0.1), 
    seq(from = 0.9, to = 1, by = 0.01)),
    digits = 2)
  ) # output indicates 90% of the genes have a median TPM value inferior to 1254
  ) 
# keeping only the “very highly expressed genes” (> 99th percentile)
genes_selected <- df_expr_tidy %>% 
  group_by(gene_id) %>% 
  summarise(median_tpm = median(tpm)) %>% 
  ungroup() %>% 
  filter(median_tpm > 112) %>% 
  dplyr::pull(gene_id)

length(genes_selected)

# subset exp matrix using list of selected genes
genes_mat <- subset(x = mat_expr_scaled, 
                    subset = rownames(mat_expr_scaled) %in% genes_selected)
dim(genes_mat)

distance_genes = dist(x = genes_mat, method = "euclidean")

# The AGNES clustering method coupled to Ward's cluster dissimilarity estimation method
hcl_genes_ward <- cluster::agnes(x = distance_genes, method = "ward")

# You can already create a dendrogram from this hierarchical cluster object
plot(hcl_genes_ward,
     labels = FALSE,  # remove unreadable gene names
     which.plots = 2,
     main = "Gene hierarchical clustering (AGNES, Ward method)")

# Tissue-specific genes through feature engineering ----
df_expr_tidy <- df %>%
  rownames_to_column('gene_id') %>%
  pivot_longer(- gene_id, names_to = "tissue", values_to = "tpm")

threshold <- df_expr_tidy %>% 
  group_by(gene_id) %>% 
  summarise(median_tpm = median(tpm)) %>% 
  with(., round(x = quantile(x = median_tpm, probs = 0.9)))

genes_selected <- df_expr_tidy %>% 
  group_by(gene_id) %>% 
  summarise(median_tpm = median(tpm)) %>% 
  ungroup() %>% 
  filter(median_tpm > threshold) %>% 
  dplyr::pull(gene_id)

df_expr_tidy_filtered <- filter(df_expr_tidy, gene_id %in% genes_selected)

# extract gene TPM value in "CB"  
fs_gene_expr <- df_expr_tidy_filtered %>% 
  filter(grepl(pattern = "FS*", x = tissue)) %>% 
  select(gene_id, tpm) %>% 
  rename(fs_tpm = tpm)

# calculate median TPM value in all other tissues
other_tissues_expr <- 
  df_expr_tidy_filtered %>% 
  filter(!grepl(pattern = "FS*", x = tissue)) %>% 
  group_by(gene_id) %>% 
  summarise(other_tissues_median_tpm = median(tpm))

## Merge the two dataframes & Calculate fold change
fs_vs_other_tissues <- inner_join(x = fs_gene_expr, 
                                       y = other_tissues_expr, 
                                       by = "gene_id") %>% 
  mutate(fc = fs_tpm / other_tissues_median_tpm) %>% 
  mutate(log2_fc = log2(fc)) 

# normal distribution of log2FC?
ggplot(fs_vs_other_tissues, aes(x = log2_fc)) +
  geom_density()

# Remove -Inf or NaN values before calculating the mean
# calculate Z-score: Z-test analysis to compute the probability
mean_of_log2fc <- with(data = fs_vs_other_tissues, mean(log2_fc[!is.infinite(fs_vs_other_tissues$log2_fc)]))
sd_of_log2fc <- with(data = fs_vs_other_tissues, sd(log2_fc[!is.infinite(fs_vs_other_tissues$log2_fc)]))

# test for normality
ks.test(x = fs_vs_other_tissues$log2_fc, y = "pnorm", mean = mean_of_log2fc, sd = sd_of_log2fc)

fs_vs_other_tissues$zscore <- map_dbl(
  fs_vs_other_tissues$log2_fc, 
  function(x) (x - mean_of_log2fc) / sd_of_log2fc
)
# filter fc > 0 + calculate one-sided p-value
fs_specific_genes <- fs_vs_other_tissues %>%  
  filter(log2_fc > 0) %>%                   # FC superior to 1
  mutate(pval = 1 - pnorm(zscore)) %>%      # one-tailed p-value
  filter(pval < 0.01) %>% 
  arrange(desc(log2_fc))

fs_specific_genes

# Filter the original dataset to keep only the 68 semen-related genes
df_fs_expr <- df[rownames(df) %in% fs_specific_genes$gene_id, ]
View(df_expr)
# conversion to matrix and scale so that gene expression values become comparable
mat_expr <- df_expr %>% 
  column_to_rownames("gene_id") %>% 
  as.matrix()

mat_expr_scaled <- mat_expr %>% 
  t() %>% 
  scale(center = TRUE, scale = TRUE) %>% 
  t() %>% na.omit()

pheatmap(mat_expr_scaled)

# favourite genes's expression ----                      
# create an object to hold our file location
expr_genes_file <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE81nnn/GSE81932/suppl/GSE81932_Dataset02.txt.gz"
# show the first 10 lines of the file (using the file location object)
read_lines(expr_genes_file, n_max = 10)

# create a data object by reading the whole file contents in
periodically_expressed_genes <- read_lines(expr_genes_file)

ribi_annotation_file <- "ribosome_biogenesis_annotations.txt"
read_lines(ribi_annotation_file, n_max = 10)

# the ! is a comment character, and lines after that should be ignored
ribi_annotation <- read_tsv(ribi_annotation_file, comment = "!")

# pull out short informal gene names column and longer systematic names
ribi_annotation <- select(ribi_annotation, Gene = "Gene/Complex", SystematicName = "Systematic Name/Complex Accession")
# check number of distinct gene name & systematic name combos
n_distinct(ribi_annotation)
ribi_genes <- distinct(ribi_annotation)

# check which ribi genes are periodically expressed?
ribi_genes %>% filter(SystematicName %in% periodically_expressed_genes)

# Find all gene names and IDs
scer_names_estimates_file <- "scer-mrna-protein-absolute-estimate.txt"
# read in file uses the # symbol to represent comments
scer_names_estimates <- read_tsv(scer_names_estimates_file, comment = "#")
# create a new data object and rename the columns 
scer_gene_names <- select(scer_names_estimates, Gene = gene, SystematicName = orf)

# arrange by `gene` in DESCENDING order and `mrna` in ASCENDING order.
arrange(scer_names_estimates, desc(gene), mrna)

# Is NOP56 on the list of periodically expressed genes?
filter(scer_gene_names, Gene=='NOP56')
filter(scer_gene_names, SystematicName %in% periodically_expressed_genes, Gene == "NOP56")

# download the mRNA data 
mRNA_data <- read_tsv("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE81nnn/GSE81932/suppl/GSE81932_Dataset01.txt.gz")

# rename ORF column to mentally (and programmatically!) compare and link datasets
names(mRNA_data)[1] <- "SystematicName"

# join short gene names onto mRNA data from dataset 01 on `SystematicName`
mRNA_data_named <- left_join(mRNA_data, scer_gene_names, by = "SystematicName")

# find out how NOP56 gene expression changes, filter its observations
filter(mRNA_data_named, Gene == "NOP56")

# How does this compare with other genes of interest? filter them out
filter(mRNA_data_named, Gene %in% c("ACT1", "NOP16", "NOP56"))

# Reshape/tidy the data
mrna_long <- pivot_longer(data = mRNA_data_named, cols = 2:9, names_to = "Vol", values_to = "log2_ratio" ) %>% 
  separate(Vol, into = "Vol_fL") %>% 
  mutate(Vol_fL = as.numeric(Vol_fL))

# filter and Plot fav genes
mrna_long %>% filter(Gene %in% c("ACT1","NOP16","NOP56")) %>% 
  ggplot() + geom_line(aes(Vol_fL, log2_ratio, colour = Gene))

mrna_long %>% filter(Gene %in% c("ACT1","NOP16","NOP56")) %>% 
  ggplot() + geom_tile(aes(Vol_fL, Gene, fill = log2_ratio)) +
  scale_color_gradient2() + theme_minimal()