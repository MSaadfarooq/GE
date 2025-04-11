library(DESeq2)
library(MLSeq)

# may need large sample size
PL_counts <- read.csv("Outputs/PL_gene_counts.csv", row.names = 1)

meta_pl <- read.csv("PL/meta_pl.csv")

# selected top 100 features having the highest gene-wise variances to decrease computation
vars <- sort(apply(PL_counts, 1, var, na.rm = TRUE), decreasing = TRUE)
data <- PL_counts[names(vars)[1:100], ]
nTest <- ceiling(ncol(data) * 0.3)
ind <- sample(ncol(data), nTest, FALSE)

# Minimum count is set to 1 in order to prevent 0 division problem within classification models
data.train <- as.matrix(data[ ,-ind] + 1)
data.test <- as.matrix(data[ ,ind] + 1)
classtr <- DataFrame(condition = meta_pl[-ind, ])
classts <- DataFrame(condition = meta_pl[ind, ])

# We do not perform a differential expression analysis to select DEGs.
# However, in practice, DE analysis might be performed before fitting classifiers. 
data.trainS4 = DESeqDataSetFromMatrix(countData = data.train, colData = classtr,
                                      design = formula(~condition.bw))
data.testS4 = DESeqDataSetFromMatrix(countData = data.test, colData = classts,
                                     design = formula(~condition.bw))

# Support Vector Machines with Radial Kernel
fit <- classify(data = data.trainS4, method = "svmRadial",
                preProcessing = "deseq-rlog", ref = "T",
                control = trainControl(method = "repeatedcv", number = 2,
                                       repeats = 2, classProbs = TRUE))
show(fit)

