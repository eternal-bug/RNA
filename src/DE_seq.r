rm(list=ls())

source ("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
biocLite("pheatmap")
biocLite("biomaRt")
biocLite("org.Rn.eg.db")
biocLite("clusterProfiler")

library(DESeq2)
library(pheatmap)
library(biomaRt)
library(org.Rn.eg.db)
library(clusterProfiler)

mart <- useDataset("rnorvegicus_gene_ensembl", useMart("ENSEMBL_MART_ENSEMBL"))

# row data and contain summary line
# 1 to 5 line should be trim
raw.data <- read.csv("merge.csv", header=TRUE, row.names = 1)
countdata <- raw.data[-(1:5),]

# phenotype data
coldata <- read.table("../phenotype/phenotype.csv", row.names = 1, header = TRUE, sep = "," )
countdata <- countdata[row.names(coldata)]

# build dds object
dds = DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design= ~ treatment)

# sample relation

# variance stabilizing transformation and not blind sample
vsdata <- vst(dds, blind=FALSE)

# PCA
plotPCA(vsdata, intgroup="treatment")

# distance
gene_data_transform <- assay(vsdata)
# transposition and calculate distance
sampleDists <- dist(t(gene_data_transform))
# transform matrix
sampleDistMatrix <- as.matrix(sampleDists)
# heatmap
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists
)

