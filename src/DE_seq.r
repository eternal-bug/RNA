rm(list=ls())

source ("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
biocLite("pheatmap")
biocLite("biomaRt")
biocLite("org.Rn.eg.db")
biocLite("clusterProfiler")
biocLite("factoextra")

library(DESeq2)
library(pheatmap)
library(biomaRt)
library(org.Rn.eg.db)
library(clusterProfiler)
library(factoextra)

# data market
mart <- useDataset("rnorvegicus_gene_ensembl", useMart("ENSEMBL_MART_ENSEMBL"))

# row data and contain summary line
# 1 to 5 line should be trim
raw.data <- read.csv("merge.csv", header=TRUE, row.names = 1)
countdata <- raw.data[-(1:5),]

# screen(discard all 0 count)
countdata <- countdata[rowSums(countdata) > 0,]

# phenotype data
coldata <- read.table("../phenotype/phenotype.csv", row.names = 1, header = TRUE, sep = "," )
countdata <- countdata[row.names(coldata)]

# build dds object
dds = DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design= ~ treatment)

# sample relation

# variance stabilizing transformation and not blind sample
vsdata <- vst(dds, blind=FALSE)

# ======== PCA ===========
plotPCA(vsdata, intgroup="treatment")

# samples distance
gene_data_transform <- assay(vsdata)
# transposition and calculate distance
sampleDists <- dist(t(gene_data_transform))

# ====== hierarchical clustering ======
res <- hcut(sampleDists, k = 2, stand = TRUE)
# Visualize
fviz_dend(res, 
          rect = TRUE,
          rect_border="cluster",
          rect_lty=2,
          lwd=0.5,
          rect_fill = T,
          cex = 0.5,
          color_labels_by_k=T,
          horiz=T)

# transform matrix
sampleDistMatrix <- as.matrix(sampleDists)
# heatmap
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists
)

# different sample
compare_samples <- function(args_list){
  # args
  samples = args_list[["samples"]]
  levels = args_list[["levels"]]
  dir    = args_list[["dir"]]
  
  # data sampling
  coldata_g = subset(coldata, row.names(coldata) %in% samples)
  countdata_g = countdata[samples]
  dds = DESeqDataSetFromMatrix(countData = countdata_g, colData = coldata_g, design= ~ condition)
  dds$condition <- factor(as.vector(dds$condition), levels = levels)
  dds <- DESeq(dds)
  # result
  result <- results(dds, pAdjustMethod = "fdr", alpha = 0.05)
  result_order <- result[order(result$pvalue),]
  
  # print result info and 
  print("====> the order result summary.")
  head(result_order)
  summary(result_order)
  table(result_order$padj<0.05)
  
  # diff gene
  # the Padj adjusted by FDR, and set threshold value to 0.05
  # set 1 as threshold value of log2FoldChange and FC is 2/1(^) or 1/2(v)
  print("====> the different genes summary.")
  diff_gene <- subset(result_order, padj < 0.05 & abs(log2FoldChange) > 1)
  summary(diff_gene)
  # get gene id
  ensembl_gene_id <- row.names(diff_gene)
  
  # id transform
  rat_symbols <- getBM(attributes=c("ensembl_gene_id","external_gene_name","entrezgene_id", "description"),
                       filters = 'ensembl_gene_id', values = ensembl_gene_id, mart = mart)
  
  # merge symbols to DE analysis data
  diff_gene$ensembl_gene_id <- ensembl_gene_id
  # as dataframe
  diff_gene_dataframe <- as.data.frame(diff_gene)
  # merge
  diff_gene_symbols <- merge(diff_gene_dataframe, rat_symbols, by = c("ensembl_gene_id"))
  
  # write data into file
  print("====> write the table.")
  dir.create(dir, recursive = TRUE)
  write.table(result, paste(dir, "/all_gene.tsv", sep = ""), sep="\t", quote = FALSE)
  write.table(diff_gene_symbols, paste(dir, "/diff_gene.tsv", sep = ""), row.names = F,sep="\t", quote = FALSE)
  
  # enrichment analysis
  for(i in c("CC", "BP", "MF")){
    print(paste("====> begin to GO ", i, " analysis", sep = ""))
    ego <- enrichGO(gene          = rat_symbols$ensembl_gene_id,
                    OrgDb         = org.Rn.eg.db,
                    keyType       = "ENSEMBL",
                    ont           = paste(i, sep = ""),
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05)
    pdf(paste(dir, "/", i, ".pdf", sep = ""))
        print(dotplot(ego, showCategory = 30, title = paste("The GO ", i, " enrichment analysis", sep = "")))
    while (!is.null(dev.list()))  dev.off()
  }
}

# compare
compare_samples( list(samples=c("NH", "LHA1", "LHA2", "LHA3"), levels = c("HC","H1"), dir = "../stat/HA") )
compare_samples( list(samples=c("NH", "LHB1", "LHB2"        ), levels = c("HC","H2"), dir = "../stat/HB") )
compare_samples( list(samples=c("NH", "LHC1", "LHC2", "LHC3"), levels = c("HC","H3"), dir = "../stat/HC") )
compare_samples( list(samples=c("NM", "LMA1", "LMA2", "LMA3"), levels = c("MC","M1"), dir = "../stat/MA") )
compare_samples( list(samples=c("NM", "LMB1", "LMB2", "LMB3"), levels = c("MC","M2"), dir = "../stat/MB") )
compare_samples( list(samples=c("NM", "LMC1", "LMC2"        ), levels = c("MC","M3"), dir = "../stat/MC") )
