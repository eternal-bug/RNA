rm(list=ls())

# the R version below 3.5 can use this function to install bio package
# otherwise use BiocManager 
# library(BiocManager)
# options(BioC_mirror='https://mirrors.tuna.tsinghua.edu.cn/bioconductor')
# BiocManager::install("package")


# This is some source at home
# 中科大: http://mirrors.ustc.edu.cn/bioc/
# 清华: https://mirrors.tuna.tsinghua.edu.cn/CRAN/

# ======== check and install package =======
source ("https://bioconductor.org/biocLite.R")
# set the download source at home
options(BioC_mirror="http://mirrors.ustc.edu.cn/bioc/") #中科大
options(repos=structure(c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")))  #清华

package_list = c("DESeq2",
                 "pheatmap",
                 "biomaRt",
                 "org.Rn.eg.db",
                 "clusterProfiler",
                 "factoextra",
                 "stringr",
                 "parallel")
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    biocLite(p)
  }
  suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
}

# ======= preset data market =======
mart <- useDataset("rnorvegicus_gene_ensembl", useMart("ENSEMBL_MART_ENSEMBL"))


# ======= read raw data ========
# row data and contain summary line
# 1 to 5 line should be trim
raw.data <- read.csv("merge.csv", header=TRUE, row.names = 1)
countdata <- raw.data[-(1:5),]

# === screen(discard all 0 count) ===
countdata <- countdata[rowSums(countdata) > 0,]

# phenotype data
coldata <- read.table("../phenotype/phenotype.csv", row.names = 1, header = TRUE, sep = "," )
countdata <- countdata[row.names(coldata)]

# ======= build dds object ========
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
# ============ heatmap ===========
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists
)

# ======== samples compare =======

compare_samples <- function(args_list){
  # args
  samples = args_list[["samples"]]
  levels = args_list[["levels"]]
  dir    = args_list[["dir"]]
  
  # data sampling
  coldata_g = subset(coldata, row.names(coldata) %in% samples)
  countdata_g = countdata[samples]
  dds = DESeqDataSetFromMatrix(countData = countdata_g, colData = coldata_g, design= ~ treatment)
  
  dds$treatment = factor(as.vector(dds$treatment), levels = levels)
  
  dds = DESeq(dds)
  
  # result
  result = results(dds, pAdjustMethod = "fdr", alpha = 0.05)
  result_order = result[order(result$pvalue),]
  
  # print result info and 
  print("====> the order result summary.")
  head(result_order)
  summary(result_order)
  table(result_order$padj<0.05)
  
  # diff gene
  # the Padj adjusted by FDR, and set threshold value to 0.05
  # set 1 as threshold value of log2FoldChange and FC is 2/1(^) or 1/2(v)
  print("====> the different genes summary.")
  diff_gene = subset(result_order, padj < 0.05 & abs(log2FoldChange) > 1)
  summary(diff_gene)
  # get gene id
  ensembl_gene_id = row.names(diff_gene)
  
  # id transform
  rat_symbols = getBM(attributes=c("ensembl_gene_id","external_gene_name","entrezgene_id", "description"),
                       filters = 'ensembl_gene_id', values = ensembl_gene_id, mart = mart)
  
  # merge symbols to DE analysis data
  diff_gene$ensembl_gene_id = ensembl_gene_id
  # as dataframe
  diff_gene_dataframe = as.data.frame(diff_gene)
  # merge
  diff_gene_symbols = merge(diff_gene_dataframe, rat_symbols, by = c("ensembl_gene_id"))
  
  # write data into file
  print("====> write the table.")
  dir.create(dir, recursive = TRUE)
  write.table(result, paste(dir, "/all_gene.tsv", sep = ""), sep="\t", quote = FALSE)
  write.table(diff_gene_symbols, paste(dir, "/diff_gene.tsv", sep = ""), row.names = F,sep="\t", quote = FALSE)
  
  # enrichment analysis
  # ============= GO ====================
  for(i in c("CC", "BP", "MF")){
    print(paste("====> begin to GO ", i, " analysis", sep = ""))
    ego = enrichGO(gene          = rat_symbols$ensembl_gene_id,
                    OrgDb         = org.Rn.eg.db,
                    keyType       = "ENSEMBL",
                    ont           = paste(i, sep = ""),
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05)
    # bubble pic
    pdf(paste(dir, "/", i, "_bubble", ".pdf", sep = ""))
        print(
          dotplot(ego, showCategory = 30,
                  title = paste("The GO ", i, " enrichment analysis", sep = "")) +
          scale_y_discrete(labels = str_wrap(ego@result$Description, width = 60)) + 
          theme(axis.text.y = element_text(size = 8))
        )
    while (!is.null(dev.list())){
        dev.off()
    }
    # graph pic
    pdf(paste(dir, "/", i, "_graph", ".pdf", sep = ""))
        print(plotGOgraph(ego))
    while (!is.null(dev.list())){
      dev.off()
    }
    # store data
    write.table(ego@result, paste(dir, "/", i, ".tsv", sep = ""),
                row.names = F,sep="\t", quote = FALSE)
  }
  
  # =========== KEGG =================
  kk = enrichKEGG(gene = rat_symbols$entrezgene_id, 
                   organism ='rno',
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   minGSSize = 1,
                   #readable = TRUE ,
                   use_internal_data = FALSE)
  # bubble pic
  pdf(paste(dir, "/", "kegg_graph", ".pdf", sep = ""))
  print(
    dotplot(kk, showCategory = 30,
            title = "The KEGG enrichment analysis") +
      scale_y_discrete(labels = str_wrap(kk@result$Description, width = 60)) + 
      theme(axis.text.y = element_text(size = 8))
  )
  while (!is.null(dev.list())){
    dev.off()
  }
  # store data
  write.table(kk@result, paste(dir, "/","kegg.tsv", sep = ""),
              row.names = F,sep="\t", quote = FALSE)
}


compare_samples( list(samples=c("NH", "LHA1", "LHA2", "LHA3"), levels = c("H_contr","H_exper"), dir = "../stat/HA") )
compare_samples( list(samples=c("NH", "LHB1", "LHB2"        ), levels = c("H_contr","H_exper"), dir = "../stat/HB") )
compare_samples( list(samples=c("NH", "LHC1", "LHC2", "LHC3"), levels = c("H_contr","H_exper"), dir = "../stat/HC") )
compare_samples( list(samples=c("NM", "LMA1", "LMA2", "LMA3"), levels = c("M_contr","M_exper"), dir = "../stat/MA") )
compare_samples( list(samples=c("NM", "LMB1", "LMB2", "LMB3"), levels = c("M_contr","M_exper"), dir = "../stat/MB") )
compare_samples( list(samples=c("NM", "LMC1", "LMC2"        ), levels = c("M_contr","M_exper"), dir = "../stat/MC") )


compare_samples( list(samples=c("NM", "LMB1", "LMC1"        ), levels = c("M_contr","M_exper"), dir = "../stat/MB1_MC1") )
compare_samples( list(samples=c("NM", "LMB2", "LMB3", "LMC2"), levels = c("M_contr","M_exper"), dir = "../stat/MB2_3_MC2") )
