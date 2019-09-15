args <- commandArgs(trailingOnly = TRUE)

file <- args[1]
dir <- args[2]

if(is.na(file)){
    stop("lack the enter gene id file!")
}

if(is.na(dir)){
    dir = "./"
}

library(org.Rn.eg.db)
library(clusterProfiler)
library(biomaRt)
library(ggplot2)
library(stringr)

data <- read.table(file, header = F, sep="\t")

mart <- useDataset("rnorvegicus_gene_ensembl", useMart("ENSEMBL_MART_ENSEMBL"))

ensembl_gene_id <- data[1]

rat_symbols <- getBM(attributes=c("ensembl_gene_id","external_gene_name","entrezgene_id", "description"), filters = 'ensembl_gene_id', values = ensembl_gene_id, mart = mart)

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
