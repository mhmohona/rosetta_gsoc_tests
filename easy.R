.libPaths(c("~/R/libs", .libPaths()))
library(DESeq2)
library(airway)

data(airway)
dds <- DESeqDataSet(airway, design = ~ cell + dex)
dds <- DESeq(dds)
res <- results(dds)
resOrdered <- res[order(res$padj),]

print(head(resOrdered, 5))
