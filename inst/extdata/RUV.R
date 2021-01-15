#RUV with positive controls - lambda reads
library(RUVSeq)
#merge the counts
h3k27ac.filtered.counts <- assay(h3k27ac.filtered.data)
colnames(h3k27ac.filtered.counts) <- metadata.h3k27ac$samples
rownames(h3k27ac.filtered.counts) <- as.character(granges(h3k27ac.filtered.data))
lambda.rse <- readRDS("atac_lambda_reads.rds")
lambda.counts <- assay(lambda.rse)
colnames(lambda.counts) <- metadata.h3k27ac$samples
rownames(lambda.counts) <- as.character(granges(lambda.rse))
h3k27ac.filtered.data.for.ruvg <- rbind(h3k27ac.filtered.counts, lambda.counts)
#smarcc1
smarcc1.filtered.counts <- assay(smarcc1.filtered.data)
colnames(smarcc1.filtered.counts) <- metadata.smarcc1$samples
rownames(smarcc1.filtered.counts) <- as.character(granges(smarcc1.filtered.data))
lambda.rse <- readRDS("atac_lambda_reads.rds")
lambda.counts <- assay(lambda.rse)
colnames(lambda.counts) <- metadata.smarcc1$samples
rownames(lambda.counts) <- as.character(granges(lambda.rse))
smarcc1.filtered.data.for.ruvg <- rbind(smarcc1.filtered.counts, lambda.counts)
#RUVg
#try up to k=4
h3k27ac.ruvg.results.k1 <- RUVg(smarcc1.filtered.data.for.ruvg,
                          rownames(lambda.counts),
                          k = 1)
#k2
h3k27ac.ruvg.results.k2 <- RUVg(h3k27ac.filtered.data.for.ruvg,
                          rownames(lambda.counts),
                          k = 2)
#k3
h3k27ac.ruvg.results.k3 <- RUVg(h3k27ac.filtered.data.for.ruvg,
                          rownames(lambda.counts),
                          k = 3)
#k4
h3k27ac.ruvg.results.k4 <- RUVg(h3k27ac.filtered.data.for.ruvg,
                          rownames(lambda.counts),
                          k = 4)
