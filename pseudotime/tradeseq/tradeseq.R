## ----starttime---------------------------------------------------------------------------------------------------------------
# save start time for script
start_tm <- Sys.time()
start_tm


## ----setup, include=FALSE----------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo=TRUE, warning=TRUE, message=TRUE, cache=TRUE, cache.lazy = FALSE, dev=c('png','pdf'), fig.width=8, fig.height=8)



## ----make_outdir-------------------------------------------------------------------------------------------------------------
outdir <- "./tradeseq_out_files/"

dir.create(outdir, recursive=TRUE)


## ----loadlibs, echo=TRUE, warning=TRUE, message=TRUE, cache=FALSE------------------------------------------------------------
library(dplyr)
library(stringr)
library(readr)
library(slingshot)
library(tradeSeq)
library(scater)



## ----readsce-----------------------------------------------------------------------------------------------------------------
sce <- readRDS("../pseudotime_out_files/slingshot_sce.rds")


## ----tradeseq----------------------------------------------------------------------------------------------------------------
BPPARAM <- BiocParallel::bpparam()
BPPARAM$workers <- 8

#keepGenes <- getTopHVGs(sce, n=1000, row.names=FALSE)

keep_feature <- nexprs(
  sce, 
  byrow = TRUE, 
  detection_limit = 3 # 1 less than min recommended here at https://github.com/statOmics/tradeSeq/issues/129
) >= 100
table(keep_feature)
keepGenes <- which(keep_feature)

batch <- factor(colData(sce)$Sample)
U <- model.matrix(~batch)

# fit negative binomial GAM
sce <- fitGAM(sce, U=U, genes=keepGenes, parallel=TRUE, BPPARAM=BPPARAM)
write_rds(sce, paste0(outdir, "tradeseq_sce.rds"))

# test for dynamic expression
ATres <- associationTest(sce, lineages = TRUE)#, l2fc = log2(2))
write_rds(ATres, paste0(outdir, "ATres.rds"))

# topgenes <- rownames(ATres[order(ATres$pvalue), ])[1:250]
# pst.ord <- order(sce$slingPseudotime_1, na.last = NA)
# heatdata <- assays(sce)$logcounts[topgenes, pst.ord]
# #heatclus <- sce$GMM[pst.ord]
# 
# heatmap(heatdata, Colv = NA)



## ----sessioninfo-------------------------------------------------------------------------------------------------------------
sessionInfo()


## ----endtime-----------------------------------------------------------------------------------------------------------------
# output time taken to run script
end_tm <- Sys.time()
end_tm
end_tm - start_tm


