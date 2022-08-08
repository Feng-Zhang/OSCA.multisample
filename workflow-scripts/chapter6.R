# set path ------------------------------------------
rm(list=ls());options(stringsAsFactors=FALSE)
project_dir <- rprojroot::find_rstudio_root_file()
raw_dir <- file.path(project_dir,"raw-data")
process_dir <- file.path(project_dir,"process-data")
output_dir <- file.path(project_dir,"output-results")
temp_dir <- file.path(project_dir,"temp-data")
bin_dir <- file.path(project_dir,"bin")
script_dir <- file.path(project_dir,"workflow-scripts")
work_dir <- output_dir
setwd(work_dir)

# library   -----------------------------------------
library(tidyverse)

# set default parameter  ---------------------------
prefix <- "6-"
output_prefix <- file.path(work_dir,prefix)
process_prefix <-  file.path(process_dir,prefix)

# perpare data ------------------------------
#--- loading ---#
library(MouseGastrulationData)
sce.chimera <- WTChimeraData(samples=5:10)
sce.chimera

#--- feature-annotation ---#
library(scater)
rownames(sce.chimera) <- uniquifyFeatureNames(
  rowData(sce.chimera)$ENSEMBL, rowData(sce.chimera)$SYMBOL)

#--- quality-control ---#
drop <- sce.chimera$celltype.mapped %in% c("stripped", "Doublet")
sce.chimera <- sce.chimera[,!drop]

#--- normalization ---#
sce.chimera <- logNormCounts(sce.chimera)

#--- variance-modelling ---#
library(scran)
dec.chimera <- modelGeneVar(sce.chimera, block=sce.chimera$sample)
chosen.hvgs <- dec.chimera$bio > 0

#--- merging ---#
library(batchelor)
set.seed(01001001)
merged <- correctExperiments(sce.chimera, 
                             batch=sce.chimera$sample, 
                             subset.row=chosen.hvgs,
                             PARAM=FastMnnParam(
                               merge.order=list(
                                 list(1,3,5), # WT (3 replicates)
                                 list(2,4,6)  # td-Tomato (3 replicates)
                               )
                             )
)

#--- clustering ---#
g <- buildSNNGraph(merged, use.dimred="corrected")
clusters <- igraph::cluster_louvain(g)
colLabels(merged) <- factor(clusters$membership)

#--- dimensionality-reduction ---#
merged <- runTSNE(merged, dimred="corrected", external_neighbors=TRUE)
merged <- runUMAP(merged, dimred="corrected", external_neighbors=TRUE)
save(merged,file=paste0(output_prefix,"merged.RData"))

# analysis --------------------
abundances <- table(merged$celltype.mapped, merged$sample) 
abundances <- unclass(abundances) 
head(abundances)

library(edgeR)
# Attaching some column metadata.
extra.info <- colData(merged)[match(colnames(abundances), merged$sample),]
y.ab <- DGEList(abundances, samples=extra.info)
y.ab

keep <- filterByExpr(y.ab, group=y.ab$samples$tomato)
y.ab <- y.ab[keep,]
summary(keep)

design <- model.matrix(~factor(pool) + factor(tomato), y.ab$samples)
y.ab <- estimateDisp(y.ab, design, trend="none")
summary(y.ab$common.dispersion)
fit.ab <- glmQLFit(y.ab, design, robust=TRUE, abundance.trend=FALSE)
summary(fit.ab$var.prior)
res <- glmQLFTest(fit.ab, coef=ncol(design))
summary(decideTests(res))
topTags(res)
