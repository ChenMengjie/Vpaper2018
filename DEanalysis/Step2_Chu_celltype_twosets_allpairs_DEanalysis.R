############################################################### 
############# running DE analysis on raw count data ###########
###############################################################  
i <- 1 #1:10
j <- 1 #1:21

load("GSE75748_sc_cell_type.rda")
xx <- gene.expression
aa <- gsub("\\_[^<>]*", "", colnames(xx))
celltype <- unique(aa)

pairs <- matrix(c("EC", "DEC", "HFF", "DEC", "NPC", "DEC", "TB", "DEC",
"H9", "EC", "TB", "H1", "H9", "HFF", "H9", "NPC", "H9", "TB", "H9", "H1",  
"EC", "H1", "EC", "TB", "EC", "NPC",  "HFF", "TB", "HFF", "NPC",
"TB",  "NPC", "H1", "DEC",  "H9", "DEC",  "EC", "HFF",  "NPC", "TB",
"H1", "NPC"), byrow = T, ncol = 2)

library(VIPER)

load(paste0("Chu_cell_type_twosets_copy", i))

type1 <- aa[set1]
type2 <- aa[set2]

flag1.1 <- which(type1%in%pairs[j, 1])
flag1.2 <- which(type1%in%pairs[j, 2])
selected.x1 <- cbind(gene.expression1[, flag1.1], gene.expression1[, flag1.2])
X1 <- c(rep(0, length(flag1.1)), rep(1, length(flag1.2)))

flag2.1 <- which(type2%in%pairs[j, 1])
flag2.2 <- which(type2%in%pairs[j, 2])
selected.x2 <- cbind(gene.expression2[, flag2.1], gene.expression2[, flag2.2])
X2 <- c(rep(0, length(flag2.1)), rep(1, length(flag2.2)))



edgeR.LRT.res1 <- edgeR_others(selected.x1, matrix(X1, length(X1), 1), method = "LRT")
edgeR.LRT.res2 <- edgeR_others(selected.x2, matrix(X2, length(X2), 1), method = "LRT")
save(edgeR.LRT.res1, edgeR.LRT.res2, file = paste0("cell_type_copy", i, "_pairs",  j, "_expression_summary_edgeR_LRT"))


edgeR.QLF.res1 <- edgeR_others(selected.x1, matrix(X1, length(X1), 1), method = "QLF")
edgeR.QLF.res2 <- edgeR_others(selected.x2, matrix(X2, length(X2), 1), method = "QLF")
save(edgeR.QLF.res1, edgeR.QLF.res2, file = paste0("cell_type_copy", i, "_pairs",  j, "_expression_summary_edgeR_QLF"))




library(scde)
sg <- factor(X1)
names(sg) <- colnames(selected.x1)
cd <- clean.counts(selected.x1, min.lib.size=100, min.reads = 1, min.detected = 1) # which is a gene by sample count matrix
counts <- apply(cd, 2, function(x) {storage.mode(x) <- 'integer'; x})
o.ifm <- scde.error.models(counts = counts, groups = sg, n.cores = 1, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1) 
valid.cells <- o.ifm$corr.a > 0
o.ifm <- o.ifm[valid.cells, ]
o.prior <- scde.expression.prior(models = o.ifm, counts = counts, length.out = 400, show.plot = FALSE)
groups <- factor(X1)
names(groups) <- row.names(o.ifm)
scde.res1 <- scde.expression.difference(o.ifm, cd, o.prior, groups = groups, n.randomizations = 100, n.cores = 1, verbose = 1)

sg <- factor(X2)
names(sg) <- colnames(selected.x2)
cd <- clean.counts(selected.x2, min.lib.size=100, min.reads = 1, min.detected = 1) # which is a gene by sample count matrix
counts <- apply(cd, 2, function(x) {storage.mode(x) <- 'integer'; x})
o.ifm <- scde.error.models(counts = counts, groups = sg, n.cores = 1, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1) 
valid.cells <- o.ifm$corr.a > 0
o.ifm <- o.ifm[valid.cells, ]
o.prior <- scde.expression.prior(models = o.ifm, counts = counts, length.out = 400, show.plot = FALSE)
groups <- factor(X2)
names(groups) <- row.names(o.ifm)
scde.res2 <- scde.expression.difference(o.ifm, cd, o.prior, groups = groups, n.randomizations = 100, n.cores = 1, verbose = 1)

save(scde.res1, scde.res2, file = paste0("cell_type_copy", i, "_pairs",  j, "_expression_summary_scde"))


###################################################################### 
############# running DE analysis on drImpute imputed data ###########
######################################################################  

i <- 1 #1:10
j <- 1 #1:21


load("GSE75748_sc_cell_type.rda")
xx <- gene.expression
aa <- gsub("\\_[^<>]*", "", colnames(xx))
celltype <- unique(aa)


pairs <- matrix(c("EC", "DEC", "HFF", "DEC", "NPC", "DEC", "TB", "DEC",
"H9", "EC", "TB", "H1", "H9", "HFF", "H9", "NPC", "H9", "TB", "H9", "H1",  
"EC", "H1", "EC", "TB", "EC", "NPC",  "HFF", "TB", "HFF", "NPC",
"TB",  "NPC", "H1", "DEC",  "H9", "DEC",  "EC", "HFF",  "NPC", "TB",
"H1", "NPC"), byrow = T, ncol = 2)



load(paste0("Chu_cell_type_twosets_copy", i))
type1 <- aa[set1]
type2 <- aa[set2]

cellnames1 <- colnames(gene.expression1)
cellnames2 <- colnames(gene.expression2)

library(VIPER)

load(paste0("Chu_drImpute_cell_type_twosets_tc_copy", i))
		
gene.expression1 <- apply(lnpX_imp1, 2, function(x){
	exp(x) - 1
})

gene.expression2 <- apply(lnpX_imp2, 2, function(x){
	exp(x) - 1
})

colnames(gene.expression1) <- cellnames1
colnames(gene.expression2) <- cellnames2

flag1.1 <- which(type1%in%pairs[j, 1])
flag1.2 <- which(type1%in%pairs[j, 2])
selected.x1 <- cbind(gene.expression1[, flag1.1], gene.expression1[, flag1.2])
X1 <- c(rep(0, length(flag1.1)), rep(1, length(flag1.2)))

flag2.1 <- which(type2%in%pairs[j, 1])
flag2.2 <- which(type2%in%pairs[j, 2])
selected.x2 <- cbind(gene.expression2[, flag2.1], gene.expression2[, flag2.2])
X2 <- c(rep(0, length(flag2.1)), rep(1, length(flag2.2)))


DESeq2.res1 <- DESeq2(selected.x1, matrix(X1, length(X1), 1))
DESeq2.res2 <- DESeq2(selected.x2, matrix(X2, length(X2), 1))
save(DESeq2.res1, DESeq2.res2, file = paste0("cell_type_copy", i, "_pairs",  j, "_drImpute_summary_DESeq2"))


edgeR.LRT.res1 <- edgeR_others(selected.x1, matrix(X1, length(X1), 1), method = "LRT")
edgeR.LRT.res2 <- edgeR_others(selected.x2, matrix(X2, length(X2), 1), method = "LRT")
save(edgeR.LRT.res1, edgeR.LRT.res2, file = paste0("cell_type_copy", i, "_pairs",  j, "_drImpute_summary_edgeR_LRT"))

edgeR.QLF.res1 <- edgeR_others(selected.x1, matrix(X1, length(X1), 1), method = "QLF")
edgeR.QLF.res2 <- edgeR_others(selected.x2, matrix(X2, length(X2), 1), method = "QLF")
save(edgeR.QLF.res1, edgeR.QLF.res2, file = paste0("cell_type_copy", i, "_pairs",  j, "_drImpute_summary_edgeR_QLF"))


library(scde)
sg <- factor(X1)
names(sg) <- colnames(selected.x1)
cd <- clean.counts(selected.x1, min.lib.size=100, min.reads = 1, min.detected = 1) # which is a gene by sample count matrix
counts <- apply(cd, 2, function(x) {storage.mode(x) <- 'integer'; x})
o.ifm <- scde.error.models(counts = counts, groups = sg, n.cores = 1, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1) 
valid.cells <- o.ifm$corr.a > 0
o.ifm <- o.ifm[valid.cells, ]
o.prior <- scde.expression.prior(models = o.ifm, counts = counts, length.out = 400, show.plot = FALSE)
groups <- factor(X1)
names(groups) <- row.names(o.ifm)
scde.res1 <- scde.expression.difference(o.ifm, cd, o.prior, groups = groups, n.randomizations = 100, n.cores = 1, verbose = 1)

sg <- factor(X2)
names(sg) <- colnames(selected.x2)
cd <- clean.counts(selected.x2, min.lib.size=100, min.reads = 1, min.detected = 1) # which is a gene by sample count matrix
counts <- apply(cd, 2, function(x) {storage.mode(x) <- 'integer'; x})
o.ifm <- scde.error.models(counts = counts, groups = sg, n.cores = 1, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1) 
valid.cells <- o.ifm$corr.a > 0
o.ifm <- o.ifm[valid.cells, ]
o.prior <- scde.expression.prior(models = o.ifm, counts = counts, length.out = 400, show.plot = FALSE)
groups <- factor(X2)
names(groups) <- row.names(o.ifm)
scde.res2 <- scde.expression.difference(o.ifm, cd, o.prior, groups = groups, n.randomizations = 100, n.cores = 1, verbose = 1)

save(scde.res1, scde.res2, file = paste0("cell_type_copy", i, "_pairs",  j, "_drImpute_summary_scde"))



###################################################################
############# running DE analysis on SAVER imputed data ###########
###################################################################

i <- 1 #1:10
j <- 1 #1:21


load("GSE75748_sc_cell_type.rda")
xx <- gene.expression
aa <- gsub("\\_[^<>]*", "", colnames(xx))
celltype <- unique(aa)

pairs <- matrix(c("EC", "DEC", "HFF", "DEC", "NPC", "DEC", "TB", "DEC",
"H9", "EC", "TB", "H1", "H9", "HFF", "H9", "NPC", "H9", "TB", "H9", "H1",  
"EC", "H1", "EC", "TB", "EC", "NPC",  "HFF", "TB", "HFF", "NPC",
"TB",  "NPC", "H1", "DEC",  "H9", "DEC",  "EC", "HFF",  "NPC", "TB",
"H1", "NPC"), byrow = T, ncol = 2)


load(paste0("Chu_cell_type_twosets_copy", i))
type1 <- aa[set1]
type2 <- aa[set2]

library(VIPER)

load(paste0("Chu_SAVER_cell_type_twosets_tc_copy", i))
		
gene.expression1 <- saver_imp1$estimate
gene.expression2 <- saver_imp2$estimate

flag1.1 <- which(type1%in%pairs[j, 1])
flag1.2 <- which(type1%in%pairs[j, 2])
selected.x1 <- cbind(gene.expression1[, flag1.1], gene.expression1[, flag1.2])
X1 <- c(rep(0, length(flag1.1)), rep(1, length(flag1.2)))

flag2.1 <- which(type2%in%pairs[j, 1])
flag2.2 <- which(type2%in%pairs[j, 2])
selected.x2 <- cbind(gene.expression2[, flag2.1], gene.expression2[, flag2.2])
X2 <- c(rep(0, length(flag2.1)), rep(1, length(flag2.2)))


DESeq2.res1 <- DESeq2(selected.x1, matrix(X1, length(X1), 1))
DESeq2.res2 <- DESeq2(selected.x2, matrix(X2, length(X2), 1))
save(DESeq2.res1, DESeq2.res2, file = paste0("cell_type_copy", i, "_pairs",  j, "_SAVER_summary_DESeq2"))

edgeR.LRT.res1 <- edgeR_others(selected.x1, matrix(X1, length(X1), 1), method = "LRT")
edgeR.LRT.res2 <- edgeR_others(selected.x2, matrix(X2, length(X2), 1), method = "LRT")
save(edgeR.LRT.res1, edgeR.LRT.res2, file = paste0("cell_type_copy", i, "_pairs",  j, "_SAVER_summary_edgeR_LRT"))

edgeR.QLF.res1 <- edgeR_others(selected.x1, matrix(X1, length(X1), 1), method = "QLF")
edgeR.QLF.res2 <- edgeR_others(selected.x2, matrix(X2, length(X2), 1), method = "QLF")
save(edgeR.QLF.res1, edgeR.QLF.res2, file = paste0("cell_type_copy", i, "_pairs",  j, "_SAVER_summary_edgeR_QLF"))


library(scde)
sg <- factor(X1)
names(sg) <- colnames(selected.x1)
cd <- clean.counts(selected.x1, min.lib.size=100, min.reads = 1, min.detected = 1) # which is a gene by sample count matrix
counts <- apply(cd, 2, function(x) {storage.mode(x) <- 'integer'; x})
o.ifm <- scde.error.models(counts = counts, groups = sg, n.cores = 1, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1) 
valid.cells <- o.ifm$corr.a > 0
o.ifm <- o.ifm[valid.cells, ]
o.prior <- scde.expression.prior(models = o.ifm, counts = counts, length.out = 400, show.plot = FALSE)
groups <- factor(X1)
names(groups) <- row.names(o.ifm)
scde.res1 <- scde.expression.difference(o.ifm, cd, o.prior, groups = groups, n.randomizations = 100, n.cores = 1, verbose = 1)

sg <- factor(X2)
names(sg) <- colnames(selected.x2)
cd <- clean.counts(selected.x2, min.lib.size=100, min.reads = 1, min.detected = 1) # which is a gene by sample count matrix
counts <- apply(cd, 2, function(x) {storage.mode(x) <- 'integer'; x})
o.ifm <- scde.error.models(counts = counts, groups = sg, n.cores = 1, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1) 
valid.cells <- o.ifm$corr.a > 0
o.ifm <- o.ifm[valid.cells, ]
o.prior <- scde.expression.prior(models = o.ifm, counts = counts, length.out = 400, show.plot = FALSE)
groups <- factor(X2)
names(groups) <- row.names(o.ifm)
scde.res2 <- scde.expression.difference(o.ifm, cd, o.prior, groups = groups, n.randomizations = 100, n.cores = 1, verbose = 1)

save(scde.res1, scde.res2, file = paste0("cell_type_copy", i, "_pairs",  j, "_SAVER_summary_scde"))


###################################################################
############# running DE analysis on scImpute imputed data ########
###################################################################


load("GSE75748_sc_cell_type.rda")
xx <- gene.expression
aa <- gsub("\\_[^<>]*", "", colnames(xx))
celltype <- unique(aa)

pairs <- matrix(c("EC", "DEC", "HFF", "DEC", "NPC", "DEC", "TB", "DEC",
"H9", "EC", "TB", "H1", "H9", "HFF", "H9", "NPC", "H9", "TB", "H9", "H1",  
"EC", "H1", "EC", "TB", "EC", "NPC",  "HFF", "TB", "HFF", "NPC",
"TB",  "NPC", "H1", "DEC",  "H9", "DEC",  "EC", "HFF",  "NPC", "TB",
"H1", "NPC"), byrow = T, ncol = 2)



load(paste0("Chu_cell_type_twosets_copy", i))
type1 <- aa[set1]
type2 <- aa[set2]


library(VIPER)

scimpute.res1 <- read.csv(paste0("twosets_copy", i, "_1scimpute_count.csv"))
rownames(scimpute.res1) <- scimpute.res1[, 1]
gene.expression1 <- scimpute.res1[, -1]

scimpute.res2 <- read.csv(paste0("twosets_copy", i, "_2scimpute_count.csv"))
rownames(scimpute.res2) <- scimpute.res2[, 1]
gene.expression2 <- scimpute.res2[, -1]

flag1.1 <- which(type1%in%pairs[j, 1])
flag1.2 <- which(type1%in%pairs[j, 2])
selected.x1 <- cbind(gene.expression1[, flag1.1], gene.expression1[, flag1.2])
X1 <- c(rep(0, length(flag1.1)), rep(1, length(flag1.2)))

flag2.1 <- which(type2%in%pairs[j, 1])
flag2.2 <- which(type2%in%pairs[j, 2])
selected.x2 <- cbind(gene.expression2[, flag2.1], gene.expression2[, flag2.2])
X2 <- c(rep(0, length(flag2.1)), rep(1, length(flag2.2)))


DESeq2.res1 <- DESeq2(selected.x1, matrix(X1, length(X1), 1))
DESeq2.res2 <- DESeq2(selected.x2, matrix(X2, length(X2), 1))
save(DESeq2.res1, DESeq2.res2, file = paste0("cell_type_copy", i, "_pairs",  j, "_scImpute_summary_DESeq2"))

edgeR.LRT.res1 <- edgeR_others(selected.x1, matrix(X1, length(X1), 1), method = "LRT")
edgeR.LRT.res2 <- edgeR_others(selected.x2, matrix(X2, length(X2), 1), method = "LRT")
save(edgeR.LRT.res1, edgeR.LRT.res2, file = paste0("cell_type_copy", i, "_pairs",  j, "_scImpute_summary_edgeR_LRT"))

edgeR.QLF.res1 <- edgeR_others(selected.x1, matrix(X1, length(X1), 1), method = "QLF")
edgeR.QLF.res2 <- edgeR_others(selected.x2, matrix(X2, length(X2), 1), method = "QLF")
save(edgeR.QLF.res1, edgeR.QLF.res2, file = paste0("cell_type_copy", i, "_pairs",  j, "_scImpute_summary_edgeR_QLF"))

library(scde)
sg <- factor(X1)
names(sg) <- colnames(selected.x1)
cd <- clean.counts(selected.x1, min.lib.size=100, min.reads = 1, min.detected = 1) # which is a gene by sample count matrix
counts <- apply(cd, 2, function(x) {storage.mode(x) <- 'integer'; x})
o.ifm <- scde.error.models(counts = counts, groups = sg, n.cores = 1, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1) 
valid.cells <- o.ifm$corr.a > 0
o.ifm <- o.ifm[valid.cells, ]
o.prior <- scde.expression.prior(models = o.ifm, counts = counts, length.out = 400, show.plot = FALSE)
groups <- factor(X1)
names(groups) <- row.names(o.ifm)
scde.res1 <- scde.expression.difference(o.ifm, cd, o.prior, groups = groups, n.randomizations = 100, n.cores = 1, verbose = 1)

sg <- factor(X2)
names(sg) <- colnames(selected.x2)
cd <- clean.counts(selected.x2, min.lib.size=100, min.reads = 1, min.detected = 1) # which is a gene by sample count matrix
counts <- apply(cd, 2, function(x) {storage.mode(x) <- 'integer'; x})
o.ifm <- scde.error.models(counts = counts, groups = sg, n.cores = 1, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1) 
valid.cells <- o.ifm$corr.a > 0
o.ifm <- o.ifm[valid.cells, ]
o.prior <- scde.expression.prior(models = o.ifm, counts = counts, length.out = 400, show.plot = FALSE)
groups <- factor(X2)
names(groups) <- row.names(o.ifm)
scde.res2 <- scde.expression.difference(o.ifm, cd, o.prior, groups = groups, n.randomizations = 100, n.cores = 1, verbose = 1)

save(scde.res1, scde.res2, file = paste0("cell_type_copy", i, "_pairs",  j, "_scImpute_summary_scde"))




###################################################################
############# running DE analysis on MAGIC imputed data ###########
###################################################################



load("GSE75748_sc_cell_type.rda")
xx <- gene.expression
aa <- gsub("\\_[^<>]*", "", colnames(xx))
celltype <- unique(aa)


pairs <- matrix(c("EC", "DEC", "HFF", "DEC", "NPC", "DEC", "TB", "DEC",
"H9", "EC", "TB", "H1", "H9", "HFF", "H9", "NPC", "H9", "TB", "H9", "H1",  
"EC", "H1", "EC", "TB", "EC", "NPC",  "HFF", "TB", "HFF", "NPC",
"TB",  "NPC", "H1", "DEC",  "H9", "DEC",  "EC", "HFF",  "NPC", "TB",
"H1", "NPC"), byrow = T, ncol = 2)

load(paste0("Chu_cell_type_twosets_copy", i))
type1 <- aa[set1]
type2 <- aa[set2]


library(VIPER)

magic.res1 <- read.csv(paste0("Chu_cell_type_twosets_copy", i, "_t1_magic.csv"), header = T)
rownames(magic.res1) <- magic.res1[, 1]
gene.expression1 <- t(magic.res1[, -1])

magic.res2 <- read.csv(paste0("Chu_cell_type_twosets_copy", i, "_t2_magic.csv"), header = T)
rownames(magic.res2) <- magic.res2[, 1]
gene.expression2 <- t(magic.res2[, -1])


flag1.1 <- which(type1%in%pairs[j, 1])
flag1.2 <- which(type1%in%pairs[j, 2])
selected.x1 <- cbind(gene.expression1[, flag1.1], gene.expression1[, flag1.2])
X1 <- c(rep(0, length(flag1.1)), rep(1, length(flag1.2)))

flag2.1 <- which(type2%in%pairs[j, 1])
flag2.2 <- which(type2%in%pairs[j, 2])
selected.x2 <- cbind(gene.expression2[, flag2.1], gene.expression2[, flag2.2])
X2 <- c(rep(0, length(flag2.1)), rep(1, length(flag2.2)))


DESeq2.res1 <- DESeq2(selected.x1, matrix(X1, length(X1), 1))
DESeq2.res2 <- DESeq2(selected.x2, matrix(X2, length(X2), 1))
save(DESeq2.res1, DESeq2.res2, file = paste0("cell_type_copy", i, "_pairs",  j, "_magic_summary_DESeq2"))


edgeR.LRT.res1 <- edgeR_others(selected.x1, matrix(X1, length(X1), 1), method = "LRT")
edgeR.LRT.res2 <- edgeR_others(selected.x2, matrix(X2, length(X2), 1), method = "LRT")
save(edgeR.LRT.res1, edgeR.LRT.res2, file = paste0("cell_type_copy", i, "_pairs",  j, "_magic_summary_edgeR_LRT"))

edgeR.QLF.res1 <- edgeR_others(selected.x1, matrix(X1, length(X1), 1), method = "QLF")
edgeR.QLF.res2 <- edgeR_others(selected.x2, matrix(X2, length(X2), 1), method = "QLF")
save(edgeR.QLF.res1, edgeR.QLF.res2, file = paste0("cell_type_copy", i, "_pairs",  j, "_magic_summary_edgeR_QLF"))

library(scde)
sg <- factor(X1)
names(sg) <- colnames(selected.x1)
cd <- clean.counts(selected.x1, min.lib.size=100, min.reads = 1, min.detected = 1) # which is a gene by sample count matrix
counts <- apply(cd, 2, function(x) {storage.mode(x) <- 'integer'; x})
o.ifm <- scde.error.models(counts = counts, groups = sg, n.cores = 1, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1) 
valid.cells <- o.ifm$corr.a > 0
o.ifm <- o.ifm[valid.cells, ]
o.prior <- scde.expression.prior(models = o.ifm, counts = counts, length.out = 400, show.plot = FALSE)
groups <- factor(X1)
names(groups) <- row.names(o.ifm)
scde.res1 <- scde.expression.difference(o.ifm, cd, o.prior, groups = groups, n.randomizations = 100, n.cores = 1, verbose = 1)

sg <- factor(X2)
names(sg) <- colnames(selected.x2)
cd <- clean.counts(selected.x2, min.lib.size=100, min.reads = 1, min.detected = 1) # which is a gene by sample count matrix
counts <- apply(cd, 2, function(x) {storage.mode(x) <- 'integer'; x})
o.ifm <- scde.error.models(counts = counts, groups = sg, n.cores = 1, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1) 
valid.cells <- o.ifm$corr.a > 0
o.ifm <- o.ifm[valid.cells, ]
o.prior <- scde.expression.prior(models = o.ifm, counts = counts, length.out = 400, show.plot = FALSE)
groups <- factor(X2)
names(groups) <- row.names(o.ifm)
scde.res2 <- scde.expression.difference(o.ifm, cd, o.prior, groups = groups, n.randomizations = 100, n.cores = 1, verbose = 1)

save(scde.res1, scde.res2, file = paste0("cell_type_copy", i, "_pairs",  j, "_magic_summary_scde"))









###################################################################
############# running DE analysis on VIPER-lasso imputed data #####
###################################################################



load("GSE75748_sc_cell_type.rda")
xx <- gene.expression
aa <- gsub("\\_[^<>]*", "", colnames(xx))
celltype <- unique(aa)


pairs <- matrix(c("EC", "DEC", "HFF", "DEC", "NPC", "DEC", "TB", "DEC",
"H9", "EC", "TB", "H1", "H9", "HFF", "H9", "NPC", "H9", "TB", "H9", "H1",  
"EC", "H1", "EC", "TB", "EC", "NPC",  "HFF", "TB", "HFF", "NPC",
"TB",  "NPC", "H1", "DEC",  "H9", "DEC",  "EC", "HFF",  "NPC", "TB",
"H1", "NPC"), byrow = T, ncol = 2)


load(paste0("Chu_cell_type_twosets_copy", i))
type1 <- aa[set1]
type2 <- aa[set2]


library(VIPER)

load(paste0("Chu_cell_type_twosets_copy_lasso_copy", i))

gene.expression1 <- apply(imp_ee1$imputed, 2, function(x){
	exp(x) - 0.1
})

gene.expression2 <- apply(imp_ee2$imputed, 2, function(x){
	exp(x) - 0.1
})

flag1.1 <- which(type1%in%pairs[j, 1])
flag1.2 <- which(type1%in%pairs[j, 2])
selected.x1 <- cbind(gene.expression1[, flag1.1], gene.expression1[, flag1.2])
X1 <- c(rep(0, length(flag1.1)), rep(1, length(flag1.2)))

flag2.1 <- which(type2%in%pairs[j, 1])
flag2.2 <- which(type2%in%pairs[j, 2])
selected.x2 <- cbind(gene.expression2[, flag2.1], gene.expression2[, flag2.2])
X2 <- c(rep(0, length(flag2.1)), rep(1, length(flag2.2)))


DESeq2.res1 <- DESeq2(selected.x1, matrix(X1, length(X1), 1))
DESeq2.res2 <- DESeq2(selected.x2, matrix(X2, length(X2), 1))
save(DESeq2.res1, DESeq2.res2, file = paste0("cell_type_copy", i, "_pairs",  j, "_lasso_summary_DESeq2"))


edgeR.LRT.res1 <- edgeR_others(selected.x1, matrix(X1, length(X1), 1), method = "LRT")
edgeR.LRT.res2 <- edgeR_others(selected.x2, matrix(X2, length(X2), 1), method = "LRT")
save(edgeR.LRT.res1, edgeR.LRT.res2, file = paste0("cell_type_copy", i, "_pairs",  j, "_lasso_summary_edgeR_LRT"))

edgeR.QLF.res1 <- edgeR_others(selected.x1, matrix(X1, length(X1), 1), method = "QLF")
edgeR.QLF.res2 <- edgeR_others(selected.x2, matrix(X2, length(X2), 1), method = "QLF")
save(edgeR.QLF.res1, edgeR.QLF.res2, file = paste0("cell_type_copy", i, "_pairs",  j, "_lasso_summary_edgeR_QLF"))

library(scde)
sg <- factor(X1)
names(sg) <- colnames(selected.x1)
cd <- clean.counts(selected.x1, min.lib.size=100, min.reads = 1, min.detected = 1) # which is a gene by sample count matrix
counts <- apply(cd, 2, function(x) {storage.mode(x) <- 'integer'; x})
o.ifm <- scde.error.models(counts = counts, groups = sg, n.cores = 1, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1) 
valid.cells <- o.ifm$corr.a > 0
o.ifm <- o.ifm[valid.cells, ]
o.prior <- scde.expression.prior(models = o.ifm, counts = counts, length.out = 400, show.plot = FALSE)
groups <- factor(X1)
names(groups) <- row.names(o.ifm)
scde.res1 <- scde.expression.difference(o.ifm, cd, o.prior, groups = groups, n.randomizations = 100, n.cores = 1, verbose = 1)

sg <- factor(X2)
names(sg) <- colnames(selected.x2)
cd <- clean.counts(selected.x2, min.lib.size=100, min.reads = 1, min.detected = 1) # which is a gene by sample count matrix
counts <- apply(cd, 2, function(x) {storage.mode(x) <- 'integer'; x})
o.ifm <- scde.error.models(counts = counts, groups = sg, n.cores = 1, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1) 
valid.cells <- o.ifm$corr.a > 0
o.ifm <- o.ifm[valid.cells, ]
o.prior <- scde.expression.prior(models = o.ifm, counts = counts, length.out = 400, show.plot = FALSE)
groups <- factor(X2)
names(groups) <- row.names(o.ifm)
scde.res2 <- scde.expression.difference(o.ifm, cd, o.prior, groups = groups, n.randomizations = 100, n.cores = 1, verbose = 1)

save(scde.res1, scde.res2, file = paste0("cell_type_copy", i, "_pairs",  j, "_lasso_summary_scde"))




#########################################################################
############# running DE analysis on VIPER-Elastic net imputed data #####
#########################################################################


load("GSE75748_sc_cell_type.rda")
xx <- gene.expression
aa <- gsub("\\_[^<>]*", "", colnames(xx))
celltype <- unique(aa)


pairs <- matrix(c("EC", "DEC", "HFF", "DEC", "NPC", "DEC", "TB", "DEC",
"H9", "EC", "TB", "H1", "H9", "HFF", "H9", "NPC", "H9", "TB", "H9", "H1",  
"EC", "H1", "EC", "TB", "EC", "NPC",  "HFF", "TB", "HFF", "NPC",
"TB",  "NPC", "H1", "DEC",  "H9", "DEC",  "EC", "HFF",  "NPC", "TB",
"H1", "NPC"), byrow = T, ncol = 2)


load(paste0("Chu_cell_type_twosets_copy", i))
type1 <- aa[set1]
type2 <- aa[set2]


library(VIPER)

load(paste0("Chu_cell_type_twosets_copy_el_copy", i))

gene.expression1 <- apply(imp_el_ee1$imputed, 2, function(x){
	exp(x) - 0.1
})

gene.expression2 <- apply(imp_el_ee2$imputed, 2, function(x){
	exp(x) - 0.1
})

flag1.1 <- which(type1%in%pairs[j, 1])
flag1.2 <- which(type1%in%pairs[j, 2])
selected.x1 <- cbind(gene.expression1[, flag1.1], gene.expression1[, flag1.2])
X1 <- c(rep(0, length(flag1.1)), rep(1, length(flag1.2)))

flag2.1 <- which(type2%in%pairs[j, 1])
flag2.2 <- which(type2%in%pairs[j, 2])
selected.x2 <- cbind(gene.expression2[, flag2.1], gene.expression2[, flag2.2])
X2 <- c(rep(0, length(flag2.1)), rep(1, length(flag2.2)))


DESeq2.res1 <- DESeq2(selected.x1, matrix(X1, length(X1), 1))
DESeq2.res2 <- DESeq2(selected.x2, matrix(X2, length(X2), 1))
save(DESeq2.res1, DESeq2.res2, file = paste0("cell_type_copy", i, "_pairs",  j, "_el_summary_DESeq2"))


edgeR.LRT.res1 <- edgeR_others(selected.x1, matrix(X1, length(X1), 1), method = "LRT")
edgeR.LRT.res2 <- edgeR_others(selected.x2, matrix(X2, length(X2), 1), method = "LRT")
save(edgeR.LRT.res1, edgeR.LRT.res2, file = paste0("cell_type_copy", i, "_pairs",  j, "_el_summary_edgeR_LRT"))

edgeR.QLF.res1 <- edgeR_others(selected.x1, matrix(X1, length(X1), 1), method = "QLF")
edgeR.QLF.res2 <- edgeR_others(selected.x2, matrix(X2, length(X2), 1), method = "QLF")
save(edgeR.QLF.res1, edgeR.QLF.res2, file = paste0("cell_type_copy", i, "_pairs",  j, "_el_summary_edgeR_QLF"))

library(scde)
sg <- factor(X1)
names(sg) <- colnames(selected.x1)
cd <- clean.counts(selected.x1, min.lib.size=100, min.reads = 1, min.detected = 1) # which is a gene by sample count matrix
counts <- apply(cd, 2, function(x) {storage.mode(x) <- 'integer'; x})
o.ifm <- scde.error.models(counts = counts, groups = sg, n.cores = 1, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1) 
valid.cells <- o.ifm$corr.a > 0
o.ifm <- o.ifm[valid.cells, ]
o.prior <- scde.expression.prior(models = o.ifm, counts = counts, length.out = 400, show.plot = FALSE)
groups <- factor(X1)
names(groups) <- row.names(o.ifm)
scde.res1 <- scde.expression.difference(o.ifm, cd, o.prior, groups = groups, n.randomizations = 100, n.cores = 1, verbose = 1)

sg <- factor(X2)
names(sg) <- colnames(selected.x2)
cd <- clean.counts(selected.x2, min.lib.size=100, min.reads = 1, min.detected = 1) # which is a gene by sample count matrix
counts <- apply(cd, 2, function(x) {storage.mode(x) <- 'integer'; x})
o.ifm <- scde.error.models(counts = counts, groups = sg, n.cores = 1, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1) 
valid.cells <- o.ifm$corr.a > 0
o.ifm <- o.ifm[valid.cells, ]
o.prior <- scde.expression.prior(models = o.ifm, counts = counts, length.out = 400, show.plot = FALSE)
groups <- factor(X2)
names(groups) <- row.names(o.ifm)
scde.res2 <- scde.expression.difference(o.ifm, cd, o.prior, groups = groups, n.randomizations = 100, n.cores = 1, verbose = 1)

save(scde.res1, scde.res2, file = paste0("cell_type_copy", i, "_pairs",  j, "_el_summary_scde"))

