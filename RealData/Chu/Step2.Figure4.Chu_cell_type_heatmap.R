magic.res <- read.csv("~/Downloads/Chu/GSE75748_sc_cell_type_magic.csv", header = T)
rownames(magic.res) <- magic.res[, 1]
magic.res <- t(magic.res[, -1])
scimpute.res <- read.csv("~/Downloads/Chu/scImpute2/scimpute_count.csv")
rownames(scimpute.res) <- scimpute.res[, 1]
scimpute.res <- scimpute.res[, -1]

log.scimpute.res <- apply(scimpute.res, 2, function(y){log(y + 0.1)})
log.magic.res <- apply(magic.res, 2, function(y){log(y + 0.1)})

load("~/Downloads/Chu/Chu_impute_ct_drImpute")
load("~/Downloads/Chu/Chu_impute_ct_SAVER")

log.drimpute.res <- lnpX_imp
log.saver.res <- apply(saver_imp$estimate, 2, function(y){log(y + 0.1)})

load("~/Downloads/Chu/Chu_impute_ct_ee_copy1")
res <- imp_ee$imputed

load("~/Downloads/Chu/Chu_impute_ct_el_ee_copy1")
log.res2 <- imp_ee$imputed

load("~/Downloads/Chu/GSE75748_sc_cell_type.rda")

logxx <- apply(gene.expression, 2, function(y){log(y + 0.1)})
nopredict <- logxx
nopredict[gene.expression==0] <- res[gene.expression==0]
nopredict2 <- logxx
nopredict2[gene.expression==0] <- log.res2[gene.expression==0]

aa <- gsub("_[^<>]*", "", colnames(res))


colorlist <- c("turquoise4", "cyan", "lavender",  "gold", "slateblue1", "violet", "skyblue1", "khaki", "pink", "salmon", "limegreen", "chocolate", "maroon", "purple", "blue", "yellow", "red",  "brown", '#FFAAD4', '#00CC00', '#66B2FF', '#B266FF', '#FFD4AA', '#AAFFD4', '#CC0000', "#B266CC")

col.color <- colorlist[as.numeric(as.factor(aa))]

library(gplots)

p <- nrow(gene.expression)
set.seed(937583495)
selected <- round(runif(2500)*p)
fll <- c(1:p)[selected]
palette.gr.marray <- colorRampPalette(c("blue", "white", "red"))(56)


png("/Users/mchen12/Downloads/Chu/Chu_ct_raw_heatmap2500.png", res=400, height = 10, width=8, unit = "in")
aa <- heatmap.2(as.matrix(logxx[fll, ]), trace = "none", col = palette.gr.marray,
                symbreaks = T, labCol = NA, dendrogram = "none", ColSideColors = col.color,
                labRow = NA, key = T, Colv = FALSE, Rowv = TRUE)
dev.off()

bb <- rev(aa$rowInd)
png("/Users/mchen12/Downloads/Chu/Chu_ct_nopredict_heatmap2500.png", res=400, height = 10, width=8, unit = "in")
heatmap.2(as.matrix(nopredict[fll, ][bb, ]), trace = "none", col = palette.gr.marray,
          symbreaks = T, labCol = NA, dendrogram = "none", ColSideColors = col.color,
          labRow = NA, key = T, Colv = FALSE, Rowv = FALSE)
dev.off()
png("/Users/mchen12/Downloads/Chu/Chu_ct_nopredict2_heatmap2500.png", res=400, height = 10, width=8, unit = "in")
heatmap.2(as.matrix(nopredict2[fll, ][bb, ]), trace = "none", col = palette.gr.marray,
          symbreaks = T, labCol = NA, dendrogram = "none", ColSideColors = col.color,
          labRow = NA, key = T, Colv = FALSE, Rowv = FALSE)
dev.off()
png("/Users/mchen12/Downloads/Chu/Chu_ct_scimpute_heatmap2500.png", res=400, height = 10, width=8, unit = "in")
heatmap.2(as.matrix(log.scimpute.res[fll, ][bb, ]), trace = "none", col = palette.gr.marray,
          symbreaks = T, labCol = NA, dendrogram = "none", ColSideColors = col.color,
          labRow = NA, key = T, Colv = FALSE, Rowv = FALSE)
dev.off()
png("/Users/mchen12/Downloads/Chu/Chu_ct_magic_heatmap2500.png", res=400, height = 10, width=8, unit = "in")
heatmap.2(as.matrix(log.magic.res[fll, ][bb, ]), trace = "none", col = palette.gr.marray,
          symbreaks = T, labCol = NA, dendrogram = "none", ColSideColors = col.color,
          labRow = NA, key = T, Colv = FALSE, Rowv = FALSE)
dev.off()
png("/Users/mchen12/Downloads/Chu/Chu_ct_saver_heatmap2500.png", res=400, height = 10, width=8, unit = "in")
heatmap.2(as.matrix(log.saver.res[fll, ][bb, ]), trace = "none", col = palette.gr.marray,
          symbreaks = T, labCol = NA, dendrogram = "none", ColSideColors = col.color,
          labRow = NA, key = T, Colv = FALSE, Rowv = FALSE)
dev.off()
png("/Users/mchen12/Downloads/Chu/Chu_ct_drimpute_heatmap2500.png", res=400, height = 10, width=8, unit = "in")
heatmap.2(as.matrix(log.drimpute.res[fll, ][bb, ]), trace = "none", col = palette.gr.marray,
          symbreaks = T, labCol = NA, dendrogram = "none", ColSideColors = col.color,
          labRow = NA, key = T, Colv = FALSE, Rowv = FALSE)
dev.off()



tt <- read.csv("~/Downloads/Chu/GSE75748_bulk_cell_type_ec.csv")
bulk <- tt[, -1]
rownames(bulk) <- tt[, 1]

selected.genes <- rownames(res)[fll]
sub.bulk <- bulk[selected.genes, ]

aa <- gsub("_re[^<>]*", "", colnames(bulk))
col.color <- colorlist[as.numeric(as.factor(aa))]

log.sub.bulk <- apply(sub.bulk, 2, function(y){log(y + 0.1)})


png(paste0("/Users/mchen12/Downloads/Chu/Chu_ct_bulk_heatmap.png"), res=400, height = 10, width=8, unit = "in")
heatmap.2(as.matrix(log.sub.bulk[bb, ]), trace = "none", col = palette.gr.marray,
          symbreaks = T, dendrogram = "none", ColSideColors = col.color,
          labRow = NA, key = T, Colv = FALSE, Rowv = FALSE)
dev.off()

log.bulk <- apply(bulk, 2, function(y){log(y + 0.1)})

save(log.drimpute.res, log.scimpute.res, log.magic.res, log.saver.res, log.bulk, res, log.res2, file = "/Users/mchen12/Downloads/Chu/impute_ct_res_summary")

