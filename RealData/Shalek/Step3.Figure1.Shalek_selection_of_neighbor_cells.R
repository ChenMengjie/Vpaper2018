
load("/project/mengjiechen/singlecell/summary/Shalek_impute_res")
res <- imp.res3$imputed

pheno <- read.delim("/project/mengjiechen/singlecell/data/Shalek/GSE48968-GPL13112-phenotype.txt", header = F)
rownames(pheno) <- pheno[, 1]	
aa <- pheno[colnames(res), 2]	

weights.mat <- apply(imp.res3$sample_weights, 2, function(x){
	  y <- x
	  y[y>0.01] <- 1
	  y[y<=0.01] <- 0
	  y
})

colorlist <- c("turquoise4", "cyan", "lavender",  "gold", "slateblue1", "violet", "skyblue1", "khaki", "pink", "salmon", "limegreen", "chocolate", "maroon", "purple", "blue", "yellow", "red",  "brown", '#FFAAD4', '#00CC00', '#66B2FF', '#B266FF', '#FFD4AA', '#AAFFD4', '#CC0000', "#B266CC")

col.color <- colorlist[as.numeric(as.factor(aa))]

library(gplots)
palette.gr.marray2 <- colorRampPalette(c("white", "red"))(4)
	png(paste0("/project/mengjiechen/singlecell/summary/Shalek_tc_selection_.png"), res=400, height = 8, width = 8, unit = "in")
	heatmap.2(as.matrix(weights.mat), trace = "none", col = palette.gr.marray2,
	          symbreaks = F, labCol = NA, dendrogram = "none", ColSideColors = col.color, labRow = NA, key = F, Colv = FALSE, Rowv = FALSE)
	dev.off()




colorlist <- c("blue","turquoise4", "cyan", "lavender",  "slateblue1", "purple", "gold", "violet", "skyblue1", "khaki", "pink", "salmon", "limegreen", "chocolate", "maroon")
pdf("testt.pdf")
plot(1:5, 1:5, col = colorlist[1:5], ylim=c(0,30), xlim = c(0,30))
legend("bottomright", names(table(aa)), col = colorlist[1:15], fill = colorlist, border = NA)
dev.off()




##########


library(gplots)
load("/project/mengjiechen/singlecell/summary/Shalek_impute_res")
res <- imp.res3$imputed

pheno <- read.delim("/project/mengjiechen/singlecell/data/Shalek/GSE48968-GPL13112-phenotype.txt", header = F)
rownames(pheno) <- pheno[, 1]	
aa <- pheno[colnames(res), 2]	

nn <- ncol(imp.res3$sample_weights)	
all.res <- imp.res3$sample_weights
all.res[all.res>-10] <- 0
colnames(all.res) <- rownames(all.res) <- colnames(imp.res3$imputed)
 
load("/project/mengjiechen/singlecell/data/Shalek/scImpute/Shalek_K10/Shalek_scImpute_2_weight_K10")

len <- length(scimpute_weight[[1]])

for(i in 1:len){
	
	mat <- scimpute_weight[[1]][i][[1]]
	sub.mat <- mat[-1, ]
	rownames(sub.mat) <- colnames(sub.mat)	
	all.res[rownames(sub.mat), rownames(sub.mat)] <- sub.mat

}


weights.mat <- apply(all.res, 2, function(x){
	  y <- x
	  y[y>0.01] <- 1
	  y[y<=0.01] <- 0
	  y
	})

	colorlist <- c("turquoise4", "cyan", "lavender",  "gold", "slateblue1", "violet", "skyblue1", "khaki", "pink", "salmon", "limegreen", "chocolate", "maroon", "purple", "blue", "yellow", "red",  "brown", '#FFAAD4', '#00CC00', '#66B2FF', '#B266FF', '#FFD4AA', '#AAFFD4', '#CC0000', "#B266CC")

	col.color <- colorlist[as.numeric(as.factor(aa))]

	palette.gr.marray2 <- colorRampPalette(c("white", "red"))(4)
	png("/project/mengjiechen/singlecell/summary/Shalek_scImpute_K10_selection.png", res=400, height = 8, width = 8, unit = "in")
	heatmap.2(as.matrix(weights.mat), trace = "none", col = palette.gr.marray2,
	          symbreaks = F, labCol = NA, dendrogram = "none", ColSideColors = col.color, labRow = NA, key = F, Colv = FALSE, Rowv = FALSE)
	dev.off()
	
	
	
	
	