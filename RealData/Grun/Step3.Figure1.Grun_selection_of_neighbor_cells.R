library(gplots)
load("grun_impute_res")
res <- imp.res3$imputed

aa <- gsub("_[0-9]*", "", colnames(res))
	

	weights.mat <- apply(imp.res3$sample_weights, 2, function(x){
	  y <- x
	  y[y>0.01] <- 1
	  y[y<=0.01] <- 0
	  y
	})

	colorlist <- c("turquoise4", "cyan", "lavender",  "gold", "slateblue1", "violet", "skyblue1", "khaki", "pink", "salmon", "limegreen", "chocolate", "maroon", "purple", "blue", "yellow", "red",  "brown", '#FFAAD4', '#00CC00', '#66B2FF', '#B266FF', '#FFD4AA', '#AAFFD4', '#CC0000', "#B266CC")

	col.color <- colorlist[as.numeric(as.factor(aa))]

	palette.gr.marray2 <- colorRampPalette(c("white", "red"))(4)
	png("grun_impute_selection_.png", res=400, height = 8, width = 8, unit = "in")
	heatmap.2(as.matrix(weights.mat), trace = "none", col = palette.gr.marray2,
	          symbreaks = F, labCol = NA, dendrogram = "none", ColSideColors = col.color, labRow = NA, key = F, Colv = FALSE, Rowv = FALSE)
	dev.off()
	
	

colorlist <- c("turquoise4", "cyan", "lavender",  "gold")
plot(1:4, 1:4, col = colorlist)
legend("bottomright", c("SC 2i", "RNA 2i", "SC serum", "RNA serum"), col = colorlist, fill = colorlist, border = NA)



