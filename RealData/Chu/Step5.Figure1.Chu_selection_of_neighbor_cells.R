library(gplots)

load("Chu_impute_tc_el")
res <- imp_ee$imputed
aa <- gsub("_[^<>]*", "", colnames(res))
	
weights.mat <- apply(imp_ee$sample_weights[!aa%in%"H9.00hb4s", !aa%in%"H9.00hb4s"], 2, function(x){
	 y <- abs(x)
	 y[y>0.01] <- 1
	 y
})

colorlist <- c("turquoise4", "cyan", "lavender",  "gold", "slateblue1", "violet", "skyblue1", "khaki", "pink", "salmon", "limegreen", "chocolate", "maroon", "purple", "blue", "yellow", "red",  "brown", '#FFAAD4', '#00CC00', '#66B2FF', '#B266FF', '#FFD4AA', '#AAFFD4', '#CC0000', "#B266CC")

col.color <- colorlist[as.numeric(as.factor(aa))][!aa%in%"H9.00hb4s"]

palette.gr.marray2 <- colorRampPalette(c("white", "red"))(4)
png("Chu_tc_el_selection.png", res=400, height = 8, width=8, unit = "in")
heatmap.2(as.matrix(weights.mat), trace = "none", col = palette.gr.marray2,
	          symbreaks = F, labCol = NA, dendrogram = "none", ColSideColors = col.color, labRow = NA, key = F, Colv = FALSE, Rowv = FALSE)
dev.off()

	
		
###########
load("Chu_impute_ct_el")
res <- imp_ee$imputed
aa <- gsub("_[^<>]*", "", colnames(res))

weights.mat <- apply(imp_ee$sample_weights, 2, function(x){
    y <- abs(x)
    y[y>0.1] <- 0.1
    y
})

colorlist <- c("turquoise4", "cyan", "lavender",  "gold", "slateblue1", "violet", "skyblue1", "khaki", "pink", "salmon", "limegreen", "chocolate", "maroon", "purple", "blue", "yellow", "red",  "brown", '#FFAAD4', '#00CC00', '#66B2FF', '#B266FF', '#FFD4AA', '#AAFFD4', '#CC0000', "#B266CC")

col.color <- colorlist[as.numeric(as.factor(aa))]

palette.gr.marray2 <- colorRampPalette(c("white", "red"))(4)
png("Chu_ct_el_selection.png",  res=400, height = 10, width=8, unit = "in")
heatmap.2(as.matrix(weights.mat), trace = "none", col = palette.gr.marray2,
            symbreaks = F, labCol = NA, dendrogram = "none", ColSideColors = col.color, 
            labRow = NA, key = T, Colv = FALSE, Rowv = FALSE)
dev.off()


