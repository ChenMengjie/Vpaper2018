library(ggplot2)
library(easyGgplot2)

for(j in 1:5){
	
	load(paste0("cell_type_pairs",  j, "_allmethods_DEsummary"))
	df <- res
	aa <- as.character(df[, 1])
	df$Data <- aa
	df$Num <- apply(df, 1, function(x){
		(as.numeric(x[3])+as.numeric(x[4]))/2
	})
	
	gg <- ggplot2.stripchart(data = df, xName = 'Data', yName = 'Num', groupName = 'Data', backgroundColor="white",  groupColors = c('#CC0000', '#66B2FF', '#00CC00', '#FFAAD4', '#B266FF', '#FFD4AA', '#6363EE'), xtitle="Imputation Methods", ytitle="Number of DE genes", mainTitle = paste0(pairs[j, 1], " vs. ", pairs[j, 2]), removePanelGrid=TRUE, removePanelBorder=TRUE, setShapeByGroupName=TRUE, addBoxplot=TRUE, boxplotFill="white", showLegend=TRUE,  legendTitle = "Method", legendTitleFont = c(15, "bold", "black"), legendTextFont = c(15, "bold", "black"),  faceting=TRUE, facetingVarNames = c("Method"),  facetingDirection="horizontal", mainTitleFont = c(15, "bold", "black"), xtitleFont = c(15, "bold", "black"), size = 4, ytitleFont = c(15, "bold", "black"), xTickLabelFont = c(15, "bold", "white"), yTickLabelFont = c(15, "bold", "black"), facetingFont = c(15, "bold", "black"), facetingRect =  list(background = NULL, lineType = NULL, lineColor = NULL, lineSize = NULL)) 
		    
	ggsave(gg, file = paste0("cell_type_pairs",  j, "_comparison_DE_Number.pdf"), width = 12, height = 4)
	
}
