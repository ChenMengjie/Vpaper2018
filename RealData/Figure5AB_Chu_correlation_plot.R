
library(ggplot2)
library(easyGgplot2)


load("impute_ct_heatmap_cor_summary")

df1 <- data.frame(case = rep(1:7, 7), method = c(rep("SAVER", 7), rep("VIPER-Lasso", 7), rep("VIPER-Elastic Net", 7), rep("MAGIC", 7), rep("scImpute", 7), rep("DrImpute", 7), rep("Before Imputation", 7)), value = c(cor.sum[, 1], cor.sum[, 2], cor.sum[, 3], cor.sum[, 4], cor.sum[, 5], cor.sum[, 6], cor.sum[, 7]))


gg <-  ggplot2.boxplot(data = df1, xName = 'method', yName = 'value', groupName = 'method', backgroundColor="white",  groupColors = c('#CC0000', '#66B2FF', '#00CC00', '#FFAAD4', '#B266FF', '#FFD4AA', '#6363EE'), xtitle="Imputation Methods", ytitle="Correlation", mainTitle = "Cell Type", removePanelGrid=TRUE, removePanelBorder=TRUE,  showLegend=TRUE, legendTitle = "", legendTitleFont = c(15, "bold", "black"), legendTextFont = c(15, "bold", "black"), mainTitleFont = c(15, "bold", "black"), xtitleFont = c(15, "bold", "black"), size = 4, ytitleFont = c(15, "bold", "black"), xTickLabelFont = c(15, "bold", "white"), yTickLabelFont = c(15, "bold", "black")) 
		    
ggsave(gg, file = paste0("impute_ct_heatmap_cor_summary.pdf"), width = 7, height = 4)
	



load("impute_tc_heatmap_cor_summary")

df2 <- data.frame(case = rep(1:7, 5), method = c(rep("SAVER", 5), rep("VIPER-Lasso", 5), rep("VIPER-Elastic Net", 5), rep("MAGIC", 5), rep("scImpute", 5), rep("DrImpute", 5), rep("Before Imputation", 5)), value = c(cor.sum[, 1], cor.sum[, 2], cor.sum[, 3], cor.sum[, 4], cor.sum[, 5], cor.sum[, 6], cor.sum[, 7]))


gg2 <-  ggplot2.boxplot(data = df2, xName = 'method', yName = 'value', groupName = 'method', backgroundColor="white",  groupColors = c('#CC0000', '#66B2FF', '#00CC00', '#FFAAD4', '#B266FF', '#FFD4AA', '#6363EE'), xtitle="Imputation Methods", ytitle="Correlation", mainTitle = "Time course", removePanelGrid=TRUE, removePanelBorder=TRUE,  showLegend=TRUE, legendTitle = "", legendTitleFont = c(15, "bold", "black"), legendTextFont = c(15, "bold", "black"),  mainTitleFont = c(15, "bold", "black"), xtitleFont = c(15, "bold", "black"), size = 4, ytitleFont = c(15, "bold", "black"), xTickLabelFont = c(15, "bold", "white"), yTickLabelFont = c(15, "bold", "black")) 
		    
ggsave(gg2, file = paste0("impute_tc_heatmap_cor_summary.pdf"), width = 7, height = 4)
	

