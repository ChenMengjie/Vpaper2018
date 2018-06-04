library(ggplot2)
library(easyGgplot2)


load("Shalek_imputation_with_missing2_summary")

df1 <- data.frame(case = rep(1:10, 6), method = c(rep("scImpute", 10), rep("DrImpute", 10), rep("SAVER", 10), rep("VIPER-Lasso", 10), rep("MAGIC", 10), rep("VIPER-Elastic Net", 10)), value = c(correlation.median.summary[, 3], correlation.median.summary[, 4], correlation.median.summary[, 5], correlation.median.summary[, 6], correlation.median.summary[, 7], correlation.median.summary[, 9]))

load("Shalek_imputation_with_missing5_summary")
df2 <- data.frame(case = rep(1:10, 6), method = c(rep("scImpute", 10), rep("DrImpute", 10), rep("SAVER", 10), rep("VIPER-Lasso", 10), rep("MAGIC", 10), rep("VIPER-Elastic Net", 10)), value = c(correlation.median.summary[, 3], correlation.median.summary[, 4], correlation.median.summary[, 5], correlation.median.summary[, 6], correlation.median.summary[, 7], correlation.median.summary[, 9]))

load("Shalek_imputation_with_missing10_summary")
df3 <- data.frame(case = rep(1:10, 6), method = c(rep("scImpute", 10), rep("DrImpute", 10), rep("SAVER", 10), rep("VIPER-Lasso", 10), rep("MAGIC", 10), rep("VIPER-Elastic Net", 10)), value = c(correlation.median.summary[, 3], correlation.median.summary[, 4], correlation.median.summary[, 5], correlation.median.summary[, 6], correlation.median.summary[, 7], correlation.median.summary[, 9]))

df <- cbind(missing = c(rep("Missing 2%", 60), rep("Missing 5%", 60), rep("Missing 10%", 60)), rbind(df1, df2, df3))
df$missing <- factor(df$missing, levels = c("Missing 2%", "Missing 5%", "Missing 10%"))
gg <- ggplot2.stripchart(data = df, xName = 'method', yName = 'value', groupName = 'method', backgroundColor="white",  groupColors = c('#66B2FF', '#00CC00', '#FFAAD4', '#B266FF', '#FFD4AA', '#6363EE'), xtitle="Imputation Methods", ytitle="Median Correlation", mainTitle = "Grun", removePanelGrid=TRUE, removePanelBorder=TRUE, boxplotFill="white", showLegend=TRUE, legendTitle = "Method", legendTitleFont = c(20, "bold", "black"), legendTextFont = c(20, "bold", "black"), setShapeByGroupName=TRUE, addBoxplot=TRUE, faceting=TRUE, facetingVarNames = c("missing"),  facetingDirection="horizontal", mainTitleFont = c(20, "bold", "black"), xtitleFont = c(20, "bold", "black"), size = 4, ytitleFont = c(20, "bold", "black"), xTickLabelFont = c(20, "bold", "white"), yTickLabelFont = c(20, "bold", "black"), facetingFont = c(20, "bold", "black"), facetingRect =  list(background = NULL, lineType = NULL, lineColor = NULL, lineSize = NULL)) 
		  	    
ggsave(gg, file = paste0("Shalek_missing_correlation_summary.pdf"), width = 12, height = 4)
	


#########
library(ggplot2)
library(easyGgplot2)

load("Shalek_imputation_with_missing2_summary")

df1 <- data.frame(case = rep(1:10, 6), method = c(rep("scImpute", 10), rep("DrImpute", 10), rep("SAVER", 10), rep("VIPER-Lasso", 10), rep("MAGIC", 10), rep("VIPER-Elastic Net", 10)), value = c(square.loss.median.summary[, 3], square.loss.median.summary[, 4], square.loss.median.summary[, 5], square.loss.median.summary[, 6], square.loss.median.summary[, 7], square.loss.median.summary[, 9]))

load("Shalek_imputation_with_missing5_summary")
df2 <- data.frame(case = rep(1:10, 6), method = c(rep("scImpute", 10), rep("DrImpute", 10), rep("SAVER", 10), rep("VIPER-Lasso", 10), rep("MAGIC", 10), rep("VIPER-Elastic Net", 10)), value = c(square.loss.median.summary[, 3], square.loss.median.summary[, 4], square.loss.median.summary[, 5], square.loss.median.summary[, 6], square.loss.median.summary[, 7], square.loss.median.summary[, 9]))

load("Shalek_imputation_with_missing10_summary")
df3 <- data.frame(case = rep(1:10, 6), method = c(rep("scImpute", 10), rep("DrImpute", 10), rep("SAVER", 10), rep("VIPER-Lasso", 10), rep("MAGIC", 10), rep("VIPER-Elastic Net", 10)), value = c(square.loss.median.summary[, 3], square.loss.median.summary[, 4], square.loss.median.summary[, 5], square.loss.median.summary[, 6], square.loss.median.summary[, 7], square.loss.median.summary[, 9]))

df <- cbind(missing = c(rep("Missing 2%", 60), rep("Missing 5%", 60), rep("Missing 10%", 60)), rbind(df1, df2, df3))
df$missing <- factor(df$missing, levels = c("Missing 2%", "Missing 5%", "Missing 10%"))
gg <- ggplot2.stripchart(data = df, xName = 'method', yName = 'value', groupName = 'method', backgroundColor="white",  groupColors = c('#66B2FF', '#00CC00', '#FFAAD4', '#B266FF', '#FFD4AA', '#6363EE'), xtitle="Imputation Methods", ytitle="Median Square Loss", mainTitle = "Grun", removePanelGrid=TRUE, removePanelBorder=TRUE, boxplotFill="white", showLegend=TRUE, legendTitle = "Method", legendTitleFont = c(20, "bold", "black"), legendTextFont = c(20, "bold", "black"), setShapeByGroupName=TRUE, addBoxplot=TRUE, faceting=TRUE, facetingVarNames = c("missing"),  facetingDirection="horizontal", mainTitleFont = c(20, "bold", "black"), xtitleFont = c(20, "bold", "black"), size = 4, ytitleFont = c(20, "bold", "black"), xTickLabelFont = c(20, "bold", "white"), yTickLabelFont = c(20, "bold", "black"), facetingFont = c(20, "bold", "black"), facetingRect =  list(background = NULL, lineType = NULL, lineColor = NULL, lineSize = NULL)) 
		  	    
ggsave(gg, file = paste0("Shalek_missing_square_loss_summary.pdf"), width = 12, height = 4)
	
	
	
#########
library(ggplot2)
library(easyGgplot2)

load("Shalek_imputation_with_missing2_summary")

df1 <- data.frame(case = rep(1:10, 6), method = c(rep("scImpute", 10), rep("DrImpute", 10), rep("SAVER", 10), rep("VIPER-Lasso", 10), rep("MAGIC", 10), rep("VIPER-Elastic Net", 10)), value = c(absolute.loss.median.summary[, 3], absolute.loss.median.summary[, 4], absolute.loss.median.summary[, 5], absolute.loss.median.summary[, 6], absolute.loss.median.summary[, 7], absolute.loss.median.summary[, 9]))

load("Shalek_imputation_with_missing5_summary")
df2 <- data.frame(case = rep(1:10, 6), method = c(rep("scImpute", 10), rep("DrImpute", 10), rep("SAVER", 10), rep("VIPER-Lasso", 10), rep("MAGIC", 10), rep("VIPER-Elastic Net", 10)), value = c(absolute.loss.median.summary[, 3], absolute.loss.median.summary[, 4], absolute.loss.median.summary[, 5], absolute.loss.median.summary[, 6], absolute.loss.median.summary[, 7], absolute.loss.median.summary[, 9]))

load("Shalek_imputation_with_missing10_summary")
df3 <- data.frame(case = rep(1:10, 6), method = c(rep("scImpute", 10), rep("DrImpute", 10), rep("SAVER", 10), rep("VIPER-Lasso", 10), rep("MAGIC", 10), rep("VIPER-Elastic Net", 10)), value = c(absolute.loss.median.summary[, 3], absolute.loss.median.summary[, 4], absolute.loss.median.summary[, 5], absolute.loss.median.summary[, 6], absolute.loss.median.summary[, 7], absolute.loss.median.summary[, 9]))

df <- cbind(missing = c(rep("Missing 2%", 60), rep("Missing 5%", 60), rep("Missing 10%", 60)), rbind(df1, df2, df3))
df$missing <- factor(df$missing, levels = c("Missing 2%", "Missing 5%", "Missing 10%"))
gg <- ggplot2.stripchart(data = df, xName = 'method', yName = 'value', groupName = 'method', backgroundColor="white",  groupColors = c('#66B2FF', '#00CC00', '#FFAAD4', '#B266FF', '#FFD4AA', '#6363EE'), xtitle="Imputation Methods", ytitle="Median L1 Loss", mainTitle = "Grun", removePanelGrid=TRUE, removePanelBorder=TRUE, boxplotFill="white", showLegend=TRUE, legendTitle = "Method", legendTitleFont = c(20, "bold", "black"), legendTextFont = c(20, "bold", "black"), setShapeByGroupName=TRUE, addBoxplot=TRUE, faceting=TRUE, facetingVarNames = c("missing"),  facetingDirection="horizontal", mainTitleFont = c(20, "bold", "black"), xtitleFont = c(20, "bold", "black"), size = 4, ytitleFont = c(20, "bold", "black"), xTickLabelFont = c(20, "bold", "white"), yTickLabelFont = c(20, "bold", "black"), facetingFont = c(20, "bold", "black"), facetingRect =  list(background = NULL, lineType = NULL, lineColor = NULL, lineSize = NULL)) 
			    
ggsave(gg, file = paste0("Shalek_missing_l1_loss_summary.pdf"), width = 12, height = 4)
		
	
