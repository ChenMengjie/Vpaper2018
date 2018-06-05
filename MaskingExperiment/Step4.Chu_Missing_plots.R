library(ggplot2)
library(easyGgplot2)


load("Chu_imputation_with_missing2_summary")

df1 <- data.frame(case = rep(1:10, 6), method = c(rep("scImpute", 10), rep("DrImpute", 10), rep("SAVER", 10), rep("VIPER-Lasso", 10), rep("MAGIC", 10), rep("VIPER-Elastic Net", 10)), value = c(correlation.median.summary[1:10, 3], correlation.median.summary[1:10, 4], correlation.median.summary[1:10, 5], correlation.median.summary[1:10, 6], correlation.median.summary[1:10, 7], correlation.median.summary[1:10, 8]))

load("Chu_imputation_with_missing5_summary")
df2 <- data.frame(case = rep(1:10, 6), method = c(rep("scImpute", 10), rep("DrImpute", 10), rep("SAVER", 10), rep("VIPER-Lasso", 10), rep("MAGIC", 10), rep("VIPER-Elastic Net", 10)), value = c(correlation.median.summary[1:10, 3], correlation.median.summary[1:10, 4], correlation.median.summary[1:10, 5], correlation.median.summary[1:10, 6], correlation.median.summary[1:10, 7], correlation.median.summary[1:10, 8]))

load("Chu_imputation_with_missing10_summary")
df3 <- data.frame(case = rep(1:10, 6), method = c(rep("scImpute", 10), rep("DrImpute", 10), rep("SAVER", 10), rep("VIPER-Lasso", 10), rep("MAGIC", 10), rep("VIPER-Elastic Net", 10)), value = c(correlation.median.summary[1:10, 3], correlation.median.summary[1:10, 4], correlation.median.summary[1:10, 5], correlation.median.summary[1:10, 6], correlation.median.summary[1:10, 7], correlation.median.summary[1:10, 8]))

df <- cbind(missing = c(rep("Missing 2%", 60), rep("Missing 5%", 60), rep("Missing 10%", 60)), rbind(df1, df2, df3))
df$missing <- factor(df$missing, levels = c("Missing 2%", "Missing 5%", "Missing 10%"))
gg <- ggplot2.stripchart(data = df, xName = 'method', yName = 'value', groupName = 'method', backgroundColor="white",  groupColors = c('#66B2FF', '#00CC00', '#FFAAD4', '#B266FF', '#FFD4AA', '#6363EE'), xtitle="Imputation Methods", ytitle="Median Correlation", mainTitle = "Grun", removePanelGrid=TRUE, removePanelBorder=TRUE, boxplotFill="white", showLegend=TRUE, legendTitle = "Method", legendTitleFont = c(20, "bold", "black"), legendTextFont = c(20, "bold", "black"), setShapeByGroupName=TRUE, addBoxplot=TRUE, faceting=TRUE, facetingVarNames = c("missing"),  facetingDirection="horizontal", mainTitleFont = c(20, "bold", "black"), xtitleFont = c(20, "bold", "black"), size = 4, ytitleFont = c(20, "bold", "black"), xTickLabelFont = c(20, "bold", "white"), yTickLabelFont = c(20, "bold", "black"), facetingFont = c(20, "bold", "black"), facetingRect =  list(background = NULL, lineType = NULL, lineColor = NULL, lineSize = NULL)) 
		  
		    
ggsave(gg, file = paste0("Chu_time_course_missing_correlation_summary.pdf"), width = 12, height = 4)
	


#########
library(ggplot2)
library(easyGgplot2)

load("Chu_imputation_with_missing2_summary")

df1 <- data.frame(case = rep(1:10, 6), method = c(rep("scImpute", 10), rep("DrImpute", 10), rep("SAVER", 10), rep("VIPER-Lasso", 10), rep("MAGIC", 10), rep("VIPER-Elastic Net", 10)), value = c(square.loss.median.summary[1:10, 3], square.loss.median.summary[1:10, 4], square.loss.median.summary[1:10, 5], square.loss.median.summary[1:10, 6], square.loss.median.summary[1:10, 7], square.loss.median.summary[1:10, 8]))

load("Chu_imputation_with_missing5_summary")
df2 <- data.frame(case = rep(1:10, 6), method = c(rep("scImpute", 10), rep("DrImpute", 10), rep("SAVER", 10), rep("VIPER-Lasso", 10), rep("MAGIC", 10), rep("VIPER-Elastic Net", 10)), value = c(square.loss.median.summary[1:10, 3], square.loss.median.summary[1:10, 4], square.loss.median.summary[1:10, 5], square.loss.median.summary[1:10, 6], square.loss.median.summary[1:10, 7], square.loss.median.summary[1:10, 8]))

load("Chu_imputation_with_missing10_summary")
df3 <- data.frame(case = rep(1:10, 6), method = c(rep("scImpute", 10), rep("DrImpute", 10), rep("SAVER", 10), rep("VIPER-Lasso", 10), rep("MAGIC", 10), rep("VIPER-Elastic Net", 10)), value = c(square.loss.median.summary[1:10, 3], square.loss.median.summary[1:10, 4], square.loss.median.summary[1:10, 5], square.loss.median.summary[1:10, 6], square.loss.median.summary[1:10, 7], square.loss.median.summary[1:10, 8]))

df <- cbind(missing = c(rep("Missing 2%", 60), rep("Missing 5%", 60), rep("Missing 10%", 60)), rbind(df1, df2, df3))
df$missing <- factor(df$missing, levels = c("Missing 2%", "Missing 5%", "Missing 10%"))
gg <- ggplot2.stripchart(data = df, xName = 'method', yName = 'value', groupName = 'method', backgroundColor="white",  groupColors = c('#66B2FF', '#00CC00', '#FFAAD4', '#B266FF', '#FFD4AA', '#6363EE'), xtitle="Imputation Methods", ytitle="Median Square Loss", mainTitle = "Grun", removePanelGrid=TRUE, removePanelBorder=TRUE, boxplotFill="white", showLegend=TRUE, legendTitle = "Method", legendTitleFont = c(20, "bold", "black"), legendTextFont = c(20, "bold", "black"), setShapeByGroupName=TRUE, addBoxplot=TRUE, faceting=TRUE, facetingVarNames = c("missing"),  facetingDirection="horizontal", mainTitleFont = c(20, "bold", "black"), xtitleFont = c(20, "bold", "black"), size = 4, ytitleFont = c(20, "bold", "black"), xTickLabelFont = c(20, "bold", "white"), yTickLabelFont = c(20, "bold", "black"), facetingFont = c(20, "bold", "black"), facetingRect =  list(background = NULL, lineType = NULL, lineColor = NULL, lineSize = NULL)) 
		  		    
ggsave(gg, file = paste0("Chu_time_course_missing_square_loss_summary.pdf"), width = 12, height = 4)
	
	
	
#########
library(ggplot2)
library(easyGgplot2)

load("Chu_imputation_with_missing2_summary")

df1 <- data.frame(case = rep(1:10, 6), method = c(rep("scImpute", 10), rep("DrImpute", 10), rep("SAVER", 10), rep("VIPER-Lasso", 10), rep("MAGIC", 10), rep("VIPER-Elastic Net", 10)), value = c(absolute.loss.median.summary[1:10, 3], absolute.loss.median.summary[1:10, 4], absolute.loss.median.summary[1:10, 5], absolute.loss.median.summary[1:10, 6], absolute.loss.median.summary[1:10, 7], absolute.loss.median.summary[1:10, 8]))

load("Chu_imputation_with_missing5_summary")
df2 <- data.frame(case = rep(1:10, 6), method = c(rep("scImpute", 10), rep("DrImpute", 10), rep("SAVER", 10), rep("VIPER-Lasso", 10), rep("MAGIC", 10), rep("VIPER-Elastic Net", 10)), value = c(absolute.loss.median.summary[1:10, 3], absolute.loss.median.summary[1:10, 4], absolute.loss.median.summary[1:10, 5], absolute.loss.median.summary[1:10, 6], absolute.loss.median.summary[1:10, 7], absolute.loss.median.summary[1:10, 8]))

load("Chu_imputation_with_missing10_summary")
df3 <- data.frame(case = rep(1:10, 6), method = c(rep("scImpute", 10), rep("DrImpute", 10), rep("SAVER", 10), rep("VIPER-Lasso", 10), rep("MAGIC", 10), rep("VIPER-Elastic Net", 10)), value = c(absolute.loss.median.summary[1:10, 3], absolute.loss.median.summary[1:10, 4], absolute.loss.median.summary[1:10, 5], absolute.loss.median.summary[1:10, 6], absolute.loss.median.summary[1:10, 7], absolute.loss.median.summary[1:10, 8]))

df <- cbind(missing = c(rep("Missing 2%", 60), rep("Missing 5%", 60), rep("Missing 10%", 60)), rbind(df1, df2, df3))
df$missing <- factor(df$missing, levels = c("Missing 2%", "Missing 5%", "Missing 10%"))
gg <- ggplot2.stripchart(data = df, xName = 'method', yName = 'value', groupName = 'method', backgroundColor="white",  groupColors = c('#66B2FF', '#00CC00', '#FFAAD4', '#B266FF', '#FFD4AA', '#6363EE'), xtitle="Imputation Methods", ytitle="Median L1 Loss", mainTitle = "Grun", removePanelGrid=TRUE, removePanelBorder=TRUE, boxplotFill="white", showLegend=TRUE, legendTitle = "Method", legendTitleFont = c(20, "bold", "black"), legendTextFont = c(20, "bold", "black"), setShapeByGroupName=TRUE, addBoxplot=TRUE, faceting=TRUE, facetingVarNames = c("missing"),  facetingDirection="horizontal", mainTitleFont = c(20, "bold", "black"), xtitleFont = c(20, "bold", "black"), size = 4, ytitleFont = c(20, "bold", "black"), xTickLabelFont = c(20, "bold", "white"), yTickLabelFont = c(20, "bold", "black"), facetingFont = c(20, "bold", "black"), facetingRect =  list(background = NULL, lineType = NULL, lineColor = NULL, lineSize = NULL)) 
		  	  
ggsave(gg, file = paste0("Chu_time_course_missing_l1_loss_summary.pdf"), width = 12, height = 4)
		
	
	
library(ggplot2)
library(easyGgplot2)


load("Chu_imputation_with_missing2_summary")

df1 <- data.frame(case = rep(11:20, 6), method = c(rep("scImpute", 10), rep("DrImpute", 10), rep("SAVER", 10), rep("VIPER-Lasso", 10), rep("MAGIC", 10), rep("VIPER-Elastic Net", 10)), value = c(correlation.median.summary[11:20, 3], correlation.median.summary[11:20, 4], correlation.median.summary[11:20, 5], correlation.median.summary[11:20, 6], correlation.median.summary[11:20, 7], correlation.median.summary[11:20, 8]))

load("Chu_imputation_with_missing5_summary")
df2 <- data.frame(case = rep(11:20, 6), method = c(rep("scImpute", 10), rep("DrImpute", 10), rep("SAVER", 10), rep("VIPER-Lasso", 10), rep("MAGIC", 10), rep("VIPER-Elastic Net", 10)), value = c(correlation.median.summary[11:20, 3], correlation.median.summary[11:20, 4], correlation.median.summary[11:20, 5], correlation.median.summary[11:20, 6], correlation.median.summary[11:20, 7], correlation.median.summary[11:20, 8]))

load("Chu_imputation_with_missing10_summary")
df3 <- data.frame(case = rep(11:20, 6), method = c(rep("scImpute", 10), rep("DrImpute", 10), rep("SAVER", 10), rep("VIPER-Lasso", 10), rep("MAGIC", 10), rep("VIPER-Elastic Net", 10)), value = c(correlation.median.summary[11:20, 3], correlation.median.summary[11:20, 4], correlation.median.summary[11:20, 5], correlation.median.summary[11:20, 6], correlation.median.summary[11:20, 7], correlation.median.summary[11:20, 8]))

df <- cbind(missing = c(rep("Missing 2%", 60), rep("Missing 5%", 60), rep("Missing 10%", 60)), rbind(df1, df2, df3))
df$missing <- factor(df$missing, levels = c("Missing 2%", "Missing 5%", "Missing 10%"))
gg <- ggplot2.stripchart(data = df, xName = 'method', yName = 'value', groupName = 'method', backgroundColor="white",  groupColors = c('#66B2FF', '#00CC00', '#FFAAD4', '#B266FF', '#FFD4AA', '#6363EE'), xtitle="Imputation Methods", ytitle="Median Correlation", mainTitle = "Grun", removePanelGrid=TRUE, removePanelBorder=TRUE, boxplotFill="white", showLegend=TRUE, legendTitle = "Method", legendTitleFont = c(20, "bold", "black"), legendTextFont = c(20, "bold", "black"), setShapeByGroupName=TRUE, addBoxplot=TRUE, faceting=TRUE, facetingVarNames = c("missing"),  facetingDirection="horizontal", mainTitleFont = c(20, "bold", "black"), xtitleFont = c(20, "bold", "black"), size = 4, ytitleFont = c(20, "bold", "black"), xTickLabelFont = c(20, "bold", "white"), yTickLabelFont = c(20, "bold", "black"), facetingFont = c(20, "bold", "black"), facetingRect =  list(background = NULL, lineType = NULL, lineColor = NULL, lineSize = NULL)) 
		  	  
ggsave(gg, file = paste0("Chu_cell_type_missing_correlation_summary.pdf"), width = 12, height = 4)
	


#########
library(ggplot2)
library(easyGgplot2)

load("Chu_imputation_with_missing2_summary")

df1 <- data.frame(case = rep(11:20, 6), method = c(rep("scImpute", 10), rep("DrImpute", 10), rep("SAVER", 10), rep("VIPER-Lasso", 10), rep("MAGIC", 10), rep("VIPER-Elastic Net", 10)), value = c(square.loss.median.summary[11:20, 3], square.loss.median.summary[11:20, 4], square.loss.median.summary[11:20, 5], square.loss.median.summary[11:20, 6], square.loss.median.summary[11:20, 7], square.loss.median.summary[11:20, 8]))

load("Chu_imputation_with_missing5_summary")
df2 <- data.frame(case = rep(11:20, 6), method = c(rep("scImpute", 10), rep("DrImpute", 10), rep("SAVER", 10), rep("VIPER-Lasso", 10), rep("MAGIC", 10), rep("VIPER-Elastic Net", 10)), value = c(square.loss.median.summary[11:20, 3], square.loss.median.summary[11:20, 4], square.loss.median.summary[11:20, 5], square.loss.median.summary[11:20, 6], square.loss.median.summary[11:20, 7], square.loss.median.summary[11:20, 8]))

load("Chu_imputation_with_missing10_summary")
df3 <- data.frame(case = rep(11:20, 6), method = c(rep("scImpute", 10), rep("DrImpute", 10), rep("SAVER", 10), rep("VIPER-Lasso", 10), rep("MAGIC", 10), rep("VIPER-Elastic Net", 10)), value = c(square.loss.median.summary[11:20, 3], square.loss.median.summary[11:20, 4], square.loss.median.summary[11:20, 5], square.loss.median.summary[11:20, 6], square.loss.median.summary[11:20, 7], square.loss.median.summary[11:20, 8]))

df <- cbind(missing = c(rep("Missing 2%", 60), rep("Missing 5%", 60), rep("Missing 10%", 60)), rbind(df1, df2, df3))
df$missing <- factor(df$missing, levels = c("Missing 2%", "Missing 5%", "Missing 10%"))
gg <- ggplot2.stripchart(data = df, xName = 'method', yName = 'value', groupName = 'method', backgroundColor="white",  groupColors = c('#66B2FF', '#00CC00', '#FFAAD4', '#B266FF', '#FFD4AA', '#6363EE'), xtitle="Imputation Methods", ytitle="Median Square Loss", mainTitle = "Grun", removePanelGrid=TRUE, removePanelBorder=TRUE, boxplotFill="white", showLegend=TRUE, legendTitle = "Method", legendTitleFont = c(20, "bold", "black"), legendTextFont = c(20, "bold", "black"), setShapeByGroupName=TRUE, addBoxplot=TRUE, faceting=TRUE, facetingVarNames = c("missing"),  facetingDirection="horizontal", mainTitleFont = c(20, "bold", "black"), xtitleFont = c(20, "bold", "black"), size = 4, ytitleFont = c(20, "bold", "black"), xTickLabelFont = c(20, "bold", "white"), yTickLabelFont = c(20, "bold", "black"), facetingFont = c(20, "bold", "black"), facetingRect =  list(background = NULL, lineType = NULL, lineColor = NULL, lineSize = NULL)) 
		  	
ggsave(gg, file = paste0("Chu_cell_type_missing_square_loss_summary.pdf"), width = 12, height = 4)
	
	
	
#########
library(ggplot2)
library(easyGgplot2)

load("Chu_imputation_with_missing2_summary")

df1 <- data.frame(case = rep(11:20, 6), method = c(rep("scImpute", 10), rep("DrImpute", 10), rep("SAVER", 10), rep("VIPER-Lasso", 10), rep("MAGIC", 10), rep("VIPER-Elastic Net", 10)), value = c(absolute.loss.median.summary[11:20, 3], absolute.loss.median.summary[11:20, 4], absolute.loss.median.summary[11:20, 5], absolute.loss.median.summary[11:20, 6], absolute.loss.median.summary[11:20, 7], absolute.loss.median.summary[11:20, 8]))

load("Chu_imputation_with_missing5_summary")
df2 <- data.frame(case = rep(11:20, 6), method = c(rep("scImpute", 10), rep("DrImpute", 10), rep("SAVER", 10), rep("VIPER-Lasso", 10), rep("MAGIC", 10), rep("VIPER-Elastic Net", 10)), value = c(absolute.loss.median.summary[11:20, 3], absolute.loss.median.summary[11:20, 4], absolute.loss.median.summary[11:20, 5], absolute.loss.median.summary[11:20, 6], absolute.loss.median.summary[11:20, 7], absolute.loss.median.summary[11:20, 8]))

load("Chu_imputation_with_missing10_summary")
df3 <- data.frame(case = rep(11:20, 6), method = c(rep("scImpute", 10), rep("DrImpute", 10), rep("SAVER", 10), rep("VIPER-Lasso", 10), rep("MAGIC", 10), rep("VIPER-Elastic Net", 10)), value = c(absolute.loss.median.summary[11:20, 3], absolute.loss.median.summary[11:20, 4], absolute.loss.median.summary[11:20, 5], absolute.loss.median.summary[11:20, 6], absolute.loss.median.summary[11:20, 7], absolute.loss.median.summary[11:20, 8]))

df <- cbind(missing = c(rep("Missing 2%", 60), rep("Missing 5%", 60), rep("Missing 10%", 60)), rbind(df1, df2, df3))
df$missing <- factor(df$missing, levels = c("Missing 2%", "Missing 5%", "Missing 10%"))
gg <- ggplot2.stripchart(data = df, xName = 'method', yName = 'value', groupName = 'method', backgroundColor="white",  groupColors = c('#66B2FF', '#00CC00', '#FFAAD4', '#B266FF', '#FFD4AA', '#6363EE'), xtitle="Imputation Methods", ytitle="Median L1 Loss", mainTitle = "Grun", removePanelGrid=TRUE, removePanelBorder=TRUE, boxplotFill="white", showLegend=TRUE, legendTitle = "Method", legendTitleFont = c(20, "bold", "black"), legendTextFont = c(20, "bold", "black"), setShapeByGroupName=TRUE, addBoxplot=TRUE, faceting=TRUE, facetingVarNames = c("missing"),  facetingDirection="horizontal", mainTitleFont = c(20, "bold", "black"), xtitleFont = c(20, "bold", "black"), size = 4, ytitleFont = c(20, "bold", "black"), xTickLabelFont = c(20, "bold", "white"), yTickLabelFont = c(20, "bold", "black"), facetingFont = c(20, "bold", "black"), facetingRect =  list(background = NULL, lineType = NULL, lineColor = NULL, lineSize = NULL)) 
   
ggsave(gg, file = paste0("Chu_cell_type_missing_l1_loss_summary.pdf"), width = 12, height = 4)
		
	

