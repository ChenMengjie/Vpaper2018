load("Chu_impute_ct_res_summary")


load("GSE75748_sc_cell_type.rda")
celltype <- gsub("_[^<>]*", "", colnames(gene.expression))

cell.category <- unique(celltype)

logxx <- apply(gene.expression, 2, function(y){log(y + 1)})
logxx[gene.expression==0] <- 0

variance.nonzero <- function(data){

	V.per.gene <- apply(data, 1, function(x){
		var(x[x!=0])
	})
	V.per.gene
}	

variance.all <- function(data){
	apply(data, 1, var)
}	

mean.all <- function(data){
	apply(data, 1, mean)
}	

mean.nonzero <- function(data){
	apply(data, 1, function(x){
		mean(x[x!=0])
	})
}	

within.variance <- function(data, data.mean){
	
	y <- apply(data, 2, function(x){
		x - data.mean
	})
	apply(y, 1, function(z){
		sum(z^2)
	})
}


ff <- which(celltype%in%c("DEC", "H1"))

raw.v <- variance.all(logxx[, ff])
drimpute.v <- variance.all(log.drimpute.res[, ff])
scimpute.v <- variance.all(log.scimpute.res[, ff])
magic.v <- variance.all(log.magic.res[, ff])
saver.v <- variance.all(log.saver.res[, ff])
imputed.v <- variance.all(res[, ff])
imputed.el.v <- variance.all(log.res2[, ff])

set.seed(5)
show.3000 <- sample(1:length(drimpute.v), 3000)

library(ggplot2)
library(easyGgplot2)



	flag <- which(celltype%in%"H1")
	n1 <- length(flag)
	drimpute.v1 <- within.variance(log.drimpute.res[, flag], mean.all(log.drimpute.res[, flag]))
	scimpute.v1 <- within.variance(log.scimpute.res[, flag], mean.all(log.scimpute.res[, flag]))
	magic.v1 <- within.variance(log.magic.res[, flag], mean.all(log.magic.res[, flag]))
	saver.v1 <- within.variance(log.saver.res[, flag], mean.all(log.saver.res[, flag]))
	imputed.v1 <- within.variance(res[, flag], mean.all(res[, flag]))
	imputed.el.v1 <- within.variance(log.res2[, flag], mean.all(log.res2[, flag]))
	raw.v1 <- within.variance(logxx[, flag], mean.all(logxx[, flag]))
	
	

	flag <- which(celltype%in%"DEC")
	n2 <- length(flag)
	drimpute.v2 <- within.variance(log.drimpute.res[, flag], mean.all(log.drimpute.res[, flag]))
	scimpute.v2 <- within.variance(log.scimpute.res[, flag], mean.all(log.scimpute.res[, flag]))
	magic.v2 <- within.variance(log.magic.res[, flag], mean.all(log.magic.res[, flag]))
	saver.v2 <- within.variance(log.saver.res[, flag], mean.all(log.saver.res[, flag]))
	imputed.v2 <- within.variance(res[, flag], mean.all(res[, flag]))
	imputed.el.v2 <- within.variance(log.res2[, flag], mean.all(log.res2[, flag]))
	raw.v2 <- within.variance(logxx[, flag], mean.all(logxx[, flag]))
	
	
drimpute.between <- 1-(drimpute.v1+drimpute.v2)/(drimpute.v*(n1+n2))
scimpute.between <- 1-(scimpute.v1+scimpute.v2)/(scimpute.v*(n1+n2))
magic.between <- 1-(magic.v1+magic.v2)/(magic.v*(n1+n2))
saver.between <- 1-(saver.v1+saver.v2)/(saver.v*(n1+n2))
imputed.between <- 1-(imputed.v1+imputed.v2)/(imputed.v*(n1+n2))
imputed.el.between <- 1-(imputed.el.v1+imputed.el.v2)/(imputed.el.v*(n1+n2))
raw.between <- 1-(raw.v1+raw.v2)/(raw.v*(n1+n2))


	
	df <- data.frame(value = c(drimpute.between[show.3000], scimpute.between[show.3000], magic.between[show.3000], saver.between[show.3000], imputed.between[show.3000], imputed.el.between[show.3000]), method=c(rep("drImpute", 3000), rep("scImpute", 3000), rep("MAGIC", 3000), rep("SAVER", 3000), rep("VIPER-Lasso", 3000), rep("VIPER-Elastic Net", 3000)), compr = rep("Without-imputation", 3000*6), raw = rep(raw.between[show.3000], 6))
	
	ff <- ggplot2.scatterplot(data = df, xName = 'raw', yName = 'value', size = 0.3, backgroundColor = "white", xtitle="Between cluster sum of squares (Before Imputation)", ytitle="Between cluster sum of squares (after Imputation)", mainTitle = "H1 DEC", removePanelGrid=TRUE, color = "blue", removePanelBorder=FALSE, showLegend=FALSE, mainTitleFont = c(10, "bold", "black"), xtitleFont = c(15, "bold", "black"),  ytitleFont = c(15, "bold", "black"), xTickLabelFont = c(15, "bold", "white"), yTickLabelFont = c(15, "bold", "black"), facetingFont = c(15, "bold", "black"), facetingRect =  list(background = NULL, lineType = NULL, lineColor = NULL, lineSize = NULL)) + ylim(0, 1)  + xlim(0, 1) + facet_wrap(~method, ncol=3)  + theme(strip.text.x = element_text(size = 15, colour = "black", face = "bold")) + geom_abline(col = "brown", linetype = "dashed", size=1)

   	ggsave(ff, file = paste0("Chu_H1_DEC_sum_of_squares_comparison.pdf"), width = 12, height = 7)
	




ff <- which(celltype%in%c("EC", "HFF"))

raw.v <- variance.all(logxx[, ff])
drimpute.v <- variance.all(log.drimpute.res[, ff])
scimpute.v <- variance.all(log.scimpute.res[, ff])
magic.v <- variance.all(log.magic.res[, ff])
saver.v <- variance.all(log.saver.res[, ff])
imputed.v <- variance.all(res[, ff])
imputed.el.v <- variance.all(log.res2[, ff])

set.seed(5)
show.3000 <- sample(1:length(drimpute.v), 3000)

library(ggplot2)
library(easyGgplot2)



	flag <- which(celltype%in%"EC")
	n1 <- length(flag)
	drimpute.v1 <- within.variance(log.drimpute.res[, flag], mean.all(log.drimpute.res[, flag]))
	scimpute.v1 <- within.variance(log.scimpute.res[, flag], mean.all(log.scimpute.res[, flag]))
	magic.v1 <- within.variance(log.magic.res[, flag], mean.all(log.magic.res[, flag]))
	saver.v1 <- within.variance(log.saver.res[, flag], mean.all(log.saver.res[, flag]))
	imputed.v1 <- within.variance(res[, flag], mean.all(res[, flag]))
	imputed.el.v1 <- within.variance(log.res2[, flag], mean.all(log.res2[, flag]))
	raw.v1 <- within.variance(logxx[, flag], mean.all(logxx[, flag]))
	
	

	flag <- which(celltype%in%"HFF")
	n2 <- length(flag)
	drimpute.v2 <- within.variance(log.drimpute.res[, flag], mean.all(log.drimpute.res[, flag]))
	scimpute.v2 <- within.variance(log.scimpute.res[, flag], mean.all(log.scimpute.res[, flag]))
	magic.v2 <- within.variance(log.magic.res[, flag], mean.all(log.magic.res[, flag]))
	saver.v2 <- within.variance(log.saver.res[, flag], mean.all(log.saver.res[, flag]))
	imputed.v2 <- within.variance(res[, flag], mean.all(res[, flag]))
	imputed.el.v2 <- within.variance(log.res2[, flag], mean.all(log.res2[, flag]))
	raw.v2 <- within.variance(logxx[, flag], mean.all(logxx[, flag]))
	
	
drimpute.between <- 1-(drimpute.v1+drimpute.v2)/(drimpute.v*(n1+n2))
scimpute.between <- 1-(scimpute.v1+scimpute.v2)/(scimpute.v*(n1+n2))
magic.between <- 1-(magic.v1+magic.v2)/(magic.v*(n1+n2))
saver.between <- 1-(saver.v1+saver.v2)/(saver.v*(n1+n2))
imputed.between <- 1-(imputed.v1+imputed.v2)/(imputed.v*(n1+n2))
imputed.el.between <- 1-(imputed.el.v1+imputed.el.v2)/(imputed.el.v*(n1+n2))
raw.between <- 1-(raw.v1+raw.v2)/(raw.v*(n1+n2))


	
	df <- data.frame(value = c(drimpute.between[show.3000], scimpute.between[show.3000], magic.between[show.3000], saver.between[show.3000], imputed.between[show.3000], imputed.el.between[show.3000]), method=c(rep("drImpute", 3000), rep("scImpute", 3000), rep("MAGIC", 3000), rep("SAVER", 3000), rep("VIPER-Lasso", 3000), rep("VIPER-Elastic Net", 3000)), compr = rep("Without-imputation", 3000*6), raw = rep(raw.between[show.3000], 6))
	
	ff <- ggplot2.scatterplot(data = df, xName = 'raw', yName = 'value', size = 0.3, backgroundColor = "white", xtitle="Between cluster sum of squares (Before Imputation)", ytitle="Between cluster sum of squares (after Imputation)", mainTitle = "EC HFF", removePanelGrid=TRUE, color = "blue", removePanelBorder=FALSE, showLegend=FALSE, mainTitleFont = c(10, "bold", "black"), xtitleFont = c(15, "bold", "black"),  ytitleFont = c(15, "bold", "black"), xTickLabelFont = c(15, "bold", "white"), yTickLabelFont = c(15, "bold", "black"), facetingFont = c(15, "bold", "black"), facetingRect =  list(background = NULL, lineType = NULL, lineColor = NULL, lineSize = NULL)) + ylim(0, 1)  + xlim(0, 1) + facet_wrap(~method, ncol=3)  + theme(strip.text.x = element_text(size = 15, colour = "black", face = "bold")) + geom_abline(col = "brown", linetype = "dashed", size=1)

   	ggsave(ff, file = paste0("Chu_EC_HFF_sum_of_squares_comparison.pdf"), width = 12, height = 7)
	






ff <- which(celltype%in%c("NPC", "TB"))

raw.v <- variance.all(logxx[, ff])
drimpute.v <- variance.all(log.drimpute.res[, ff])
scimpute.v <- variance.all(log.scimpute.res[, ff])
magic.v <- variance.all(log.magic.res[, ff])
saver.v <- variance.all(log.saver.res[, ff])
imputed.v <- variance.all(res[, ff])
imputed.el.v <- variance.all(log.res2[, ff])

set.seed(5)
show.3000 <- sample(1:length(drimpute.v), 3000)

library(ggplot2)
library(easyGgplot2)



	flag <- which(celltype%in%"NPC")
	n1 <- length(flag)
	drimpute.v1 <- within.variance(log.drimpute.res[, flag], mean.all(log.drimpute.res[, flag]))
	scimpute.v1 <- within.variance(log.scimpute.res[, flag], mean.all(log.scimpute.res[, flag]))
	magic.v1 <- within.variance(log.magic.res[, flag], mean.all(log.magic.res[, flag]))
	saver.v1 <- within.variance(log.saver.res[, flag], mean.all(log.saver.res[, flag]))
	imputed.v1 <- within.variance(res[, flag], mean.all(res[, flag]))
	imputed.el.v1 <- within.variance(log.res2[, flag], mean.all(log.res2[, flag]))
	raw.v1 <- within.variance(logxx[, flag], mean.all(logxx[, flag]))
	
	

	flag <- which(celltype%in%"TB")
	n2 <- length(flag)
	drimpute.v2 <- within.variance(log.drimpute.res[, flag], mean.all(log.drimpute.res[, flag]))
	scimpute.v2 <- within.variance(log.scimpute.res[, flag], mean.all(log.scimpute.res[, flag]))
	magic.v2 <- within.variance(log.magic.res[, flag], mean.all(log.magic.res[, flag]))
	saver.v2 <- within.variance(log.saver.res[, flag], mean.all(log.saver.res[, flag]))
	imputed.v2 <- within.variance(res[, flag], mean.all(res[, flag]))
	imputed.el.v2 <- within.variance(log.res2[, flag], mean.all(log.res2[, flag]))
	raw.v2 <- within.variance(logxx[, flag], mean.all(logxx[, flag]))
	
	
drimpute.between <- 1-(drimpute.v1+drimpute.v2)/(drimpute.v*(n1+n2))
scimpute.between <- 1-(scimpute.v1+scimpute.v2)/(scimpute.v*(n1+n2))
magic.between <- 1-(magic.v1+magic.v2)/(magic.v*(n1+n2))
saver.between <- 1-(saver.v1+saver.v2)/(saver.v*(n1+n2))
imputed.between <- 1-(imputed.v1+imputed.v2)/(imputed.v*(n1+n2))
imputed.el.between <- 1-(imputed.el.v1+imputed.el.v2)/(imputed.el.v*(n1+n2))
raw.between <- 1-(raw.v1+raw.v2)/(raw.v*(n1+n2))


	
	df <- data.frame(value = c(drimpute.between[show.3000], scimpute.between[show.3000], magic.between[show.3000], saver.between[show.3000], imputed.between[show.3000], imputed.el.between[show.3000]), method=c(rep("drImpute", 3000), rep("scImpute", 3000), rep("MAGIC", 3000), rep("SAVER", 3000), rep("VIPER-Lasso", 3000), rep("VIPER-Elastic Net", 3000)), compr = rep("Without-imputation", 3000*6), raw = rep(raw.between[show.3000], 6))
	
	ff <- ggplot2.scatterplot(data = df, xName = 'raw', yName = 'value', size = 0.3, backgroundColor = "white", xtitle="Between cluster sum of squares (Before Imputation)", ytitle="Between cluster sum of squares (after Imputation)", mainTitle = "NPC TB", removePanelGrid=TRUE, color = "blue", removePanelBorder=FALSE, showLegend=FALSE, mainTitleFont = c(10, "bold", "black"), xtitleFont = c(15, "bold", "black"),  ytitleFont = c(15, "bold", "black"), xTickLabelFont = c(15, "bold", "white"), yTickLabelFont = c(15, "bold", "black"), facetingFont = c(15, "bold", "black"), facetingRect =  list(background = NULL, lineType = NULL, lineColor = NULL, lineSize = NULL)) + ylim(0, 1)  + xlim(0, 1) + facet_wrap(~method, ncol=3)  + theme(strip.text.x = element_text(size = 15, colour = "black", face = "bold")) + geom_abline(col = "brown", linetype = "dashed", size=1)

   	ggsave(ff, file = paste0("Chu_NPC_TB_sum_of_squares_comparison.pdf"), width = 12, height = 7)
	
