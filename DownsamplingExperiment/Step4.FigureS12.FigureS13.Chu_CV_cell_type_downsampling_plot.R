le <- 0.8

xx <- read.csv("GSE75748_sc_cell_type_ec.csv", as.is = T)
rownames(xx) <- xx[, 1]
xx <- xx[, -1]
aa <- gsub("\\_[^<>]*", "", colnames(xx))
types <- unique(aa)
gene.sum <- apply(xx, 1, sum)
xx <- xx[gene.sum > 0, ]

load(paste0("Chu_cell_type_filtered_uniform_rate", le))
sub.xx <- xx[rownames(gene.expression), colnames(gene.expression)]
	
load(paste0("Chu_cell_type_uniform_flag_rate", le))
	
sub.all.samples <- all.samples[rownames(gene.expression), colnames(gene.expression)]
all.samples.before.dr <- all.samples.before.dr[rownames(gene.expression), colnames(gene.expression)]
sub.all.ori <- all.ori[rownames(gene.expression), colnames(gene.expression)]
	
log.sub.xx <- apply(sub.xx, 2, function(x){log(x+0.1)})
	
load(paste0("Chu_cell_type_filtered_uniform_rate", le, "_comparison_imputation.summary"))

nonzero.flag <- sub.all.samples == 0 & sub.all.ori != 0  & all.samples.before.dr != 0
zero.flag <- sub.all.samples == 0 & sub.all.ori != 0  & all.samples.before.dr == 0

magic.res <- read.csv(paste0("Chu_cell_type_filtered_uniform_rate", le, "_magic.csv"), header = T)
rownames(magic.res) <- magic.res[, 1]
magic.res <- t(magic.res[, -1])
scimpute.res <- read.csv(paste0("Chu_cell_type_filtered_uniform_rate", le, ".csv"))
rownames(scimpute.res) <- scimpute.res[, 1]
scimpute.res <- scimpute.res[, -1]

log.scimpute.res <- apply(scimpute.res, 2, function(y){log(y + 0.1)})
log.magic.res <- apply(magic.res, 2, function(y){log(y + 0.1)})

load(paste0("Chu_drImpute_cell_type_filtered_uniform_rate", le))
load(paste0("Chu_SAVER_cell_type_filtered_uniform_rate", le))

log.drimpute.res <- lnpX_imp
log.saver.res <- apply(saver_imp$estimate, 2, function(y){log(y + 0.1)})

load(paste0("revision_Chu_cell_type_filtered_uniform_ee_rate", le))
res <- imp_ee$imputed

load(paste0("revision_Chu_cell_type_filtered_uniform_el_ee_rate", le))
log.res2 <- imp_ee$imputed

save(log.drimpute.res, log.scimpute.res, log.magic.res, log.saver.res, res, log.res2, file = paste0("revision_Chu_cell_type_uniform_summary", le))


load(paste0("revision_Chu_cell_type_uniform_summary", le))

CV.all <- function(data){
	CV <- function(mean, sd){
      abs(sd/mean)
	}
	CV.per.gene <- apply(data, 1, function(x){
		CV(mean(x), sd(x)) 
	})
	CV.per.gene
}	



log.sub.xx[sub.all.ori==0] <- 0
log.drimpute.res[sub.all.ori==0] <- 0
log.scimpute.res[sub.all.ori==0] <- 0
log.magic.res[sub.all.ori==0] <- 0
log.saver.res[sub.all.ori==0] <- 0
res[sub.all.ori==0] <- 0
log.res2[sub.all.ori==0] <- 0

CV.nonzero <- function(data){
	CV <- function(mean, sd){
      abs(sd/mean)
	}
	CV.per.gene <- apply(data, 1, function(x){
		CV(mean(x[x!=0]), sd(x[x!=0])) 
	})
	CV.per.gene
}	

celltype <- gsub("_[^<>]*", "", colnames(sub.all.samples))

cell.category <- unique(celltype)
set.seed(5)
show.3000 <- sample(1:nrow(log.sub.xx), 3000)

library(ggplot2)
library(easyGgplot2)


for(i in 1:7){
	
	flag <- which(celltype%in%cell.category[i])
	raw.data <- zero.flag[, flag]
	zero.num <- apply(raw.data, 1, function(x){
		length(x[x!=0])
	})
	zero.rate <- round(zero.num/ncol(raw.data), 2)*100
	raw.data2 <- nonzero.flag[, flag]
	dropout.num <- apply(raw.data2, 1, function(x){
		length(x[x!=0])
	})
	dropout.rate <- round(dropout.num/ncol(raw.data), 2)*100
	drimpute.CV <- CV.nonzero(log.drimpute.res[, flag])[show.3000]
	scimpute.CV <- CV.nonzero(log.scimpute.res[, flag])[show.3000]
	magic.CV <- CV.nonzero(log.magic.res[, flag])[show.3000]
	saver.CV <- CV.nonzero(log.saver.res[, flag])[show.3000]
	VIPER.CV <- CV.nonzero(res[, flag])[show.3000]
	VIPER.el.CV <- CV.nonzero(log.res2[, flag])[show.3000]
	raw.cv <- CV.nonzero(log.sub.xx[, flag])[show.3000]
	zero.rate.selected <- zero.rate[show.3000]
	dropout.rate.selected <- dropout.rate[show.3000]
	
	df <- data.frame(value = c(drimpute.CV, scimpute.CV, magic.CV, saver.CV, VIPER.CV, VIPER.el.CV), method=c(rep("drImpute", 3000), rep("scImpute", 3000), rep("MAGIC", 3000), rep("SAVER", 3000), rep("VIPER-Lasso", 3000), rep("VIPER-Elastic Net", 3000)), compr = rep("Without-imputation", 3000*6), raw = rep(raw.cv, 6), zero = rep(zero.rate.selected, 6), dropout = rep(dropout.rate.selected, 6))

	gg <- ggplot2.scatterplot(data = df, xName = 'raw', yName = 'value', size = 0.3, backgroundColor = "white", xtitle="CV (Original)", ytitle="CV (Imputed on Downsampled)", mainTitle = cell.category[i], removePanelGrid=TRUE,   removePanelBorder=FALSE, showLegend=TRUE, legendTitle = "Percentage \n of Down-sampling", legendTitleFont = c(15, "bold", "black"), legendTextFont = c(15, "bold", "black"), mainTitleFont = c(10, "bold", "black"), xtitleFont = c(15, "bold", "black"),  ytitleFont = c(15, "bold", "black"), xTickLabelFont = c(15, "bold", "white"), yTickLabelFont = c(15, "bold", "black"), facetingFont = c(15, "bold", "black"), facetingRect =  list(background = NULL, lineType = NULL, lineColor = NULL, lineSize = NULL)) + geom_point(aes(colour = zero), size = 0.3) + ylim(0, 5)  + xlim(0, 5) + facet_wrap(~method, ncol=3)  + theme(strip.text.x = element_text(size = 15, colour = "black", face = "bold")) + geom_abline(col = "brown", linetype = "dashed", size=1) 

   	ggsave(gg, file = paste0("Chu_cell_type_",  cell.category[i], "_CV_zero_downsampling_comparison.pdf"), width = 12, height = 7)
	

	gg <- ggplot2.scatterplot(data = df, xName = 'raw', yName = 'value', size = 0.3, backgroundColor = "white", xtitle="CV (Original)", ytitle="CV (Imputed on Downsampled)", mainTitle = cell.category[i], removePanelGrid=TRUE,   removePanelBorder=FALSE, showLegend=TRUE, legendTitle = "Percentage \n of Drop-out", legendTitleFont = c(15, "bold", "black"), legendTextFont = c(15, "bold", "black"), mainTitleFont = c(10, "bold", "black"), xtitleFont = c(15, "bold", "black"),  ytitleFont = c(15, "bold", "black"), xTickLabelFont = c(15, "bold", "white"), yTickLabelFont = c(15, "bold", "black"), facetingFont = c(15, "bold", "black"), facetingRect =  list(background = NULL, lineType = NULL, lineColor = NULL, lineSize = NULL)) + geom_point(aes(colour = dropout), size = 0.3) + ylim(0, 5)  + xlim(0, 5) + facet_wrap(~method, ncol=3)  + theme(strip.text.x = element_text(size = 15, colour = "black", face = "bold")) + geom_abline(col = "brown", linetype = "dashed", size=1) 

   	ggsave(gg, file = paste0("Chu_cell_type_",  cell.category[i], "_CV_zero_dropout_comparison.pdf"), width = 12, height = 7)
}



	  