###############################################################
############# DE analysis summary on raw count data ###########
###############################################################

jaccard <- function(x, y){
	length(intersect(x, y))/length(union(x, y))
}

load("GSE75748_sc_cell_type.rda")
xx <- gene.expression
aa <- gsub("\\_[^<>]*", "", colnames(xx))
celltype <- unique(aa)

pairs <- matrix(c("EC", "DEC", "HFF", "DEC", "NPC", "DEC", "TB", "DEC",
"H9", "EC", "TB", "H1", "H9", "HFF", "H9", "NPC", "H9", "TB", "H9", "H1",  
"EC", "H1", "EC", "TB", "EC", "NPC",  "HFF", "TB", "HFF", "NPC",
"TB",  "NPC", "H1", "DEC",  "H9", "DEC",  "EC", "HFF",  "NPC", "TB",
"H1", "NPC"), byrow = T, ncol = 2)

library(VIPER)

for(i in 1:10){
	for(j in 1:21){	
		

		load(paste0("cell_type_copy", i, "_pairs",  j, "_expression_summary_DESeq2"))
		load(paste0("cell_type_copy", i, "_pairs",  j, "_expression_summary_edgeR_LRT"))
		load(paste0("cell_type_copy", i, "_pairs",  j, "_expression_summary_edgeR_QLF"))
		load(paste0("cell_type_copy", i, "_pairs",  j, "_expression_summary_scde"))
				
		cutoff1.100 <- sort(DESeq2.res1$p.adjust)[100]
		DE.id1.100 <- which(DESeq2.res1$p.adjust <= cutoff1.100)		
		cutoff2.100 <- sort(DESeq2.res2$p.adjust)[100]
		DE.id2.100 <- which(DESeq2.res2$p.adjust <= cutoff2.100)
		
		cutoff1.200 <- sort(DESeq2.res1$p.adjust)[200]
		DE.id1.200 <- which(DESeq2.res1$p.adjust <= cutoff1.200)		
		cutoff2.200 <- sort(DESeq2.res2$p.adjust)[200]
		DE.id2.200 <- which(DESeq2.res2$p.adjust <= cutoff2.200)
		
		cutoff1.500 <- sort(DESeq2.res1$p.adjust)[500]
		DE.id1.500 <- which(DESeq2.res1$p.adjust <= cutoff1.500)		
		cutoff2.500 <- sort(DESeq2.res2$p.adjust)[500]
		DE.id2.500 <- which(DESeq2.res2$p.adjust <= cutoff2.500)
		
		cutoff1.1000 <- sort(DESeq2.res1$p.adjust)[1000]
		DE.id1.1000 <- which(DESeq2.res1$p.adjust <= cutoff1.1000)		
		cutoff2.1000 <- sort(DESeq2.res2$p.adjust)[1000]
		DE.id2.1000 <- which(DESeq2.res2$p.adjust <= cutoff2.1000)
		
		res1 <- c(jaccard(DE.id1.100, DE.id2.100), jaccard(DE.id1.200, DE.id2.200), jaccard(DE.id1.500, DE.id2.500), jaccard(DE.id1.1000, DE.id2.1000))
		
		
		cutoff1.100 <- sort(edgeR.LRT.res1$p.adjust)[100]
		edgeR.LRT.id1.100 <- which(edgeR.LRT.res1$p.adjust <= cutoff1.100)		
		cutoff2.100 <- sort(edgeR.LRT.res2$p.adjust)[100]
		edgeR.LRT.id2.100 <- which(edgeR.LRT.res2$p.adjust <= cutoff2.100)
		
		cutoff1.200 <- sort(edgeR.LRT.res1$p.adjust)[200]
		edgeR.LRT.id1.200 <- which(edgeR.LRT.res1$p.adjust <= cutoff1.200)		
		cutoff2.200 <- sort(edgeR.LRT.res2$p.adjust)[200]
		edgeR.LRT.id2.200 <- which(edgeR.LRT.res2$p.adjust <= cutoff2.200)
		
		cutoff1.500 <- sort(edgeR.LRT.res1$p.adjust)[500]
		edgeR.LRT.id1.500 <- which(edgeR.LRT.res1$p.adjust <= cutoff1.500)		
		cutoff2.500 <- sort(edgeR.LRT.res2$p.adjust)[500]
		edgeR.LRT.id2.500 <- which(edgeR.LRT.res2$p.adjust <= cutoff2.500)
		
		cutoff1.1000 <- sort(edgeR.LRT.res1$p.adjust)[1000]
		edgeR.LRT.id1.1000 <- which(edgeR.LRT.res1$p.adjust <= cutoff1.1000)		
		cutoff2.1000 <- sort(edgeR.LRT.res2$p.adjust)[1000]
		edgeR.LRT.id2.1000 <- which(edgeR.LRT.res2$p.adjust <= cutoff2.1000)
		
		res2 <- c(jaccard(edgeR.LRT.id1.100, edgeR.LRT.id2.100), jaccard(edgeR.LRT.id1.200, edgeR.LRT.id2.200), jaccard(edgeR.LRT.id1.500, edgeR.LRT.id2.500), jaccard(edgeR.LRT.id1.1000, edgeR.LRT.id2.1000))		

        cutoff1.100 <- sort(edgeR.QLF.res1$p.adjust)[100]
		edgeR.QLF.id1.100 <- which(edgeR.QLF.res1$p.adjust <= cutoff1.100)		
		cutoff2.100 <- sort(edgeR.QLF.res2$p.adjust)[100]
		edgeR.QLF.id2.100 <- which(edgeR.QLF.res2$p.adjust <= cutoff2.100)
		
		cutoff1.200 <- sort(edgeR.QLF.res1$p.adjust)[200]
		edgeR.QLF.id1.200 <- which(edgeR.QLF.res1$p.adjust <= cutoff1.200)		
		cutoff2.200 <- sort(edgeR.QLF.res2$p.adjust)[200]
		edgeR.QLF.id2.200 <- which(edgeR.QLF.res2$p.adjust <= cutoff2.200)
		
		cutoff1.500 <- sort(edgeR.QLF.res1$p.adjust)[500]
		edgeR.QLF.id1.500 <- which(edgeR.QLF.res1$p.adjust <= cutoff1.500)		
		cutoff2.500 <- sort(edgeR.QLF.res2$p.adjust)[500]
		edgeR.QLF.id2.500 <- which(edgeR.QLF.res2$p.adjust <= cutoff2.500)
		
		cutoff1.1000 <- sort(edgeR.QLF.res1$p.adjust)[1000]
		edgeR.QLF.id1.1000 <- which(edgeR.QLF.res1$p.adjust <= cutoff1.1000)		
		cutoff2.1000 <- sort(edgeR.QLF.res2$p.adjust)[1000]
		edgeR.QLF.id2.1000 <- which(edgeR.QLF.res2$p.adjust <= cutoff2.1000)
		
		res3 <- c(jaccard(edgeR.QLF.id1.100, edgeR.QLF.id2.100), jaccard(edgeR.QLF.id1.200, edgeR.QLF.id2.200), jaccard(edgeR.QLF.id1.500, edgeR.QLF.id2.500), jaccard(edgeR.QLF.id1.1000, edgeR.QLF.id2.1000))
		
		pvalue1 <- 2*(1-pnorm(abs(scde.res1$Z)))
		pvalue2 <- 2*(1-pnorm(abs(scde.res2$Z)))
		
		cutoff1.100 <- sort(pvalue1)[100]
		scde.id1.100 <- which(pvalue1 <= cutoff1.100)		
		cutoff2.100 <- sort(pvalue2)[100]
		scde.id2.100 <- which(pvalue2 <= cutoff2.100)
		
		cutoff1.200 <- sort(pvalue1)[200]
		scde.id1.200 <- which(pvalue1 <= cutoff1.200)		
		cutoff2.200 <- sort(pvalue2)[200]
		scde.id2.200 <- which(pvalue2 <= cutoff2.200)	
				
		cutoff1.500 <- sort(pvalue1)[500]
		scde.id1.500 <- which(pvalue1 <= cutoff1.500)		
		cutoff2.500 <- sort(pvalue2)[500]
		scde.id2.500 <- which(pvalue2 <= cutoff2.500)

       	cutoff1.1000 <- sort(pvalue1)[1000]
		scde.id1.1000 <- which(pvalue1 <= cutoff1.1000)		
		cutoff2.1000 <- sort(pvalue2)[1000]
		scde.id2.1000 <- which(pvalue2 <= cutoff2.1000)
		
		res4 <- c(jaccard(scde.id1.100, scde.id2.100), jaccard(scde.id1.200, scde.id2.200), jaccard(scde.id1.500, scde.id2.500), jaccard(scde.id1.1000, scde.id2.1000))
		
		DEsummary <- data.frame(Dataset = rep("Before Imputation", 4), Method = c("DESeq2", "edgeR_LRT", "edgeR_QLF", "SCDE"), rbind(res1, res2, res3, res4))
		colnames(DEsummary)[-c(1:2)] <- c("top 100", "top 200", "top 500", "top 1000")
		save(DEsummary, file = paste0("cell_type_copy", i, "_pairs",  j, "_expression_top_DEsummary"))
	
	}
}





#########################################################################
############# DE analysis summary on VIPER-lasso imputed data ###########
#########################################################################

load("GSE75748_sc_cell_type.rda")
xx <- gene.expression
aa <- gsub("\\_[^<>]*", "", colnames(xx))
celltype <- unique(aa)

pairs <- matrix(c("EC", "DEC", "HFF", "DEC", "NPC", "DEC", "TB", "DEC",
"H9", "EC", "TB", "H1", "H9", "HFF", "H9", "NPC", "H9", "TB", "H9", "H1",  
"EC", "H1", "EC", "TB", "EC", "NPC",  "HFF", "TB", "HFF", "NPC",
"TB",  "NPC", "H1", "DEC",  "H9", "DEC",  "EC", "HFF",  "NPC", "TB",
"H1", "NPC"), byrow = T, ncol = 2)


library(VIPER)

for(i in 1:10){
	for(j in 1:21){	
		

		load(paste0("cell_type_copy", i, "_pairs",  j, "_lasso_summary_DESeq2"))
		load(paste0("cell_type_copy", i, "_pairs",  j, "_lasso_summary_edgeR_LRT"))
		load(paste0("cell_type_copy", i, "_pairs",  j, "_lasso_summary_edgeR_QLF"))
		load(paste0("cell_type_copy", i, "_pairs",  j, "_lasso_summary_scde"))
				
		cutoff1.100 <- sort(DESeq2.res1$p.adjust)[100]
		DE.id1.100 <- which(DESeq2.res1$p.adjust <= cutoff1.100)		
		cutoff2.100 <- sort(DESeq2.res2$p.adjust)[100]
		DE.id2.100 <- which(DESeq2.res2$p.adjust <= cutoff2.100)
		
		cutoff1.200 <- sort(DESeq2.res1$p.adjust)[200]
		DE.id1.200 <- which(DESeq2.res1$p.adjust <= cutoff1.200)		
		cutoff2.200 <- sort(DESeq2.res2$p.adjust)[200]
		DE.id2.200 <- which(DESeq2.res2$p.adjust <= cutoff2.200)
		
		cutoff1.500 <- sort(DESeq2.res1$p.adjust)[500]
		DE.id1.500 <- which(DESeq2.res1$p.adjust <= cutoff1.500)		
		cutoff2.500 <- sort(DESeq2.res2$p.adjust)[500]
		DE.id2.500 <- which(DESeq2.res2$p.adjust <= cutoff2.500)
		
		cutoff1.1000 <- sort(DESeq2.res1$p.adjust)[1000]
		DE.id1.1000 <- which(DESeq2.res1$p.adjust <= cutoff1.1000)		
		cutoff2.1000 <- sort(DESeq2.res2$p.adjust)[1000]
		DE.id2.1000 <- which(DESeq2.res2$p.adjust <= cutoff2.1000)
		
		res1 <- c(jaccard(DE.id1.100, DE.id2.100), jaccard(DE.id1.200, DE.id2.200), jaccard(DE.id1.500, DE.id2.500), jaccard(DE.id1.1000, DE.id2.1000))
		
		
		cutoff1.100 <- sort(edgeR.LRT.res1$p.adjust)[100]
		edgeR.LRT.id1.100 <- which(edgeR.LRT.res1$p.adjust <= cutoff1.100)		
		cutoff2.100 <- sort(edgeR.LRT.res2$p.adjust)[100]
		edgeR.LRT.id2.100 <- which(edgeR.LRT.res2$p.adjust <= cutoff2.100)
		
		cutoff1.200 <- sort(edgeR.LRT.res1$p.adjust)[200]
		edgeR.LRT.id1.200 <- which(edgeR.LRT.res1$p.adjust <= cutoff1.200)		
		cutoff2.200 <- sort(edgeR.LRT.res2$p.adjust)[200]
		edgeR.LRT.id2.200 <- which(edgeR.LRT.res2$p.adjust <= cutoff2.200)
		
		cutoff1.500 <- sort(edgeR.LRT.res1$p.adjust)[500]
		edgeR.LRT.id1.500 <- which(edgeR.LRT.res1$p.adjust <= cutoff1.500)		
		cutoff2.500 <- sort(edgeR.LRT.res2$p.adjust)[500]
		edgeR.LRT.id2.500 <- which(edgeR.LRT.res2$p.adjust <= cutoff2.500)
		
		cutoff1.1000 <- sort(edgeR.LRT.res1$p.adjust)[1000]
		edgeR.LRT.id1.1000 <- which(edgeR.LRT.res1$p.adjust <= cutoff1.1000)		
		cutoff2.1000 <- sort(edgeR.LRT.res2$p.adjust)[1000]
		edgeR.LRT.id2.1000 <- which(edgeR.LRT.res2$p.adjust <= cutoff2.1000)
		
		res2 <- c(jaccard(edgeR.LRT.id1.100, edgeR.LRT.id2.100), jaccard(edgeR.LRT.id1.200, edgeR.LRT.id2.200), jaccard(edgeR.LRT.id1.500, edgeR.LRT.id2.500), jaccard(edgeR.LRT.id1.1000, edgeR.LRT.id2.1000))		

        cutoff1.100 <- sort(edgeR.QLF.res1$p.adjust)[100]
		edgeR.QLF.id1.100 <- which(edgeR.QLF.res1$p.adjust <= cutoff1.100)		
		cutoff2.100 <- sort(edgeR.QLF.res2$p.adjust)[100]
		edgeR.QLF.id2.100 <- which(edgeR.QLF.res2$p.adjust <= cutoff2.100)
		
		cutoff1.200 <- sort(edgeR.QLF.res1$p.adjust)[200]
		edgeR.QLF.id1.200 <- which(edgeR.QLF.res1$p.adjust <= cutoff1.200)		
		cutoff2.200 <- sort(edgeR.QLF.res2$p.adjust)[200]
		edgeR.QLF.id2.200 <- which(edgeR.QLF.res2$p.adjust <= cutoff2.200)
		
		cutoff1.500 <- sort(edgeR.QLF.res1$p.adjust)[500]
		edgeR.QLF.id1.500 <- which(edgeR.QLF.res1$p.adjust <= cutoff1.500)		
		cutoff2.500 <- sort(edgeR.QLF.res2$p.adjust)[500]
		edgeR.QLF.id2.500 <- which(edgeR.QLF.res2$p.adjust <= cutoff2.500)
		
		cutoff1.1000 <- sort(edgeR.QLF.res1$p.adjust)[1000]
		edgeR.QLF.id1.1000 <- which(edgeR.QLF.res1$p.adjust <= cutoff1.1000)		
		cutoff2.1000 <- sort(edgeR.QLF.res2$p.adjust)[1000]
		edgeR.QLF.id2.1000 <- which(edgeR.QLF.res2$p.adjust <= cutoff2.1000)
		
		res3 <- c(jaccard(edgeR.QLF.id1.100, edgeR.QLF.id2.100), jaccard(edgeR.QLF.id1.200, edgeR.QLF.id2.200), jaccard(edgeR.QLF.id1.500, edgeR.QLF.id2.500), jaccard(edgeR.QLF.id1.1000, edgeR.QLF.id2.1000))
		
		pvalue1 <- 2*(1-pnorm(abs(scde.res1$Z)))
		pvalue2 <- 2*(1-pnorm(abs(scde.res2$Z)))
		
		cutoff1.100 <- sort(pvalue1)[100]
		scde.id1.100 <- which(pvalue1 <= cutoff1.100)		
		cutoff2.100 <- sort(pvalue2)[100]
		scde.id2.100 <- which(pvalue2 <= cutoff2.100)
		
		cutoff1.200 <- sort(pvalue1)[200]
		scde.id1.200 <- which(pvalue1 <= cutoff1.200)		
		cutoff2.200 <- sort(pvalue2)[200]
		scde.id2.200 <- which(pvalue2 <= cutoff2.200)	
				
		cutoff1.500 <- sort(pvalue1)[500]
		scde.id1.500 <- which(pvalue1 <= cutoff1.500)		
		cutoff2.500 <- sort(pvalue2)[500]
		scde.id2.500 <- which(pvalue2 <= cutoff2.500)

       	cutoff1.1000 <- sort(pvalue1)[1000]
		scde.id1.1000 <- which(pvalue1 <= cutoff1.1000)		
		cutoff2.1000 <- sort(pvalue2)[1000]
		scde.id2.1000 <- which(pvalue2 <= cutoff2.1000)
		
		res4 <- c(jaccard(scde.id1.100, scde.id2.100), jaccard(scde.id1.200, scde.id2.200), jaccard(scde.id1.500, scde.id2.500), jaccard(scde.id1.1000, scde.id2.1000))
		
		DEsummary <- data.frame(Dataset = rep("VIPER-Lasso", 4), Method = c("DESeq2", "edgeR_LRT", "edgeR_QLF", "SCDE"), rbind(res1, res2, res3, res4))
		colnames(DEsummary)[-c(1:2)] <- c("top 100", "top 200", "top 500", "top 1000")
		save(DEsummary, file = paste0("cell_type_copy", i, "_pairs",  j, "_lasso_top_DEsummary"))
	
	}
}









###############################################################################
############# DE analysis summary on VIPER-elastic net imputed data ###########
###############################################################################

load("GSE75748_sc_cell_type.rda")
xx <- gene.expression
aa <- gsub("\\_[^<>]*", "", colnames(xx))
celltype <- unique(aa)

pairs <- matrix(c("EC", "DEC", "HFF", "DEC", "NPC", "DEC", "TB", "DEC",
"H9", "EC", "TB", "H1", "H9", "HFF", "H9", "NPC", "H9", "TB", "H9", "H1",  
"EC", "H1", "EC", "TB", "EC", "NPC",  "HFF", "TB", "HFF", "NPC",
"TB",  "NPC", "H1", "DEC",  "H9", "DEC",  "EC", "HFF",  "NPC", "TB",
"H1", "NPC"), byrow = T, ncol = 2)

library(VIPER)

for(i in 1:10){
	for(j in 1:21){	
		

		load(paste0("cell_type_copy", i, "_pairs",  j, "_el_summary_DESeq2"))
		load(paste0("cell_type_copy", i, "_pairs",  j, "_el_summary_edgeR_LRT"))
		load(paste0("cell_type_copy", i, "_pairs",  j, "_el_summary_edgeR_QLF"))
		load(paste0("cell_type_copy", i, "_pairs",  j, "_el_summary_scde"))
				
		cutoff1.100 <- sort(DESeq2.res1$p.adjust)[100]
		DE.id1.100 <- which(DESeq2.res1$p.adjust <= cutoff1.100)		
		cutoff2.100 <- sort(DESeq2.res2$p.adjust)[100]
		DE.id2.100 <- which(DESeq2.res2$p.adjust <= cutoff2.100)
		
		cutoff1.200 <- sort(DESeq2.res1$p.adjust)[200]
		DE.id1.200 <- which(DESeq2.res1$p.adjust <= cutoff1.200)		
		cutoff2.200 <- sort(DESeq2.res2$p.adjust)[200]
		DE.id2.200 <- which(DESeq2.res2$p.adjust <= cutoff2.200)
		
		cutoff1.500 <- sort(DESeq2.res1$p.adjust)[500]
		DE.id1.500 <- which(DESeq2.res1$p.adjust <= cutoff1.500)		
		cutoff2.500 <- sort(DESeq2.res2$p.adjust)[500]
		DE.id2.500 <- which(DESeq2.res2$p.adjust <= cutoff2.500)
		
		cutoff1.1000 <- sort(DESeq2.res1$p.adjust)[1000]
		DE.id1.1000 <- which(DESeq2.res1$p.adjust <= cutoff1.1000)		
		cutoff2.1000 <- sort(DESeq2.res2$p.adjust)[1000]
		DE.id2.1000 <- which(DESeq2.res2$p.adjust <= cutoff2.1000)
		
		res1 <- c(jaccard(DE.id1.100, DE.id2.100), jaccard(DE.id1.200, DE.id2.200), jaccard(DE.id1.500, DE.id2.500), jaccard(DE.id1.1000, DE.id2.1000))
		
		
		cutoff1.100 <- sort(edgeR.LRT.res1$p.adjust)[100]
		edgeR.LRT.id1.100 <- which(edgeR.LRT.res1$p.adjust <= cutoff1.100)		
		cutoff2.100 <- sort(edgeR.LRT.res2$p.adjust)[100]
		edgeR.LRT.id2.100 <- which(edgeR.LRT.res2$p.adjust <= cutoff2.100)
		
		cutoff1.200 <- sort(edgeR.LRT.res1$p.adjust)[200]
		edgeR.LRT.id1.200 <- which(edgeR.LRT.res1$p.adjust <= cutoff1.200)		
		cutoff2.200 <- sort(edgeR.LRT.res2$p.adjust)[200]
		edgeR.LRT.id2.200 <- which(edgeR.LRT.res2$p.adjust <= cutoff2.200)
		
		cutoff1.500 <- sort(edgeR.LRT.res1$p.adjust)[500]
		edgeR.LRT.id1.500 <- which(edgeR.LRT.res1$p.adjust <= cutoff1.500)		
		cutoff2.500 <- sort(edgeR.LRT.res2$p.adjust)[500]
		edgeR.LRT.id2.500 <- which(edgeR.LRT.res2$p.adjust <= cutoff2.500)
		
		cutoff1.1000 <- sort(edgeR.LRT.res1$p.adjust)[1000]
		edgeR.LRT.id1.1000 <- which(edgeR.LRT.res1$p.adjust <= cutoff1.1000)		
		cutoff2.1000 <- sort(edgeR.LRT.res2$p.adjust)[1000]
		edgeR.LRT.id2.1000 <- which(edgeR.LRT.res2$p.adjust <= cutoff2.1000)
		
		res2 <- c(jaccard(edgeR.LRT.id1.100, edgeR.LRT.id2.100), jaccard(edgeR.LRT.id1.200, edgeR.LRT.id2.200), jaccard(edgeR.LRT.id1.500, edgeR.LRT.id2.500), jaccard(edgeR.LRT.id1.1000, edgeR.LRT.id2.1000))		

        cutoff1.100 <- sort(edgeR.QLF.res1$p.adjust)[100]
		edgeR.QLF.id1.100 <- which(edgeR.QLF.res1$p.adjust <= cutoff1.100)		
		cutoff2.100 <- sort(edgeR.QLF.res2$p.adjust)[100]
		edgeR.QLF.id2.100 <- which(edgeR.QLF.res2$p.adjust <= cutoff2.100)
		
		cutoff1.200 <- sort(edgeR.QLF.res1$p.adjust)[200]
		edgeR.QLF.id1.200 <- which(edgeR.QLF.res1$p.adjust <= cutoff1.200)		
		cutoff2.200 <- sort(edgeR.QLF.res2$p.adjust)[200]
		edgeR.QLF.id2.200 <- which(edgeR.QLF.res2$p.adjust <= cutoff2.200)
		
		cutoff1.500 <- sort(edgeR.QLF.res1$p.adjust)[500]
		edgeR.QLF.id1.500 <- which(edgeR.QLF.res1$p.adjust <= cutoff1.500)		
		cutoff2.500 <- sort(edgeR.QLF.res2$p.adjust)[500]
		edgeR.QLF.id2.500 <- which(edgeR.QLF.res2$p.adjust <= cutoff2.500)
		
		cutoff1.1000 <- sort(edgeR.QLF.res1$p.adjust)[1000]
		edgeR.QLF.id1.1000 <- which(edgeR.QLF.res1$p.adjust <= cutoff1.1000)		
		cutoff2.1000 <- sort(edgeR.QLF.res2$p.adjust)[1000]
		edgeR.QLF.id2.1000 <- which(edgeR.QLF.res2$p.adjust <= cutoff2.1000)
		
		res3 <- c(jaccard(edgeR.QLF.id1.100, edgeR.QLF.id2.100), jaccard(edgeR.QLF.id1.200, edgeR.QLF.id2.200), jaccard(edgeR.QLF.id1.500, edgeR.QLF.id2.500), jaccard(edgeR.QLF.id1.1000, edgeR.QLF.id2.1000))
		
		pvalue1 <- 2*(1-pnorm(abs(scde.res1$Z)))
		pvalue2 <- 2*(1-pnorm(abs(scde.res2$Z)))
		
		cutoff1.100 <- sort(pvalue1)[100]
		scde.id1.100 <- which(pvalue1 <= cutoff1.100)		
		cutoff2.100 <- sort(pvalue2)[100]
		scde.id2.100 <- which(pvalue2 <= cutoff2.100)
		
		cutoff1.200 <- sort(pvalue1)[200]
		scde.id1.200 <- which(pvalue1 <= cutoff1.200)		
		cutoff2.200 <- sort(pvalue2)[200]
		scde.id2.200 <- which(pvalue2 <= cutoff2.200)	
				
		cutoff1.500 <- sort(pvalue1)[500]
		scde.id1.500 <- which(pvalue1 <= cutoff1.500)		
		cutoff2.500 <- sort(pvalue2)[500]
		scde.id2.500 <- which(pvalue2 <= cutoff2.500)

       	cutoff1.1000 <- sort(pvalue1)[1000]
		scde.id1.1000 <- which(pvalue1 <= cutoff1.1000)		
		cutoff2.1000 <- sort(pvalue2)[1000]
		scde.id2.1000 <- which(pvalue2 <= cutoff2.1000)
		
		res4 <- c(jaccard(scde.id1.100, scde.id2.100), jaccard(scde.id1.200, scde.id2.200), jaccard(scde.id1.500, scde.id2.500), jaccard(scde.id1.1000, scde.id2.1000))
		
		DEsummary <- data.frame(Dataset = rep("VIPER-Elastic Net", 4), Method = c("DESeq2", "edgeR_LRT", "edgeR_QLF", "SCDE"), rbind(res1, res2, res3, res4))
		colnames(DEsummary)[-c(1:2)] <- c("top 100", "top 200", "top 500", "top 1000")
		save(DEsummary, file = paste0("cell_type_copy", i, "_pairs",  j, "_el_top_DEsummary"))
	
	}
}




###############################################################################
############# DE analysis summary on drImpute imputed data ####################
###############################################################################




load("GSE75748_sc_cell_type.rda")
xx <- gene.expression
aa <- gsub("\\_[^<>]*", "", colnames(xx))
celltype <- unique(aa)



pairs <- matrix(c("EC", "DEC", "HFF", "DEC", "NPC", "DEC", "TB", "DEC",
"H9", "EC", "TB", "H1", "H9", "HFF", "H9", "NPC", "H9", "TB", "H9", "H1",  
"EC", "H1", "EC", "TB", "EC", "NPC",  "HFF", "TB", "HFF", "NPC",
"TB",  "NPC", "H1", "DEC",  "H9", "DEC",  "EC", "HFF",  "NPC", "TB",
"H1", "NPC"), byrow = T, ncol = 2)


library(VIPER)

for(i in 1:10){
	for(j in 1:21){	
		

		load(paste0("cell_type_copy", i, "_pairs",  j, "_drImpute_summary_DESeq2"))
		load(paste0("cell_type_copy", i, "_pairs",  j, "_drImpute_summary_edgeR_LRT"))
		load(paste0("cell_type_copy", i, "_pairs",  j, "_drImpute_summary_edgeR_QLF"))
		load(paste0("cell_type_copy", i, "_pairs",  j, "_drImpute_summary_scde"))
				
		cutoff1.100 <- sort(DESeq2.res1$p.adjust)[100]
		DE.id1.100 <- which(DESeq2.res1$p.adjust <= cutoff1.100)		
		cutoff2.100 <- sort(DESeq2.res2$p.adjust)[100]
		DE.id2.100 <- which(DESeq2.res2$p.adjust <= cutoff2.100)
		
		cutoff1.200 <- sort(DESeq2.res1$p.adjust)[200]
		DE.id1.200 <- which(DESeq2.res1$p.adjust <= cutoff1.200)		
		cutoff2.200 <- sort(DESeq2.res2$p.adjust)[200]
		DE.id2.200 <- which(DESeq2.res2$p.adjust <= cutoff2.200)
		
		cutoff1.500 <- sort(DESeq2.res1$p.adjust)[500]
		DE.id1.500 <- which(DESeq2.res1$p.adjust <= cutoff1.500)		
		cutoff2.500 <- sort(DESeq2.res2$p.adjust)[500]
		DE.id2.500 <- which(DESeq2.res2$p.adjust <= cutoff2.500)
		
		cutoff1.1000 <- sort(DESeq2.res1$p.adjust)[1000]
		DE.id1.1000 <- which(DESeq2.res1$p.adjust <= cutoff1.1000)		
		cutoff2.1000 <- sort(DESeq2.res2$p.adjust)[1000]
		DE.id2.1000 <- which(DESeq2.res2$p.adjust <= cutoff2.1000)
		
		res1 <- c(jaccard(DE.id1.100, DE.id2.100), jaccard(DE.id1.200, DE.id2.200), jaccard(DE.id1.500, DE.id2.500), jaccard(DE.id1.1000, DE.id2.1000))
		
		
		cutoff1.100 <- sort(edgeR.LRT.res1$p.adjust)[100]
		edgeR.LRT.id1.100 <- which(edgeR.LRT.res1$p.adjust <= cutoff1.100)		
		cutoff2.100 <- sort(edgeR.LRT.res2$p.adjust)[100]
		edgeR.LRT.id2.100 <- which(edgeR.LRT.res2$p.adjust <= cutoff2.100)
		
		cutoff1.200 <- sort(edgeR.LRT.res1$p.adjust)[200]
		edgeR.LRT.id1.200 <- which(edgeR.LRT.res1$p.adjust <= cutoff1.200)		
		cutoff2.200 <- sort(edgeR.LRT.res2$p.adjust)[200]
		edgeR.LRT.id2.200 <- which(edgeR.LRT.res2$p.adjust <= cutoff2.200)
		
		cutoff1.500 <- sort(edgeR.LRT.res1$p.adjust)[500]
		edgeR.LRT.id1.500 <- which(edgeR.LRT.res1$p.adjust <= cutoff1.500)		
		cutoff2.500 <- sort(edgeR.LRT.res2$p.adjust)[500]
		edgeR.LRT.id2.500 <- which(edgeR.LRT.res2$p.adjust <= cutoff2.500)
		
		cutoff1.1000 <- sort(edgeR.LRT.res1$p.adjust)[1000]
		edgeR.LRT.id1.1000 <- which(edgeR.LRT.res1$p.adjust <= cutoff1.1000)		
		cutoff2.1000 <- sort(edgeR.LRT.res2$p.adjust)[1000]
		edgeR.LRT.id2.1000 <- which(edgeR.LRT.res2$p.adjust <= cutoff2.1000)
		
		res2 <- c(jaccard(edgeR.LRT.id1.100, edgeR.LRT.id2.100), jaccard(edgeR.LRT.id1.200, edgeR.LRT.id2.200), jaccard(edgeR.LRT.id1.500, edgeR.LRT.id2.500), jaccard(edgeR.LRT.id1.1000, edgeR.LRT.id2.1000))		

        cutoff1.100 <- sort(edgeR.QLF.res1$p.adjust)[100]
		edgeR.QLF.id1.100 <- which(edgeR.QLF.res1$p.adjust <= cutoff1.100)		
		cutoff2.100 <- sort(edgeR.QLF.res2$p.adjust)[100]
		edgeR.QLF.id2.100 <- which(edgeR.QLF.res2$p.adjust <= cutoff2.100)
		
		cutoff1.200 <- sort(edgeR.QLF.res1$p.adjust)[200]
		edgeR.QLF.id1.200 <- which(edgeR.QLF.res1$p.adjust <= cutoff1.200)		
		cutoff2.200 <- sort(edgeR.QLF.res2$p.adjust)[200]
		edgeR.QLF.id2.200 <- which(edgeR.QLF.res2$p.adjust <= cutoff2.200)
		
		cutoff1.500 <- sort(edgeR.QLF.res1$p.adjust)[500]
		edgeR.QLF.id1.500 <- which(edgeR.QLF.res1$p.adjust <= cutoff1.500)		
		cutoff2.500 <- sort(edgeR.QLF.res2$p.adjust)[500]
		edgeR.QLF.id2.500 <- which(edgeR.QLF.res2$p.adjust <= cutoff2.500)
		
		cutoff1.1000 <- sort(edgeR.QLF.res1$p.adjust)[1000]
		edgeR.QLF.id1.1000 <- which(edgeR.QLF.res1$p.adjust <= cutoff1.1000)		
		cutoff2.1000 <- sort(edgeR.QLF.res2$p.adjust)[1000]
		edgeR.QLF.id2.1000 <- which(edgeR.QLF.res2$p.adjust <= cutoff2.1000)
		
		res3 <- c(jaccard(edgeR.QLF.id1.100, edgeR.QLF.id2.100), jaccard(edgeR.QLF.id1.200, edgeR.QLF.id2.200), jaccard(edgeR.QLF.id1.500, edgeR.QLF.id2.500), jaccard(edgeR.QLF.id1.1000, edgeR.QLF.id2.1000))
		
		pvalue1 <- 2*(1-pnorm(abs(scde.res1$Z)))
		pvalue2 <- 2*(1-pnorm(abs(scde.res2$Z)))
		
		cutoff1.100 <- sort(pvalue1)[100]
		scde.id1.100 <- which(pvalue1 <= cutoff1.100)		
		cutoff2.100 <- sort(pvalue2)[100]
		scde.id2.100 <- which(pvalue2 <= cutoff2.100)
		
		cutoff1.200 <- sort(pvalue1)[200]
		scde.id1.200 <- which(pvalue1 <= cutoff1.200)		
		cutoff2.200 <- sort(pvalue2)[200]
		scde.id2.200 <- which(pvalue2 <= cutoff2.200)	
				
		cutoff1.500 <- sort(pvalue1)[500]
		scde.id1.500 <- which(pvalue1 <= cutoff1.500)		
		cutoff2.500 <- sort(pvalue2)[500]
		scde.id2.500 <- which(pvalue2 <= cutoff2.500)

       	cutoff1.1000 <- sort(pvalue1)[1000]
		scde.id1.1000 <- which(pvalue1 <= cutoff1.1000)		
		cutoff2.1000 <- sort(pvalue2)[1000]
		scde.id2.1000 <- which(pvalue2 <= cutoff2.1000)
		
		res4 <- c(jaccard(scde.id1.100, scde.id2.100), jaccard(scde.id1.200, scde.id2.200), jaccard(scde.id1.500, scde.id2.500), jaccard(scde.id1.1000, scde.id2.1000))
		
		DEsummary <- data.frame(Dataset = rep("DrImpute", 4), Method = c("DESeq2", "edgeR_LRT", "edgeR_QLF", "SCDE"), rbind(res1, res2, res3, res4))
		colnames(DEsummary)[-c(1:2)] <- c("top 100", "top 200", "top 500", "top 1000")
		save(DEsummary, file = paste0("cell_type_copy", i, "_pairs",  j, "_DrImpute_top_DEsummary"))
	
	}
}







###############################################################################
############# DE analysis summary on SAVER imputed data #######################
###############################################################################


load("GSE75748_sc_cell_type.rda")
xx <- gene.expression
aa <- gsub("\\_[^<>]*", "", colnames(xx))
celltype <- unique(aa)



pairs <- matrix(c("EC", "DEC", "HFF", "DEC", "NPC", "DEC", "TB", "DEC",
"H9", "EC", "TB", "H1", "H9", "HFF", "H9", "NPC", "H9", "TB", "H9", "H1",  
"EC", "H1", "EC", "TB", "EC", "NPC",  "HFF", "TB", "HFF", "NPC",
"TB",  "NPC", "H1", "DEC",  "H9", "DEC",  "EC", "HFF",  "NPC", "TB",
"H1", "NPC"), byrow = T, ncol = 2)


library(VIPER)

for(i in 1:10){
	for(j in 1:21){	
		

		load(paste0("cell_type_copy", i, "_pairs",  j, "_SAVER_summary_DESeq2"))
		load(paste0("cell_type_copy", i, "_pairs",  j, "_SAVER_summary_edgeR_LRT"))
		load(paste0("cell_type_copy", i, "_pairs",  j, "_SAVER_summary_edgeR_QLF"))
		load(paste0("cell_type_copy", i, "_pairs",  j, "_SAVER_summary_scde"))
				
		cutoff1.100 <- sort(DESeq2.res1$p.adjust)[100]
		DE.id1.100 <- which(DESeq2.res1$p.adjust <= cutoff1.100)		
		cutoff2.100 <- sort(DESeq2.res2$p.adjust)[100]
		DE.id2.100 <- which(DESeq2.res2$p.adjust <= cutoff2.100)
		
		cutoff1.200 <- sort(DESeq2.res1$p.adjust)[200]
		DE.id1.200 <- which(DESeq2.res1$p.adjust <= cutoff1.200)		
		cutoff2.200 <- sort(DESeq2.res2$p.adjust)[200]
		DE.id2.200 <- which(DESeq2.res2$p.adjust <= cutoff2.200)
		
		cutoff1.500 <- sort(DESeq2.res1$p.adjust)[500]
		DE.id1.500 <- which(DESeq2.res1$p.adjust <= cutoff1.500)		
		cutoff2.500 <- sort(DESeq2.res2$p.adjust)[500]
		DE.id2.500 <- which(DESeq2.res2$p.adjust <= cutoff2.500)
		
		cutoff1.1000 <- sort(DESeq2.res1$p.adjust)[1000]
		DE.id1.1000 <- which(DESeq2.res1$p.adjust <= cutoff1.1000)		
		cutoff2.1000 <- sort(DESeq2.res2$p.adjust)[1000]
		DE.id2.1000 <- which(DESeq2.res2$p.adjust <= cutoff2.1000)
		
		res1 <- c(jaccard(DE.id1.100, DE.id2.100), jaccard(DE.id1.200, DE.id2.200), jaccard(DE.id1.500, DE.id2.500), jaccard(DE.id1.1000, DE.id2.1000))
		
		
		cutoff1.100 <- sort(edgeR.LRT.res1$p.adjust)[100]
		edgeR.LRT.id1.100 <- which(edgeR.LRT.res1$p.adjust <= cutoff1.100)		
		cutoff2.100 <- sort(edgeR.LRT.res2$p.adjust)[100]
		edgeR.LRT.id2.100 <- which(edgeR.LRT.res2$p.adjust <= cutoff2.100)
		
		cutoff1.200 <- sort(edgeR.LRT.res1$p.adjust)[200]
		edgeR.LRT.id1.200 <- which(edgeR.LRT.res1$p.adjust <= cutoff1.200)		
		cutoff2.200 <- sort(edgeR.LRT.res2$p.adjust)[200]
		edgeR.LRT.id2.200 <- which(edgeR.LRT.res2$p.adjust <= cutoff2.200)
		
		cutoff1.500 <- sort(edgeR.LRT.res1$p.adjust)[500]
		edgeR.LRT.id1.500 <- which(edgeR.LRT.res1$p.adjust <= cutoff1.500)		
		cutoff2.500 <- sort(edgeR.LRT.res2$p.adjust)[500]
		edgeR.LRT.id2.500 <- which(edgeR.LRT.res2$p.adjust <= cutoff2.500)
		
		cutoff1.1000 <- sort(edgeR.LRT.res1$p.adjust)[1000]
		edgeR.LRT.id1.1000 <- which(edgeR.LRT.res1$p.adjust <= cutoff1.1000)		
		cutoff2.1000 <- sort(edgeR.LRT.res2$p.adjust)[1000]
		edgeR.LRT.id2.1000 <- which(edgeR.LRT.res2$p.adjust <= cutoff2.1000)
		
		res2 <- c(jaccard(edgeR.LRT.id1.100, edgeR.LRT.id2.100), jaccard(edgeR.LRT.id1.200, edgeR.LRT.id2.200), jaccard(edgeR.LRT.id1.500, edgeR.LRT.id2.500), jaccard(edgeR.LRT.id1.1000, edgeR.LRT.id2.1000))		

        cutoff1.100 <- sort(edgeR.QLF.res1$p.adjust)[100]
		edgeR.QLF.id1.100 <- which(edgeR.QLF.res1$p.adjust <= cutoff1.100)		
		cutoff2.100 <- sort(edgeR.QLF.res2$p.adjust)[100]
		edgeR.QLF.id2.100 <- which(edgeR.QLF.res2$p.adjust <= cutoff2.100)
		
		cutoff1.200 <- sort(edgeR.QLF.res1$p.adjust)[200]
		edgeR.QLF.id1.200 <- which(edgeR.QLF.res1$p.adjust <= cutoff1.200)		
		cutoff2.200 <- sort(edgeR.QLF.res2$p.adjust)[200]
		edgeR.QLF.id2.200 <- which(edgeR.QLF.res2$p.adjust <= cutoff2.200)
		
		cutoff1.500 <- sort(edgeR.QLF.res1$p.adjust)[500]
		edgeR.QLF.id1.500 <- which(edgeR.QLF.res1$p.adjust <= cutoff1.500)		
		cutoff2.500 <- sort(edgeR.QLF.res2$p.adjust)[500]
		edgeR.QLF.id2.500 <- which(edgeR.QLF.res2$p.adjust <= cutoff2.500)
		
		cutoff1.1000 <- sort(edgeR.QLF.res1$p.adjust)[1000]
		edgeR.QLF.id1.1000 <- which(edgeR.QLF.res1$p.adjust <= cutoff1.1000)		
		cutoff2.1000 <- sort(edgeR.QLF.res2$p.adjust)[1000]
		edgeR.QLF.id2.1000 <- which(edgeR.QLF.res2$p.adjust <= cutoff2.1000)
		
		res3 <- c(jaccard(edgeR.QLF.id1.100, edgeR.QLF.id2.100), jaccard(edgeR.QLF.id1.200, edgeR.QLF.id2.200), jaccard(edgeR.QLF.id1.500, edgeR.QLF.id2.500), jaccard(edgeR.QLF.id1.1000, edgeR.QLF.id2.1000))
		
		pvalue1 <- 2*(1-pnorm(abs(scde.res1$Z)))
		pvalue2 <- 2*(1-pnorm(abs(scde.res2$Z)))
		
		cutoff1.100 <- sort(pvalue1)[100]
		scde.id1.100 <- which(pvalue1 <= cutoff1.100)		
		cutoff2.100 <- sort(pvalue2)[100]
		scde.id2.100 <- which(pvalue2 <= cutoff2.100)
		
		cutoff1.200 <- sort(pvalue1)[200]
		scde.id1.200 <- which(pvalue1 <= cutoff1.200)		
		cutoff2.200 <- sort(pvalue2)[200]
		scde.id2.200 <- which(pvalue2 <= cutoff2.200)	
				
		cutoff1.500 <- sort(pvalue1)[500]
		scde.id1.500 <- which(pvalue1 <= cutoff1.500)		
		cutoff2.500 <- sort(pvalue2)[500]
		scde.id2.500 <- which(pvalue2 <= cutoff2.500)

       	cutoff1.1000 <- sort(pvalue1)[1000]
		scde.id1.1000 <- which(pvalue1 <= cutoff1.1000)		
		cutoff2.1000 <- sort(pvalue2)[1000]
		scde.id2.1000 <- which(pvalue2 <= cutoff2.1000)
		
		res4 <- c(jaccard(scde.id1.100, scde.id2.100), jaccard(scde.id1.200, scde.id2.200), jaccard(scde.id1.500, scde.id2.500), jaccard(scde.id1.1000, scde.id2.1000))
		
		DEsummary <- data.frame(Dataset = rep("SAVER", 4), Method = c("DESeq2", "edgeR_LRT", "edgeR_QLF", "SCDE"), rbind(res1, res2, res3, res4))
		colnames(DEsummary)[-c(1:2)] <- c("top 100", "top 200", "top 500", "top 1000")
		save(DEsummary, file = paste0("cell_type_copy", i, "_pairs",  j, "_SAVER_top_DEsummary"))
	
	}
}





###############################################################################
############# DE analysis summary on scImpute imputed data ####################
###############################################################################


load("GSE75748_sc_cell_type.rda")
xx <- gene.expression
aa <- gsub("\\_[^<>]*", "", colnames(xx))
celltype <- unique(aa)



pairs <- matrix(c("EC", "DEC", "HFF", "DEC", "NPC", "DEC", "TB", "DEC",
"H9", "EC", "TB", "H1", "H9", "HFF", "H9", "NPC", "H9", "TB", "H9", "H1",  
"EC", "H1", "EC", "TB", "EC", "NPC",  "HFF", "TB", "HFF", "NPC",
"TB",  "NPC", "H1", "DEC",  "H9", "DEC",  "EC", "HFF",  "NPC", "TB",
"H1", "NPC"), byrow = T, ncol = 2)



library(VIPER)

for(i in 1:10){
	for(j in 1:21){	
		

		load(paste0("cell_type_copy", i, "_pairs",  j, "_scImpute_summary_DESeq2"))
		load(paste0("cell_type_copy", i, "_pairs",  j, "_scImpute_summary_edgeR_LRT"))
		load(paste0("cell_type_copy", i, "_pairs",  j, "_scImpute_summary_edgeR_QLF"))
		load(paste0("cell_type_copy", i, "_pairs",  j, "_scImpute_summary_scde"))
				
		cutoff1.100 <- sort(DESeq2.res1$p.adjust)[100]
		DE.id1.100 <- which(DESeq2.res1$p.adjust <= cutoff1.100)		
		cutoff2.100 <- sort(DESeq2.res2$p.adjust)[100]
		DE.id2.100 <- which(DESeq2.res2$p.adjust <= cutoff2.100)
		
		cutoff1.200 <- sort(DESeq2.res1$p.adjust)[200]
		DE.id1.200 <- which(DESeq2.res1$p.adjust <= cutoff1.200)		
		cutoff2.200 <- sort(DESeq2.res2$p.adjust)[200]
		DE.id2.200 <- which(DESeq2.res2$p.adjust <= cutoff2.200)
		
		cutoff1.500 <- sort(DESeq2.res1$p.adjust)[500]
		DE.id1.500 <- which(DESeq2.res1$p.adjust <= cutoff1.500)		
		cutoff2.500 <- sort(DESeq2.res2$p.adjust)[500]
		DE.id2.500 <- which(DESeq2.res2$p.adjust <= cutoff2.500)
		
		cutoff1.1000 <- sort(DESeq2.res1$p.adjust)[1000]
		DE.id1.1000 <- which(DESeq2.res1$p.adjust <= cutoff1.1000)		
		cutoff2.1000 <- sort(DESeq2.res2$p.adjust)[1000]
		DE.id2.1000 <- which(DESeq2.res2$p.adjust <= cutoff2.1000)
		
		res1 <- c(jaccard(DE.id1.100, DE.id2.100), jaccard(DE.id1.200, DE.id2.200), jaccard(DE.id1.500, DE.id2.500), jaccard(DE.id1.1000, DE.id2.1000))
		
		
		cutoff1.100 <- sort(edgeR.LRT.res1$p.adjust)[100]
		edgeR.LRT.id1.100 <- which(edgeR.LRT.res1$p.adjust <= cutoff1.100)		
		cutoff2.100 <- sort(edgeR.LRT.res2$p.adjust)[100]
		edgeR.LRT.id2.100 <- which(edgeR.LRT.res2$p.adjust <= cutoff2.100)
		
		cutoff1.200 <- sort(edgeR.LRT.res1$p.adjust)[200]
		edgeR.LRT.id1.200 <- which(edgeR.LRT.res1$p.adjust <= cutoff1.200)		
		cutoff2.200 <- sort(edgeR.LRT.res2$p.adjust)[200]
		edgeR.LRT.id2.200 <- which(edgeR.LRT.res2$p.adjust <= cutoff2.200)
		
		cutoff1.500 <- sort(edgeR.LRT.res1$p.adjust)[500]
		edgeR.LRT.id1.500 <- which(edgeR.LRT.res1$p.adjust <= cutoff1.500)		
		cutoff2.500 <- sort(edgeR.LRT.res2$p.adjust)[500]
		edgeR.LRT.id2.500 <- which(edgeR.LRT.res2$p.adjust <= cutoff2.500)
		
		cutoff1.1000 <- sort(edgeR.LRT.res1$p.adjust)[1000]
		edgeR.LRT.id1.1000 <- which(edgeR.LRT.res1$p.adjust <= cutoff1.1000)		
		cutoff2.1000 <- sort(edgeR.LRT.res2$p.adjust)[1000]
		edgeR.LRT.id2.1000 <- which(edgeR.LRT.res2$p.adjust <= cutoff2.1000)
		
		res2 <- c(jaccard(edgeR.LRT.id1.100, edgeR.LRT.id2.100), jaccard(edgeR.LRT.id1.200, edgeR.LRT.id2.200), jaccard(edgeR.LRT.id1.500, edgeR.LRT.id2.500), jaccard(edgeR.LRT.id1.1000, edgeR.LRT.id2.1000))		

        cutoff1.100 <- sort(edgeR.QLF.res1$p.adjust)[100]
		edgeR.QLF.id1.100 <- which(edgeR.QLF.res1$p.adjust <= cutoff1.100)		
		cutoff2.100 <- sort(edgeR.QLF.res2$p.adjust)[100]
		edgeR.QLF.id2.100 <- which(edgeR.QLF.res2$p.adjust <= cutoff2.100)
		
		cutoff1.200 <- sort(edgeR.QLF.res1$p.adjust)[200]
		edgeR.QLF.id1.200 <- which(edgeR.QLF.res1$p.adjust <= cutoff1.200)		
		cutoff2.200 <- sort(edgeR.QLF.res2$p.adjust)[200]
		edgeR.QLF.id2.200 <- which(edgeR.QLF.res2$p.adjust <= cutoff2.200)
		
		cutoff1.500 <- sort(edgeR.QLF.res1$p.adjust)[500]
		edgeR.QLF.id1.500 <- which(edgeR.QLF.res1$p.adjust <= cutoff1.500)		
		cutoff2.500 <- sort(edgeR.QLF.res2$p.adjust)[500]
		edgeR.QLF.id2.500 <- which(edgeR.QLF.res2$p.adjust <= cutoff2.500)
		
		cutoff1.1000 <- sort(edgeR.QLF.res1$p.adjust)[1000]
		edgeR.QLF.id1.1000 <- which(edgeR.QLF.res1$p.adjust <= cutoff1.1000)		
		cutoff2.1000 <- sort(edgeR.QLF.res2$p.adjust)[1000]
		edgeR.QLF.id2.1000 <- which(edgeR.QLF.res2$p.adjust <= cutoff2.1000)
		
		res3 <- c(jaccard(edgeR.QLF.id1.100, edgeR.QLF.id2.100), jaccard(edgeR.QLF.id1.200, edgeR.QLF.id2.200), jaccard(edgeR.QLF.id1.500, edgeR.QLF.id2.500), jaccard(edgeR.QLF.id1.1000, edgeR.QLF.id2.1000))
		
		pvalue1 <- 2*(1-pnorm(abs(scde.res1$Z)))
		pvalue2 <- 2*(1-pnorm(abs(scde.res2$Z)))
		
		cutoff1.100 <- sort(pvalue1)[100]
		scde.id1.100 <- which(pvalue1 <= cutoff1.100)		
		cutoff2.100 <- sort(pvalue2)[100]
		scde.id2.100 <- which(pvalue2 <= cutoff2.100)
		
		cutoff1.200 <- sort(pvalue1)[200]
		scde.id1.200 <- which(pvalue1 <= cutoff1.200)		
		cutoff2.200 <- sort(pvalue2)[200]
		scde.id2.200 <- which(pvalue2 <= cutoff2.200)	
				
		cutoff1.500 <- sort(pvalue1)[500]
		scde.id1.500 <- which(pvalue1 <= cutoff1.500)		
		cutoff2.500 <- sort(pvalue2)[500]
		scde.id2.500 <- which(pvalue2 <= cutoff2.500)

       	cutoff1.1000 <- sort(pvalue1)[1000]
		scde.id1.1000 <- which(pvalue1 <= cutoff1.1000)		
		cutoff2.1000 <- sort(pvalue2)[1000]
		scde.id2.1000 <- which(pvalue2 <= cutoff2.1000)
		
		res4 <- c(jaccard(scde.id1.100, scde.id2.100), jaccard(scde.id1.200, scde.id2.200), jaccard(scde.id1.500, scde.id2.500), jaccard(scde.id1.1000, scde.id2.1000))
		
		DEsummary <- data.frame(Dataset = rep("scImpute", 4), Method = c("DESeq2", "edgeR_LRT", "edgeR_QLF", "SCDE"), rbind(res1, res2, res3, res4))
		colnames(DEsummary)[-c(1:2)] <- c("top 100", "top 200", "top 500", "top 1000")
		save(DEsummary, file = paste0("cell_type_copy", i, "_pairs",  j, "_scImpute_top_DEsummary"))
	
	}
}




###########################################################################
############# DE analysis summary on MAGIC imputed data ####################
############################################################################


load("GSE75748_sc_cell_type.rda")
xx <- gene.expression
aa <- gsub("\\_[^<>]*", "", colnames(xx))
celltype <- unique(aa)



pairs <- matrix(c("EC", "DEC", "HFF", "DEC", "NPC", "DEC", "TB", "DEC",
"H9", "EC", "TB", "H1", "H9", "HFF", "H9", "NPC", "H9", "TB", "H9", "H1",  
"EC", "H1", "EC", "TB", "EC", "NPC",  "HFF", "TB", "HFF", "NPC",
"TB",  "NPC", "H1", "DEC",  "H9", "DEC",  "EC", "HFF",  "NPC", "TB",
"H1", "NPC"), byrow = T, ncol = 2)

library(VIPER)

for(i in 1:10){
	for(j in 1:21){	
		

		load(paste0("cell_type_copy", i, "_pairs",  j, "_magic_summary_DESeq2"))
		load(paste0("cell_type_copy", i, "_pairs",  j, "_magic_summary_edgeR_LRT"))
		load(paste0("cell_type_copy", i, "_pairs",  j, "_magic_summary_edgeR_QLF"))
		load(paste0("cell_type_copy", i, "_pairs",  j, "_magic_summary_scde"))
				
		cutoff1.100 <- sort(DESeq2.res1$p.adjust)[100]
		DE.id1.100 <- which(DESeq2.res1$p.adjust <= cutoff1.100)		
		cutoff2.100 <- sort(DESeq2.res2$p.adjust)[100]
		DE.id2.100 <- which(DESeq2.res2$p.adjust <= cutoff2.100)
		
		cutoff1.200 <- sort(DESeq2.res1$p.adjust)[200]
		DE.id1.200 <- which(DESeq2.res1$p.adjust <= cutoff1.200)		
		cutoff2.200 <- sort(DESeq2.res2$p.adjust)[200]
		DE.id2.200 <- which(DESeq2.res2$p.adjust <= cutoff2.200)
		
		cutoff1.500 <- sort(DESeq2.res1$p.adjust)[500]
		DE.id1.500 <- which(DESeq2.res1$p.adjust <= cutoff1.500)		
		cutoff2.500 <- sort(DESeq2.res2$p.adjust)[500]
		DE.id2.500 <- which(DESeq2.res2$p.adjust <= cutoff2.500)
		
		cutoff1.1000 <- sort(DESeq2.res1$p.adjust)[1000]
		DE.id1.1000 <- which(DESeq2.res1$p.adjust <= cutoff1.1000)		
		cutoff2.1000 <- sort(DESeq2.res2$p.adjust)[1000]
		DE.id2.1000 <- which(DESeq2.res2$p.adjust <= cutoff2.1000)
		
		res1 <- c(jaccard(DE.id1.100, DE.id2.100), jaccard(DE.id1.200, DE.id2.200), jaccard(DE.id1.500, DE.id2.500), jaccard(DE.id1.1000, DE.id2.1000))
		
		
		cutoff1.100 <- sort(edgeR.LRT.res1$p.adjust)[100]
		edgeR.LRT.id1.100 <- which(edgeR.LRT.res1$p.adjust <= cutoff1.100)		
		cutoff2.100 <- sort(edgeR.LRT.res2$p.adjust)[100]
		edgeR.LRT.id2.100 <- which(edgeR.LRT.res2$p.adjust <= cutoff2.100)
		
		cutoff1.200 <- sort(edgeR.LRT.res1$p.adjust)[200]
		edgeR.LRT.id1.200 <- which(edgeR.LRT.res1$p.adjust <= cutoff1.200)		
		cutoff2.200 <- sort(edgeR.LRT.res2$p.adjust)[200]
		edgeR.LRT.id2.200 <- which(edgeR.LRT.res2$p.adjust <= cutoff2.200)
		
		cutoff1.500 <- sort(edgeR.LRT.res1$p.adjust)[500]
		edgeR.LRT.id1.500 <- which(edgeR.LRT.res1$p.adjust <= cutoff1.500)		
		cutoff2.500 <- sort(edgeR.LRT.res2$p.adjust)[500]
		edgeR.LRT.id2.500 <- which(edgeR.LRT.res2$p.adjust <= cutoff2.500)
		
		cutoff1.1000 <- sort(edgeR.LRT.res1$p.adjust)[1000]
		edgeR.LRT.id1.1000 <- which(edgeR.LRT.res1$p.adjust <= cutoff1.1000)		
		cutoff2.1000 <- sort(edgeR.LRT.res2$p.adjust)[1000]
		edgeR.LRT.id2.1000 <- which(edgeR.LRT.res2$p.adjust <= cutoff2.1000)
		
		res2 <- c(jaccard(edgeR.LRT.id1.100, edgeR.LRT.id2.100), jaccard(edgeR.LRT.id1.200, edgeR.LRT.id2.200), jaccard(edgeR.LRT.id1.500, edgeR.LRT.id2.500), jaccard(edgeR.LRT.id1.1000, edgeR.LRT.id2.1000))		

        cutoff1.100 <- sort(edgeR.QLF.res1$p.adjust)[100]
		edgeR.QLF.id1.100 <- which(edgeR.QLF.res1$p.adjust <= cutoff1.100)		
		cutoff2.100 <- sort(edgeR.QLF.res2$p.adjust)[100]
		edgeR.QLF.id2.100 <- which(edgeR.QLF.res2$p.adjust <= cutoff2.100)
		
		cutoff1.200 <- sort(edgeR.QLF.res1$p.adjust)[200]
		edgeR.QLF.id1.200 <- which(edgeR.QLF.res1$p.adjust <= cutoff1.200)		
		cutoff2.200 <- sort(edgeR.QLF.res2$p.adjust)[200]
		edgeR.QLF.id2.200 <- which(edgeR.QLF.res2$p.adjust <= cutoff2.200)
		
		cutoff1.500 <- sort(edgeR.QLF.res1$p.adjust)[500]
		edgeR.QLF.id1.500 <- which(edgeR.QLF.res1$p.adjust <= cutoff1.500)		
		cutoff2.500 <- sort(edgeR.QLF.res2$p.adjust)[500]
		edgeR.QLF.id2.500 <- which(edgeR.QLF.res2$p.adjust <= cutoff2.500)
		
		cutoff1.1000 <- sort(edgeR.QLF.res1$p.adjust)[1000]
		edgeR.QLF.id1.1000 <- which(edgeR.QLF.res1$p.adjust <= cutoff1.1000)		
		cutoff2.1000 <- sort(edgeR.QLF.res2$p.adjust)[1000]
		edgeR.QLF.id2.1000 <- which(edgeR.QLF.res2$p.adjust <= cutoff2.1000)
		
		res3 <- c(jaccard(edgeR.QLF.id1.100, edgeR.QLF.id2.100), jaccard(edgeR.QLF.id1.200, edgeR.QLF.id2.200), jaccard(edgeR.QLF.id1.500, edgeR.QLF.id2.500), jaccard(edgeR.QLF.id1.1000, edgeR.QLF.id2.1000))
		
		pvalue1 <- 2*(1-pnorm(abs(scde.res1$Z)))
		pvalue2 <- 2*(1-pnorm(abs(scde.res2$Z)))
		
		cutoff1.100 <- sort(pvalue1)[100]
		scde.id1.100 <- which(pvalue1 <= cutoff1.100)		
		cutoff2.100 <- sort(pvalue2)[100]
		scde.id2.100 <- which(pvalue2 <= cutoff2.100)
		
		cutoff1.200 <- sort(pvalue1)[200]
		scde.id1.200 <- which(pvalue1 <= cutoff1.200)		
		cutoff2.200 <- sort(pvalue2)[200]
		scde.id2.200 <- which(pvalue2 <= cutoff2.200)	
				
		cutoff1.500 <- sort(pvalue1)[500]
		scde.id1.500 <- which(pvalue1 <= cutoff1.500)		
		cutoff2.500 <- sort(pvalue2)[500]
		scde.id2.500 <- which(pvalue2 <= cutoff2.500)

       	cutoff1.1000 <- sort(pvalue1)[1000]
		scde.id1.1000 <- which(pvalue1 <= cutoff1.1000)		
		cutoff2.1000 <- sort(pvalue2)[1000]
		scde.id2.1000 <- which(pvalue2 <= cutoff2.1000)
		
		res4 <- c(jaccard(scde.id1.100, scde.id2.100), jaccard(scde.id1.200, scde.id2.200), jaccard(scde.id1.500, scde.id2.500), jaccard(scde.id1.1000, scde.id2.1000))
		
		DEsummary <- data.frame(Dataset = rep("MAGIC", 4), Method = c("DESeq2", "edgeR_LRT", "edgeR_QLF", "SCDE"), rbind(res1, res2, res3, res4))
		colnames(DEsummary)[-c(1:2)] <- c("top 100", "top 200", "top 500", "top 1000")
		save(DEsummary, file = paste0("cell_type_copy", i, "_pairs",  j, "_magic_top_DEsummary"))
	
	}
}




############################################################# 
############# DE analysis summary for all methods ###########
############################################################# 

for(i in 1:10){
	for(j in 1:21){	
		allsummary <- NULL
		load(paste0("cell_type_copy", i, "_pairs",  j, "_magic_top_DEsummary"))
		allsummary <- rbind(allsummary, DEsummary)
		load(paste0("cell_type_copy", i, "_pairs",  j, "_scImpute_top_DEsummary"))
		allsummary <- rbind(allsummary, DEsummary)
		load(paste0("cell_type_copy", i, "_pairs",  j, "_drImpute_top_DEsummary"))
		allsummary <- rbind(allsummary, DEsummary)
		load(paste0("cell_type_copy", i, "_pairs",  j, "_el_top_DEsummary"))
		allsummary <- rbind(allsummary, DEsummary)
		load(paste0("cell_type_copy", i, "_pairs",  j, "_SAVER_top_DEsummary"))
		allsummary <- rbind(allsummary, DEsummary)
		load(paste0("cell_type_copy", i, "_pairs",  j, "_expression_top_DEsummary"))
		allsummary <- rbind(allsummary, DEsummary)
		load(paste0("cell_type_copy", i, "_pairs",  j, "_lasso_top_DEsummary"))
		allsummary <- rbind(allsummary, DEsummary)
		save(allsummary, file = paste0("cell_type_copy", i, "_pairs",  j, "_allmethods_topgenes_DEsummary"))
		
	}
}	

for(j in 1:21){	
	res <- NULL
	for(i in 1:10){
		load(paste0("cell_type_copy", i, "_pairs",  j, "_allmethods_topgenes_DEsummary"))
		res <- rbind(res, data.frame(allsummary, i))
	}
	save(res, file = paste0("cell_type_pairs",  j, "_allmethods_topgenes_DEsummary"))
}


