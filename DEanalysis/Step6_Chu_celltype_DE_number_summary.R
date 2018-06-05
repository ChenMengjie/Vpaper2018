

load("GSE75748_sc_cell_type.rda")
xx <- gene.expression
aa <- gsub("\\_[^<>]*", "", colnames(xx))
celltype <- unique(aa)


pairs <- matrix(c("H1", "DEC", 
"H9", "DEC", 
"EC", "HFF", 
"NPC", "TB",
"H1", "NPC"), byrow = T, ncol = 2)

library(scPoissonMix)

for(i in 1:10){
	for(j in 1:5){	
			
		load(paste0("cell_type_copy", i, "_pairs",  j, "_expression_summary_DESeq2"))
		load(paste0("cell_type_copy", i, "_pairs",  j, "_expression_summary_edgeR_LRT"))
		load(paste0("cell_type_copy", i, "_pairs",  j, "_expression_summary_edgeR_QLF"))
		load(paste0("cell_type_copy", i, "_pairs",  j, "_expression_summary_scde"))

		DE.id1 <- which(DESeq2.res1$p.adjust < 0.01)
		DE.id2 <- which(DESeq2.res2$p.adjust < 0.01)
		res1 <- c(length(DE.id1), length(DE.id2), length(intersect(DE.id1, DE.id2)))
		
		edgeR.LRT.id1 <- which(edgeR.LRT.res1$p.adjust < 0.01)
		edgeR.LRT.id2 <- which(edgeR.LRT.res2$p.adjust < 0.01)
		res2 <- c(length(edgeR.LRT.id1), length(edgeR.LRT.id2), length(intersect(edgeR.LRT.id1, edgeR.LRT.id2)))
		
				
		edgeR.QLF.id1 <- which(edgeR.QLF.res1$p.adjust < 0.01)
		edgeR.QLF.id2 <- which(edgeR.QLF.res2$p.adjust < 0.01)
		res3 <- c(length(edgeR.QLF.id1), length(edgeR.QLF.id2), length(intersect(edgeR.QLF.id1, edgeR.QLF.id2)))
		
		pvalue1 <- 2*(1-pnorm(abs(scde.res1$Z)))
		scde.id1 <- which(pvalue1 < 0.01)
		pvalue2 <- 2*(1-pnorm(abs(scde.res2$Z)))
		scde.id2 <- which(pvalue2 < 0.01)
		res4 <- c(length(scde.id1), length(scde.id2), length(intersect(scde.id1, scde.id2)))
		
		DEsummary <- data.frame(Dataset = rep("Before Imputation", 4), Method = c("DESeq2", "edgeR-LRT", "edgeR-QLF", "SCDE"), rbind(res1, res2, res3, res4))
		colnames(DEsummary)[-c(1:2)] <- c("Sig1", "Sig2", "Overlap")
		save(DEsummary, file = paste0("cell_type_copy", i, "_pairs",  j, "_expression_DEsummary"))
	
	}
}


load("GSE75748_sc_cell_type.rda")
xx <- gene.expression
aa <- gsub("\\_[^<>]*", "", colnames(xx))
celltype <- unique(aa)


pairs <- matrix(c("H1", "DEC", 
"H9", "DEC", 
"EC", "HFF", 
"NPC", "TB",
"H1", "NPC"), byrow = T, ncol = 2)

library(scPoissonMix)

for(i in 1:10){
	for(j in 1:5){	
			
		load(paste0("cell_type_copy", i, "_pairs",  j, "_drImpute_summary_DESeq2"))
		load(paste0("cell_type_copy", i, "_pairs",  j, "_drImpute_summary_edgeR_LRT"))
		load(paste0("cell_type_copy", i, "_pairs",  j, "_drImpute_summary_edgeR_QLF"))
		load(paste0("cell_type_copy", i, "_pairs",  j, "_drImpute_summary_scde"))

		DE.id1 <- which(DESeq2.res1$p.adjust < 0.01)
		DE.id2 <- which(DESeq2.res2$p.adjust < 0.01)
		res1 <- c(length(DE.id1), length(DE.id2), length(intersect(DE.id1, DE.id2)))
		
		edgeR.LRT.id1 <- which(edgeR.LRT.res1$p.adjust < 0.01)
		edgeR.LRT.id2 <- which(edgeR.LRT.res2$p.adjust < 0.01)
		res2 <- c(length(edgeR.LRT.id1), length(edgeR.LRT.id2), length(intersect(edgeR.LRT.id1, edgeR.LRT.id2)))
		
				
		edgeR.QLF.id1 <- which(edgeR.QLF.res1$p.adjust < 0.01)
		edgeR.QLF.id2 <- which(edgeR.QLF.res2$p.adjust < 0.01)
		res3 <- c(length(edgeR.QLF.id1), length(edgeR.QLF.id2), length(intersect(edgeR.QLF.id1, edgeR.QLF.id2)))
		
		pvalue1 <- 2*(1-pnorm(abs(scde.res1$Z)))
		scde.id1 <- which(pvalue1 < 0.01)
		pvalue2 <- 2*(1-pnorm(abs(scde.res2$Z)))
		scde.id2 <- which(pvalue2 < 0.01)
		res4 <- c(length(scde.id1), length(scde.id2), length(intersect(scde.id1, scde.id2)))
		
		DEsummary <- data.frame(Dataset = rep("DrImpute", 4), Method = c("DESeq2", "edgeR-LRT", "edgeR-QLF", "SCDE"), rbind(res1, res2, res3, res4))		
		colnames(DEsummary)[-c(1:2)] <- c("Sig1", "Sig2", "Overlap")
		save(DEsummary, file = paste0("cell_type_copy", i, "_pairs",  j, "_drImpute_DEsummary"))
	
	}
}







load("GSE75748_sc_cell_type.rda")
xx <- gene.expression
aa <- gsub("\\_[^<>]*", "", colnames(xx))
celltype <- unique(aa)


pairs <- matrix(c("H1", "DEC", 
"H9", "DEC", 
"EC", "HFF", 
"NPC", "TB",
"H1", "NPC"), byrow = T, ncol = 2)

library(scPoissonMix)

for(i in 1:10){
	for(j in 1:5){	
		
			
		load(paste0("cell_type_copy", i, "_pairs",  j, "_SAVER_summary_DESeq2"))
		load(paste0("cell_type_copy", i, "_pairs",  j, "_SAVER_summary_edgeR_LRT"))
		load(paste0("cell_type_copy", i, "_pairs",  j, "_SAVER_summary_edgeR_QLF"))
		load(paste0("cell_type_copy", i, "_pairs",  j, "_SAVER_summary_scde"))

		DE.id1 <- which(DESeq2.res1$p.adjust < 0.01)
		DE.id2 <- which(DESeq2.res2$p.adjust < 0.01)
		res1 <- c(length(DE.id1), length(DE.id2), length(intersect(DE.id1, DE.id2)))
		
		edgeR.LRT.id1 <- which(edgeR.LRT.res1$p.adjust < 0.01)
		edgeR.LRT.id2 <- which(edgeR.LRT.res2$p.adjust < 0.01)
		res2 <- c(length(edgeR.LRT.id1), length(edgeR.LRT.id2), length(intersect(edgeR.LRT.id1, edgeR.LRT.id2)))
		
				
		edgeR.QLF.id1 <- which(edgeR.QLF.res1$p.adjust < 0.01)
		edgeR.QLF.id2 <- which(edgeR.QLF.res2$p.adjust < 0.01)
		res3 <- c(length(edgeR.QLF.id1), length(edgeR.QLF.id2), length(intersect(edgeR.QLF.id1, edgeR.QLF.id2)))
		
		pvalue1 <- 2*(1-pnorm(abs(scde.res1$Z)))
		scde.id1 <- which(pvalue1 < 0.01)
		pvalue2 <- 2*(1-pnorm(abs(scde.res2$Z)))
		scde.id2 <- which(pvalue2 < 0.01)
		res4 <- c(length(scde.id1), length(scde.id2), length(intersect(scde.id1, scde.id2)))
		
		DEsummary <- data.frame(Dataset = rep("SAVER", 4), Method = c("DESeq2", "edgeR-LRT", "edgeR-QLF", "SCDE"), rbind(res1, res2, res3, res4))						
		colnames(DEsummary)[-c(1:2)] <- c("Sig1", "Sig2", "Overlap")
		save(DEsummary, file = paste0("cell_type_copy", i, "_pairs",  j, "_SAVER_DEsummary"))
	
	}
}





load("GSE75748_sc_cell_type.rda")
xx <- gene.expression
aa <- gsub("\\_[^<>]*", "", colnames(xx))
celltype <- unique(aa)


pairs <- matrix(c("H1", "DEC", 
"H9", "DEC", 
"EC", "HFF", 
"NPC", "TB",
"H1", "NPC"), byrow = T, ncol = 2)

library(scPoissonMix)

for(i in 1:10){
	for(j in 1:5){	
		
		load(paste0("cell_type_copy", i, "_pairs",  j, "_scImpute_summary_DESeq2"))
		load(paste0("cell_type_copy", i, "_pairs",  j, "_scImpute_summary_edgeR_LRT"))
		load(paste0("cell_type_copy", i, "_pairs",  j, "_scImpute_summary_edgeR_QLF"))
		load(paste0("cell_type_copy", i, "_pairs",  j, "_scImpute_summary_scde"))

		DE.id1 <- which(DESeq2.res1$p.adjust < 0.01)
		DE.id2 <- which(DESeq2.res2$p.adjust < 0.01)
		res1 <- c(length(DE.id1), length(DE.id2), length(intersect(DE.id1, DE.id2)))
		
		edgeR.LRT.id1 <- which(edgeR.LRT.res1$p.adjust < 0.01)
		edgeR.LRT.id2 <- which(edgeR.LRT.res2$p.adjust < 0.01)
		res2 <- c(length(edgeR.LRT.id1), length(edgeR.LRT.id2), length(intersect(edgeR.LRT.id1, edgeR.LRT.id2)))
		
				
		edgeR.QLF.id1 <- which(edgeR.QLF.res1$p.adjust < 0.01)
		edgeR.QLF.id2 <- which(edgeR.QLF.res2$p.adjust < 0.01)
		res3 <- c(length(edgeR.QLF.id1), length(edgeR.QLF.id2), length(intersect(edgeR.QLF.id1, edgeR.QLF.id2)))
		
		pvalue1 <- 2*(1-pnorm(abs(scde.res1$Z)))
		scde.id1 <- which(pvalue1 < 0.01)
		pvalue2 <- 2*(1-pnorm(abs(scde.res2$Z)))
		scde.id2 <- which(pvalue2 < 0.01)
		res4 <- c(length(scde.id1), length(scde.id2), length(intersect(scde.id1, scde.id2)))
		
		DEsummary <- data.frame(Dataset = rep("scImpute", 4), Method = c("DESeq2", "edgeR-LRT", "edgeR-QLF", "SCDE"), rbind(res1, res2, res3, res4))					
		colnames(DEsummary)[-c(1:2)] <- c("Sig1", "Sig2", "Overlap")
		save(DEsummary, file = paste0("cell_type_copy", i, "_pairs",  j, "_scImpute_DEsummary"))
	
	}
}






load("GSE75748_sc_cell_type.rda")
xx <- gene.expression
aa <- gsub("\\_[^<>]*", "", colnames(xx))
celltype <- unique(aa)


pairs <- matrix(c("H1", "DEC", 
"H9", "DEC", 
"EC", "HFF", 
"NPC", "TB",
"H1", "NPC"), byrow = T, ncol = 2)

library(scPoissonMix)

for(i in 1:10){
	for(j in 1:5){	

		
		load(paste0("cell_type_copy", i, "_pairs",  j, "_magic_summary_DESeq2"))
		load(paste0("cell_type_copy", i, "_pairs",  j, "_magic_summary_edgeR_LRT"))
		load(paste0("cell_type_copy", i, "_pairs",  j, "_magic_summary_edgeR_QLF"))
		load(paste0("cell_type_copy", i, "_pairs",  j, "_magic_summary_scde"))

		DE.id1 <- which(DESeq2.res1$p.adjust < 0.01)
		DE.id2 <- which(DESeq2.res2$p.adjust < 0.01)
		res1 <- c(length(DE.id1), length(DE.id2), length(intersect(DE.id1, DE.id2)))
		
		edgeR.LRT.id1 <- which(edgeR.LRT.res1$p.adjust < 0.01)
		edgeR.LRT.id2 <- which(edgeR.LRT.res2$p.adjust < 0.01)
		res2 <- c(length(edgeR.LRT.id1), length(edgeR.LRT.id2), length(intersect(edgeR.LRT.id1, edgeR.LRT.id2)))
		
				
		edgeR.QLF.id1 <- which(edgeR.QLF.res1$p.adjust < 0.01)
		edgeR.QLF.id2 <- which(edgeR.QLF.res2$p.adjust < 0.01)
		res3 <- c(length(edgeR.QLF.id1), length(edgeR.QLF.id2), length(intersect(edgeR.QLF.id1, edgeR.QLF.id2)))
		
		pvalue1 <- 2*(1-pnorm(abs(scde.res1$Z)))
		scde.id1 <- which(pvalue1 < 0.01)
		pvalue2 <- 2*(1-pnorm(abs(scde.res2$Z)))
		scde.id2 <- which(pvalue2 < 0.01)
		res4 <- c(length(scde.id1), length(scde.id2), length(intersect(scde.id1, scde.id2)))
		
		DEsummary <- data.frame(Dataset = rep("MAGIC", 4), Method = c("DESeq2", "edgeR-LRT", "edgeR-QLF", "SCDE"), rbind(res1, res2, res3, res4))	
		colnames(DEsummary)[-c(1:2)] <- c("Sig1", "Sig2", "Overlap")
		save(DEsummary, file = paste0("cell_type_copy", i, "_pairs",  j, "_magic_DEsummary"))
	
	}
}




load("GSE75748_sc_cell_type.rda")
xx <- gene.expression
aa <- gsub("\\_[^<>]*", "", colnames(xx))
celltype <- unique(aa)


pairs <- matrix(c("H1", "DEC", 
"H9", "DEC", 
"EC", "HFF", 
"NPC", "TB",
"H1", "NPC"), byrow = T, ncol = 2)

library(scPoissonMix)

for(i in 1:10){
	for(j in 1:5){	
		
		load(paste0("cell_type_copy", i, "_pairs",  j, "_ee_summary_DESeq2"))
		load(paste0("cell_type_copy", i, "_pairs",  j, "_ee_summary_edgeR_LRT"))
		load(paste0("cell_type_copy", i, "_pairs",  j, "_ee_summary_edgeR_QLF"))
		load(paste0("cell_type_copy", i, "_pairs",  j, "_ee_summary_scde"))

		DE.id1 <- which(DESeq2.res1$p.adjust < 0.01)
		DE.id2 <- which(DESeq2.res2$p.adjust < 0.01)
		res1 <- c(length(DE.id1), length(DE.id2), length(intersect(DE.id1, DE.id2)))
		
		edgeR.LRT.id1 <- which(edgeR.LRT.res1$p.adjust < 0.01)
		edgeR.LRT.id2 <- which(edgeR.LRT.res2$p.adjust < 0.01)
		res2 <- c(length(edgeR.LRT.id1), length(edgeR.LRT.id2), length(intersect(edgeR.LRT.id1, edgeR.LRT.id2)))
		
				
		edgeR.QLF.id1 <- which(edgeR.QLF.res1$p.adjust < 0.01)
		edgeR.QLF.id2 <- which(edgeR.QLF.res2$p.adjust < 0.01)
		res3 <- c(length(edgeR.QLF.id1), length(edgeR.QLF.id2), length(intersect(edgeR.QLF.id1, edgeR.QLF.id2)))
		
		pvalue1 <- 2*(1-pnorm(abs(scde.res1$Z)))
		scde.id1 <- which(pvalue1 < 0.01)
		pvalue2 <- 2*(1-pnorm(abs(scde.res2$Z)))
		scde.id2 <- which(pvalue2 < 0.01)
		res4 <- c(length(scde.id1), length(scde.id2), length(intersect(scde.id1, scde.id2)))
		
		DEsummary <- data.frame(Dataset = rep("VIPER-Lasso", 4), Method = c("DESeq2", "edgeR-LRT", "edgeR-QLF", "SCDE"), rbind(res1, res2, res3, res4))					
		colnames(DEsummary)[-c(1:2)] <- c("Sig1", "Sig2", "Overlap")
		save(DEsummary, file = paste0("cell_type_copy", i, "_pairs",  j, "_ee_DEsummary"))
	
	}
}





load("GSE75748_sc_cell_type.rda")
xx <- gene.expression
aa <- gsub("\\_[^<>]*", "", colnames(xx))
celltype <- unique(aa)


pairs <- matrix(c("H1", "DEC", 
"H9", "DEC", 
"EC", "HFF", 
"NPC", "TB",
"H1", "NPC"), byrow = T, ncol = 2)

library(scPoissonMix)

for(i in 1:10){
	for(j in 1:5){	
		
		
		load(paste0("cell_type_copy", i, "_pairs",  j, "_el_ee_summary_DESeq2"))
		load(paste0("cell_type_copy", i, "_pairs",  j, "_el_ee_summary_edgeR_LRT"))
		load(paste0("cell_type_copy", i, "_pairs",  j, "_el_ee_summary_edgeR_QLF"))
		load(paste0("cell_type_copy", i, "_pairs",  j, "_el_ee_summary_scde"))

		DE.id1 <- which(DESeq2.res1$p.adjust < 0.01)
		DE.id2 <- which(DESeq2.res2$p.adjust < 0.01)
		res1 <- c(length(DE.id1), length(DE.id2), length(intersect(DE.id1, DE.id2)))
		
		edgeR.LRT.id1 <- which(edgeR.LRT.res1$p.adjust < 0.01)
		edgeR.LRT.id2 <- which(edgeR.LRT.res2$p.adjust < 0.01)
		res2 <- c(length(edgeR.LRT.id1), length(edgeR.LRT.id2), length(intersect(edgeR.LRT.id1, edgeR.LRT.id2)))
		
				
		edgeR.QLF.id1 <- which(edgeR.QLF.res1$p.adjust < 0.01)
		edgeR.QLF.id2 <- which(edgeR.QLF.res2$p.adjust < 0.01)
		res3 <- c(length(edgeR.QLF.id1), length(edgeR.QLF.id2), length(intersect(edgeR.QLF.id1, edgeR.QLF.id2)))
		
		pvalue1 <- 2*(1-pnorm(abs(scde.res1$Z)))
		scde.id1 <- which(pvalue1 < 0.01)
		pvalue2 <- 2*(1-pnorm(abs(scde.res2$Z)))
		scde.id2 <- which(pvalue2 < 0.01)
		res4 <- c(length(scde.id1), length(scde.id2), length(intersect(scde.id1, scde.id2)))
		
		DEsummary <- data.frame(Dataset = rep("VIPER-Elastic net", 4), Method = c("DESeq2", "edgeR-LRT", "edgeR-QLF", "SCDE"), rbind(res1, res2, res3, res4))					
		colnames(DEsummary)[-c(1:2)] <- c("Sig1", "Sig2", "Overlap")
		save(DEsummary, file = paste0("cell_type_copy", i, "_pairs",  j, "_el_ee_DEsummary"))
	
	}
}






for(i in 1:5){
	for(j in 1:5){	
		allsummary <- NULL
		load(paste0("cell_type_copy", i, "_pairs",  j, "_magic_DEsummary"))
		allsummary <- rbind(allsummary, DEsummary)
		load(paste0("cell_type_copy", i, "_pairs",  j, "_scImpute_DEsummary"))
		allsummary <- rbind(allsummary, DEsummary)
		load(paste0("cell_type_copy", i, "_pairs",  j, "_drImpute_DEsummary"))
		allsummary <- rbind(allsummary, DEsummary)
		load(paste0("cell_type_copy", i, "_pairs",  j, "_SAVER_DEsummary"))
		allsummary <- rbind(allsummary, DEsummary)
		load(paste0("cell_type_copy", i, "_pairs",  j, "_expression_DEsummary"))
		allsummary <- rbind(allsummary, DEsummary)
		load(paste0("cell_type_copy", i, "_pairs",  j, "_ee_DEsummary"))
		allsummary <- rbind(allsummary, DEsummary)
		load(paste0("cell_type_copy", i, "_pairs",  j, "_el_ee_DEsummary"))
		allsummary <- rbind(allsummary, DEsummary)
		
		allsummary$Prop <- apply(allsummary[, -c(1:2)], 1, function(x){
			y <- as.numeric(x)
			round(x[3]/max(x[1], x[2]), 2)
		})
		save(allsummary, file = paste0("cell_type_copy", i, "_pairs",  j, "_allmethods_DEsummary"))
		
	}
}	




###################
for(j in 1:5){	
	res <- NULL
	for(i in 1:5){
		load(paste0("cell_type_copy", i, "_pairs",  j, "_allmethods_DEsummary"))
		res <- rbind(res, data.frame(allsummary, i))
	}
	save(res, file = paste0("cell_type_pairs",  j, "_allmethods_DEsummary"))
}
























