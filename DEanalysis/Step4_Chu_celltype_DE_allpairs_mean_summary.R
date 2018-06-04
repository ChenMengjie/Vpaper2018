pairs <- matrix(c("EC", "DEC", "HFF", "DEC", "NPC", "DEC", "TB", "DEC",
"H9", "EC", "TB", "H1", "H9", "HFF", "H9", "NPC", "H9", "TB", "H9", "H1",  
"EC", "H1", "EC", "TB", "EC", "NPC",  "HFF", "TB", "HFF", "NPC",
"TB",  "NPC", "H1", "DEC",  "H9", "DEC",  "EC", "HFF",  "NPC", "TB",
"H1", "NPC"), byrow = T, ncol = 2)

total.pairs <- NULL
for(j in 1:21){	
	total <- NULL
	allsummary <- NULL
	for(i in 1:10){
		load(paste0("cell_type_copy", i, "_pairs",  j, "_magic_top_DEsummary"))
		allsummary <- rbind(allsummary, DEsummary)
	}	
	aa1 <- cbind(tapply(allsummary[, 3], allsummary[, 2], mean), tapply(allsummary[, 3], allsummary[, 2], sd))
	aa2 <- cbind(tapply(allsummary[, 4], allsummary[, 2], mean), tapply(allsummary[, 4], allsummary[, 2], sd))
	aa3 <- cbind(tapply(allsummary[, 5], allsummary[, 2], mean), tapply(allsummary[, 5], allsummary[, 2], sd))
    aa4 <- cbind(tapply(allsummary[, 6], allsummary[, 2], mean), tapply(allsummary[, 6], allsummary[, 2], sd))
	asummary <- data.frame(Dataset = allsummary[1, 1], Method = rep(rownames(aa1), 4), Mean = c(aa1[, 1], aa2[, 1], aa3[, 1], aa4[, 1]), SE = c(aa1[, 2], aa2[, 2], aa3[, 2], aa4[, 2]), Num = c(rep(100, 4), rep(200, 4), rep(500, 4), rep(1000, 4)))
	total <- rbind(total, asummary)
	
   allsummary <- NULL
   for(i in 1:10){
		load(paste0("cell_type_copy", i, "_pairs",  j, "_scImpute_top_DEsummary"))
		allsummary <- rbind(allsummary, DEsummary)
	}	
	aa1 <- cbind(tapply(allsummary[, 3], allsummary[, 2], mean), tapply(allsummary[, 3], allsummary[, 2], sd))
	aa2 <- cbind(tapply(allsummary[, 4], allsummary[, 2], mean), tapply(allsummary[, 4], allsummary[, 2], sd))
	aa3 <- cbind(tapply(allsummary[, 5], allsummary[, 2], mean), tapply(allsummary[, 5], allsummary[, 2], sd))
    aa4 <- cbind(tapply(allsummary[, 6], allsummary[, 2], mean), tapply(allsummary[, 6], allsummary[, 2], sd))
	asummary <- data.frame(Dataset = allsummary[1, 1], Method = rep(rownames(aa1), 4), Mean = c(aa1[, 1], aa2[, 1], aa3[, 1], aa4[, 1]), SE = c(aa1[, 2], aa2[, 2], aa3[, 2], aa4[, 2]), Num = c(rep(100, 4), rep(200, 4), rep(500, 4), rep(1000, 4)))
   total <- rbind(total, asummary)
   
   allsummary <- NULL
   for(i in 1:10){
		load(paste0("cell_type_copy", i, "_pairs",  j, "_drImpute_top_DEsummary"))
		allsummary <- rbind(allsummary, DEsummary)
	}	
	aa1 <- cbind(tapply(allsummary[, 3], allsummary[, 2], mean), tapply(allsummary[, 3], allsummary[, 2], sd))
	aa2 <- cbind(tapply(allsummary[, 4], allsummary[, 2], mean), tapply(allsummary[, 4], allsummary[, 2], sd))
	aa3 <- cbind(tapply(allsummary[, 5], allsummary[, 2], mean), tapply(allsummary[, 5], allsummary[, 2], sd))
    aa4 <- cbind(tapply(allsummary[, 6], allsummary[, 2], mean), tapply(allsummary[, 6], allsummary[, 2], sd))
	asummary <- data.frame(Dataset = allsummary[1, 1], Method = rep(rownames(aa1), 4), Mean = c(aa1[, 1], aa2[, 1], aa3[, 1], aa4[, 1]), SE = c(aa1[, 2], aa2[, 2], aa3[, 2], aa4[, 2]), Num = c(rep(100, 4), rep(200, 4), rep(500, 4), rep(1000, 4)))    
	total <- rbind(total, asummary)
    
   allsummary <- NULL
   for(i in 1:10){
		load(paste0("cell_type_copy", i, "_pairs",  j, "_el_top_DEsummary"))
		allsummary <- rbind(allsummary, DEsummary)
	}	
	aa1 <- cbind(tapply(allsummary[, 3], allsummary[, 2], mean), tapply(allsummary[, 3], allsummary[, 2], sd))
	aa2 <- cbind(tapply(allsummary[, 4], allsummary[, 2], mean), tapply(allsummary[, 4], allsummary[, 2], sd))
	aa3 <- cbind(tapply(allsummary[, 5], allsummary[, 2], mean), tapply(allsummary[, 5], allsummary[, 2], sd))
    aa4 <- cbind(tapply(allsummary[, 6], allsummary[, 2], mean), tapply(allsummary[, 6], allsummary[, 2], sd))
	asummary <- data.frame(Dataset = allsummary[1, 1], Method = rep(rownames(aa1), 4), Mean = c(aa1[, 1], aa2[, 1], aa3[, 1], aa4[, 1]), SE = c(aa1[, 2], aa2[, 2], aa3[, 2], aa4[, 2]), Num = c(rep(100, 4), rep(200, 4), rep(500, 4), rep(1000, 4)))  
	total <- rbind(total, asummary)
      
   allsummary <- NULL
   for(i in 1:10){
		load(paste0("cell_type_copy", i, "_pairs",  j, "_SAVER_top_DEsummary"))
		allsummary <- rbind(allsummary, DEsummary)
	}	
	aa1 <- cbind(tapply(allsummary[, 3], allsummary[, 2], mean), tapply(allsummary[, 3], allsummary[, 2], sd))
	aa2 <- cbind(tapply(allsummary[, 4], allsummary[, 2], mean), tapply(allsummary[, 4], allsummary[, 2], sd))
	aa3 <- cbind(tapply(allsummary[, 5], allsummary[, 2], mean), tapply(allsummary[, 5], allsummary[, 2], sd))
    aa4 <- cbind(tapply(allsummary[, 6], allsummary[, 2], mean), tapply(allsummary[, 6], allsummary[, 2], sd))
	asummary <- data.frame(Dataset = allsummary[1, 1], Method = rep(rownames(aa1), 4), Mean = c(aa1[, 1], aa2[, 1], aa3[, 1], aa4[, 1]), SE = c(aa1[, 2], aa2[, 2], aa3[, 2], aa4[, 2]), Num = c(rep(100, 4), rep(200, 4), rep(500, 4), rep(1000, 4))) 
	total <- rbind(total, asummary)
    
    
   allsummary <- NULL
   for(i in 1:10){
		load(paste0("cell_type_copy", i, "_pairs",  j, "_expression_top_DEsummary"))
		allsummary <- rbind(allsummary, DEsummary)
	}	
	aa1 <- cbind(tapply(allsummary[, 3], allsummary[, 2], mean), tapply(allsummary[, 3], allsummary[, 2], sd))
	aa2 <- cbind(tapply(allsummary[, 4], allsummary[, 2], mean), tapply(allsummary[, 4], allsummary[, 2], sd))
	aa3 <- cbind(tapply(allsummary[, 5], allsummary[, 2], mean), tapply(allsummary[, 5], allsummary[, 2], sd))
    aa4 <- cbind(tapply(allsummary[, 6], allsummary[, 2], mean), tapply(allsummary[, 6], allsummary[, 2], sd))
	asummary <- data.frame(Dataset = allsummary[1, 1], Method = rep(rownames(aa1), 4), Mean = c(aa1[, 1], aa2[, 1], aa3[, 1], aa4[, 1]), SE = c(aa1[, 2], aa2[, 2], aa3[, 2], aa4[, 2]), Num = c(rep(100, 4), rep(200, 4), rep(500, 4), rep(1000, 4)))    
	total <- rbind(total, asummary)
    
   allsummary <- NULL
   for(i in 1:10){
		load(paste0("cell_type_copy", i, "_pairs",  j, "_lasso_top_DEsummary"))
		allsummary <- rbind(allsummary, DEsummary)
	}	
	aa1 <- cbind(tapply(allsummary[, 3], allsummary[, 2], mean), tapply(allsummary[, 3], allsummary[, 2], sd))
	aa2 <- cbind(tapply(allsummary[, 4], allsummary[, 2], mean), tapply(allsummary[, 4], allsummary[, 2], sd))
	aa3 <- cbind(tapply(allsummary[, 5], allsummary[, 2], mean), tapply(allsummary[, 5], allsummary[, 2], sd))
    aa4 <- cbind(tapply(allsummary[, 6], allsummary[, 2], mean), tapply(allsummary[, 6], allsummary[, 2], sd))
	asummary <- data.frame(Dataset = allsummary[1, 1], Method = rep(rownames(aa1), 4), Mean = c(aa1[, 1], aa2[, 1], aa3[, 1], aa4[, 1]), SE = c(aa1[, 2], aa2[, 2], aa3[, 2], aa4[, 2]), Num = c(rep(100, 4), rep(200, 4), rep(500, 4), rep(1000, 4)))
	total <- rbind(total, asummary)
	
	total.pairs <- rbind(total.pairs, data.frame(total, pair = paste0(pairs[j, 1], "-", pairs[j, 2])))
		
}
	
save(total.pairs, file = paste0("cell_type_pairs_allmethods_topgenes_mean_se_DEsummary"))


