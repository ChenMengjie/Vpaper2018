load("Shalek_filtered.rda")
for(i in 1:10){
	set.seed(i*500+80224523)
	len <- length(gene.expression[gene.expression!=0])
	missing <- sample(c(0, 1), len, replace=T, prob=c(0.02, 0.98))
	topredict <- gene.expression[gene.expression!=0][missing==0]
	
	scimpute <- read.csv(paste0("scImpute/Shalek/missing2_copy", i, "scimpute_count.csv"), as.is = T)
	scimpute.imp <- scimpute[, -1]
	scimpute.pred <- scimpute.imp[gene.expression!=0][missing==0]

	load(paste0("Shalek_drImpute_missing2_copy", i))	
	lnpX.pred <- lnpX_imp[gene.expression!=0][missing==0]
	
	load(paste0("Shalek_SAVER_missing2_copy", i))
	saver.pred <- saver_imp$estimate[gene.expression!=0][missing==0]
	
	load(paste0("Shalek_lasso_missing2_copy", i))
	imputed.pred <- imp.ee$imputed[gene.expression!=0][missing==0]
	
	load(paste0("Shalek_el_missing2_copy", i))
	imputed.el.pred <- imp.el.ee$imputed[gene.expression!=0][missing==0]
	
	magic <- read.csv(paste0("Shalek_missing2_copy_t", i, "_magic.csv"))
	magic.imp <- t(magic[, -1])
	magic.pred <- magic.imp[gene.expression!=0][missing==0]

	save(topredict, scimpute.pred, lnpX.pred, saver.pred, imputed.pred, imputed.el.pred, magic.pred, file = paste0("Shalek_missing2_copy", i, "_imputation.summary"))
	
}



load("Shalek_filtered.rda")
for(i in 1:10){
	set.seed(i*500+80224523)
	len <- length(gene.expression[gene.expression!=0])
	missing <- sample(c(0, 1), len, replace=T, prob=c(0.05, 0.95))
	topredict <- gene.expression[gene.expression!=0][missing==0]
	
	scimpute <- read.csv(paste0("scImpute/Shalek/missing5_copy", i, "scimpute_count.csv"), as.is = T)
	scimpute.imp <- scimpute[, -1]
	scimpute.pred <- scimpute.imp[gene.expression!=0][missing==0]

	load(paste0("Shalek_drImpute_missing5_copy", i))	
	lnpX.pred <- lnpX_imp[gene.expression!=0][missing==0]
	
	load(paste0("Shalek_SAVER_missing5_copy", i))
	saver.pred <- saver_imp$estimate[gene.expression!=0][missing==0]
	
	load(paste0("Shalek_lasso_missing5_copy", i))
	imputed.pred <- imp.ee$imputed[gene.expression!=0][missing==0]
	
	load(paste0("Shalek_el_missing5_copy", i))
	imputed.el.pred <- imp.el.ee$imputed[gene.expression!=0][missing==0]
	
	magic <- read.csv(paste0("Shalek_missing5_copy_t", i, "_magic.csv"))
	magic.imp <- t(magic[, -1])
	magic.pred <- magic.imp[gene.expression!=0][missing==0]

	save(topredict, scimpute.pred, lnpX.pred, saver.pred, imputed.pred, imputed.el.pred, magic.pred, file = paste0("Shalek_missing5_copy", i, "_imputation.summary"))
	
}


load("Shalek_filtered.rda")
for(i in 1:10){
	set.seed(i*500+80224523)
	len <- length(gene.expression[gene.expression!=0])
	missing <- sample(c(0, 1), len, replace=T, prob=c(0.1, 0.9))
	topredict <- gene.expression[gene.expression!=0][missing==0]
	
	scimpute <- read.csv(paste0("scImpute/Shalek/missing10_copy", i, "scimpute_count.csv"), as.is = T)
	scimpute.imp <- scimpute[, -1]
	scimpute.pred <- scimpute.imp[gene.expression!=0][missing==0]

	load(paste0("Shalek_drImpute_missing10_copy", i))	
	lnpX.pred <- lnpX_imp[gene.expression!=0][missing==0]
	
	load(paste0("Shalek_SAVER_missing10_copy", i))
	saver.pred <- saver_imp$estimate[gene.expression!=0][missing==0]
	
	load(paste0("Shalek_lasso_missing10_copy", i))
	imputed.pred <- imp.ee$imputed[gene.expression!=0][missing==0]
	
	load(paste0("Shalek_el_missing10_copy", i))
	imputed.el.pred <- imp.el.ee$imputed[gene.expression!=0][missing==0]
	
	magic <- read.csv(paste0("Shalek_missing10_copy_t", i, "_magic.csv"))
	magic.imp <- t(magic[, -1])
	magic.pred <- magic.imp[gene.expression!=0][missing==0]

	save(topredict, scimpute.pred, lnpX.pred, saver.pred, imputed.pred, imputed.el.pred, magic.pred, file = paste0("Shalek_missing10_copy", i, "_imputation.summary"))
	
}

