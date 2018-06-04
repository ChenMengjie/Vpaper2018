load("GSE75748_sc_time_course.rda")
for(i in 1:10){
	set.seed(i*500+80224523)
	len <- length(gene.expression[gene.expression!=0])
	missing <- sample(c(0, 1), len, replace=T, prob=c(0.02, 0.98))
	topredict <- gene.expression[gene.expression!=0][missing==0]
	
	scimpute <- read.csv(paste0("scImpute/time_course/missing2_copy", i, "scimpute_count.csv"), as.is = T)
	scimpute.imp <- scimpute[, -1]
	scimpute.pred <- scimpute.imp[gene.expression!=0][missing==0]

	load(paste0("Chu_drImpute_time_course_missing2_copy", i))	
	lnpX.pred <- lnpX_imp[gene.expression!=0][missing==0]
	
	load(paste0("Chu_SAVER_time_course_missing2_copy", i))
	saver.pred <- saver_imp$estimate[gene.expression!=0][missing==0]
	
	load(paste0("Chu_lasso_time_course_missing2_copy", i))
	imputed.pred <- imp.ee$imputed[gene.expression!=0][missing==0]
	
	load(paste0("Chu_el_time_course_missing2_copy", i))
	imputed.el.pred <- imp.el.ee$imputed[gene.expression!=0][missing==0]

	magic <- read.csv(paste0("Chu_time_course_missing2_copy", i, "_magic.csv"))
	magic.imp <- t(magic[, -1])
	magic.pred <- magic.imp[gene.expression!=0][missing==0]
	
	save(topredict, scimpute.pred, lnpX.pred, saver.pred, imputed.pred, imputed.el.pred, magic.pred, file = paste0("Chu_time_course_missing2_copy", i, "_imputation.summary"))
	
}



load("GSE75748_sc_cell_type.rda")
for(i in 1:10){
	set.seed(i*500+80224523)
	len <- length(gene.expression[gene.expression!=0])
	missing <- sample(c(0, 1), len, replace=T, prob=c(0.02, 0.98))
	topredict <- gene.expression[gene.expression!=0][missing==0]
	
	scimpute <- read.csv(paste0("scImpute/cell_type/missing2_copy", i, "scimpute_count.csv"), as.is = T)
	scimpute.imp <- scimpute[, -1]
	scimpute.pred <- scimpute.imp[gene.expression!=0][missing==0]

	load(paste0("Chu_drImpute_cell_type_missing2_copy", i))	
	lnpX.pred <- lnpX_imp[gene.expression!=0][missing==0]
	
	load(paste0("Chu_SAVER_cell_type_missing2_copy", i))
	saver.pred <- saver_imp$estimate[gene.expression!=0][missing==0]
	
	load(paste0("Chu_lasso_cell_type_missing2_copy", i))
	imputed.pred <- imp.ee$imputed[gene.expression!=0][missing==0]
	
	load(paste0("Chu_el_cell_type_missing2_copy", i))
	imputed.el.pred <- imp.el.ee$imputed[gene.expression!=0][missing==0]
	
	magic <- read.csv(paste0("Chu_cell_type_missing2_copy", i, "_magic.csv"))
	magic.imp <- t(magic[, -1])
	magic.pred <- magic.imp[gene.expression!=0][missing==0]
	
	save(topredict, scimpute.pred, lnpX.pred, saver.pred, imputed.pred, imputed.el.pred, magic.pred, file = paste0("Chu_cell_type_missing2_copy", i, "_imputation.summary"))
	
}



load("GSE75748_sc_time_course.rda")
for(i in 1:10){
	set.seed(i*500+80224523)
	len <- length(gene.expression[gene.expression!=0])
	missing <- sample(c(0, 1), len, replace=T, prob=c(0.05, 0.95))
	topredict <- gene.expression[gene.expression!=0][missing==0]
	
	scimpute <- read.csv(paste0("scImpute/time_course/missing5_copy", i, "scimpute_count.csv"), as.is = T)
	scimpute.imp <- scimpute[, -1]
	scimpute.pred <- scimpute.imp[gene.expression!=0][missing==0]

	load(paste0("Chu_drImpute_time_course_missing5_copy", i))	
	lnpX.pred <- lnpX_imp[gene.expression!=0][missing==0]
	
	load(paste0("Chu_SAVER_time_course_missing5_copy", i))
	saver.pred <- saver_imp$estimate[gene.expression!=0][missing==0]
	
	load(paste0("Chu_lasso_time_course_missing5_copy", i))
	imputed.pred <- imp.ee$imputed[gene.expression!=0][missing==0]
	
	load(paste0("Chu_el_time_course_missing5_copy", i))
	imputed.el.pred <- imp.el.ee$imputed[gene.expression!=0][missing==0]

	magic <- read.csv(paste0("Chu_time_course_missing5_copy", i, "_magic.csv"))
	magic.imp <- t(magic[, -1])
	magic.pred <- magic.imp[gene.expression!=0][missing==0]
	
	save(topredict, scimpute.pred, lnpX.pred, saver.pred, imputed.pred, imputed.el.pred, magic.pred, file = paste0("Chu_time_course_missing5_copy", i, "_imputation.summary"))
	
}



load("GSE75748_sc_cell_type.rda")
for(i in 1:10){
	set.seed(i*500+80224523)
	len <- length(gene.expression[gene.expression!=0])
	missing <- sample(c(0, 1), len, replace=T, prob=c(0.05, 0.95))
	topredict <- gene.expression[gene.expression!=0][missing==0]
	
	scimpute <- read.csv(paste0("scImpute/cell_type/missing5_copy", i, "scimpute_count.csv"), as.is = T)
	scimpute.imp <- scimpute[, -1]
	scimpute.pred <- scimpute.imp[gene.expression!=0][missing==0]

	load(paste0("Chu_drImpute_cell_type_missing5_copy", i))	
	lnpX.pred <- lnpX_imp[gene.expression!=0][missing==0]
	
	load(paste0("Chu_SAVER_cell_type_missing5_copy", i))
	saver.pred <- saver_imp$estimate[gene.expression!=0][missing==0]
	
	load(paste0("Chu_lasso_cell_type_missing5_copy", i))
	imputed.pred <- imp.ee$imputed[gene.expression!=0][missing==0]
	
	load(paste0("Chu_el_cell_type_missing5_copy", i))
	imputed.el.pred <- imp.el.ee$imputed[gene.expression!=0][missing==0]
	
	magic <- read.csv(paste0("Chu_cell_type_missing5_copy", i, "_magic.csv"))
	magic.imp <- t(magic[, -1])
	magic.pred <- magic.imp[gene.expression!=0][missing==0]
	
	save(topredict, scimpute.pred, lnpX.pred, saver.pred, imputed.pred, imputed.el.pred, magic.pred, file = paste0("Chu_cell_type_missing5_copy", i, "_imputation.summary"))
	
}

	
	
load("GSE75748_sc_time_course.rda")
for(i in 1:10){
	set.seed(i*500+80224523)
	len <- length(gene.expression[gene.expression!=0])
	missing <- sample(c(0, 1), len, replace=T, prob=c(0.1, 0.9))
	topredict <- gene.expression[gene.expression!=0][missing==0]
	
	scimpute <- read.csv(paste0("scImpute/time_course/missing10_copy", i, "scimpute_count.csv"), as.is = T)
	scimpute.imp <- scimpute[, -1]
	scimpute.pred <- scimpute.imp[gene.expression!=0][missing==0]

	load(paste0("Chu_drImpute_time_course_missing10_copy", i))	
	lnpX.pred <- lnpX_imp[gene.expression!=0][missing==0]
	
	load(paste0("Chu_SAVER_time_course_missing10_copy", i))
	saver.pred <- saver_imp$estimate[gene.expression!=0][missing==0]
	
	load(paste0("Chu_lasso_time_course_missing10_copy", i))
	imputed.pred <- imp.ee$imputed[gene.expression!=0][missing==0]
	
	load(paste0("Chu_el_time_course_missing10_copy", i))
	imputed.el.pred <- imp.el.ee$imputed[gene.expression!=0][missing==0]

	magic <- read.csv(paste0("Chu_time_course_missing10_copy", i, "_magic.csv"))
	magic.imp <- t(magic[, -1])
	magic.pred <- magic.imp[gene.expression!=0][missing==0]
	
	save(topredict, scimpute.pred, lnpX.pred, saver.pred, imputed.pred, imputed.el.pred, magic.pred, file = paste0("Chu_time_course_missing10_copy", i, "_imputation.summary"))
	
}



load("GSE75748_sc_cell_type.rda")
for(i in 1:10){
	set.seed(i*500+80224523)
	len <- length(gene.expression[gene.expression!=0])
	missing <- sample(c(0, 1), len, replace=T, prob=c(0.1, 0.9))
	topredict <- gene.expression[gene.expression!=0][missing==0]
	
	scimpute <- read.csv(paste0("scImpute/cell_type/missing10_copy", i, "scimpute_count.csv"), as.is = T)
	scimpute.imp <- scimpute[, -1]
	scimpute.pred <- scimpute.imp[gene.expression!=0][missing==0]

	load(paste0("Chu_drImpute_cell_type_missing10_copy", i))	
	lnpX.pred <- lnpX_imp[gene.expression!=0][missing==0]
	
	load(paste0("Chu_SAVER_cell_type_missing10_copy", i))
	saver.pred <- saver_imp$estimate[gene.expression!=0][missing==0]
	
	load(paste0("Chu_lasso_cell_type_missing10_copy", i))
	imputed.pred <- imp.ee$imputed[gene.expression!=0][missing==0]
	
	load(paste0("Chu_el_cell_type_missing10_copy", i))
	imputed.el.pred <- imp.el.ee$imputed[gene.expression!=0][missing==0]
	
	magic <- read.csv(paste0("Chu_cell_type_missing10_copy", i, "_magic.csv"))
	magic.imp <- t(magic[, -1])
	magic.pred <- magic.imp[gene.expression!=0][missing==0]
	
	save(topredict, scimpute.pred, lnpX.pred, saver.pred, imputed.pred, imputed.el.pred, magic.pred, file = paste0("Chu_cell_type_missing10_copy", i, "_imputation.summary"))
	
}

		
	