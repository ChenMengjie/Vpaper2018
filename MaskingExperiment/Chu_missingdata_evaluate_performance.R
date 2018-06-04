all.res <- NULL

for(i in 1:10){
	 load(paste0("Chu_time_course_missing2_copy", i, "_imputation.summary"))
	 log.topredict <- log(topredict+1)
	 log.scimpute.pred <- log(scimpute.pred+1)
	 log.lnpX.pred <- lnpX.pred
	 log.imputed.pred <- log(exp(imputed.pred)-0.1+1)
	 log.imputed.el.pred <- log(exp(imputed.el.pred)-0.1+1)
	 log.magic.pred <- log(magic.pred+1)
	 log.saver.pred <- log(saver.pred+1)
 
	 aa <- c(median(abs(log.topredict - log.scimpute.pred)),
	 	median(abs(log.topredict - log.lnpX.pred)),
	 	median(abs(log.topredict - log.imputed.pred)),
	 	median(abs(log.topredict - log.magic.pred)),
	 	median(abs(log.topredict - log.saver.pred)), 
	 	median(abs(log.topredict - log.imputed.el.pred)))
	 
	 all.res <- rbind(all.res, aa)
}

all.res2 <- NULL

for(i in 1:10){
	 load(paste0("Chu_cell_type_missing2_copy", i, "_imputation.summary"))
	 log.topredict <- log(topredict+1)
	 log.scimpute.pred <- log(scimpute.pred+1)
	 log.lnpX.pred <- lnpX.pred
	 log.imputed.pred <- log(exp(imputed.pred)-0.1+1)
	 log.imputed.el.pred <- log(exp(imputed.el.pred)-0.1+1)
	 log.magic.pred <- log(magic.pred+1)
	 log.saver.pred <- log(saver.pred+1)
 
	 aa <- c(median(abs(log.topredict - log.scimpute.pred)),
	 	median(abs(log.topredict - log.lnpX.pred)),
	 	median(abs(log.topredict - log.imputed.pred)),
	 	median(abs(log.topredict - log.magic.pred)),
	 	median(abs(log.topredict - log.saver.pred)), 
	 	median(abs(log.topredict - log.imputed.el.pred)))
		 
	 all.res2 <- rbind(all.res2, aa)
}

absolute.loss.median.summary <- data.frame(dataset = c(rep("TimeCourse", nrow(all.res)), rep("CellType", nrow(all.res))), 
           scenario = c(rep("missing2", nrow(all.res)), rep("missing2", nrow(all.res))),
           scimpute = c(all.res[, 1], all.res2[, 1]),
           lnpX = c(all.res[, 2], all.res2[, 2]),
           saver = c(all.res[, 5], all.res2[, 5]),
           imputed.lasso = c(all.res[, 3], all.res2[, 3]),       
           magic = c(all.res[, 4], all.res2[, 4]),
           imputed.elasticnet = c(all.res[, 6], all.res2[, 6]))
           
        
all.res <- NULL

for(i in 1:10){
	 load(paste0("Chu_time_course_missing2_copy", i, "_imputation.summary"))
	 log.topredict <- log(topredict+1)
	 log.scimpute.pred <- log(scimpute.pred+1)
	 log.lnpX.pred <- lnpX.pred
	 log.imputed.pred <- log(exp(imputed.pred)-0.1+1)
	 log.imputed.el.pred <- log(exp(imputed.el.pred)-0.1+1)
	 log.magic.pred <- log(magic.pred+1)
	 log.saver.pred <- log(saver.pred+1)

	 aa <- c(median((log.topredict - log.scimpute.pred)^2),
	 	median((log.topredict - log.lnpX.pred)^2),
	 	median((log.topredict - log.imputed.pred)^2),
	 	median((log.topredict - log.magic.pred)^2),
	 	median((log.topredict - log.saver.pred)^2), 
	 	median((log.topredict - log.imputed.el.pred)^2))
	 
	 all.res <- rbind(all.res, aa)
}

all.res2 <- NULL

for(i in 1:10){
	 load(paste0("Chu_cell_type_missing2_copy", i, "_imputation.summary"))
	 log.topredict <- log(topredict+1)
	 log.scimpute.pred <- log(scimpute.pred+1)
	 log.lnpX.pred <- lnpX.pred
	 log.imputed.pred <- log(exp(imputed.pred)-0.1+1)
	 log.imputed.el.pred <- log(exp(imputed.el.pred)-0.1+1)
	 log.magic.pred <- log(magic.pred+1)
	 log.saver.pred <- log(saver.pred+1)

	 aa <- c(median((log.topredict - log.scimpute.pred)^2),
	 	median((log.topredict - log.lnpX.pred)^2),
	 	median((log.topredict - log.imputed.pred)^2),
	 	median((log.topredict - log.magic.pred)^2),
	 	median((log.topredict - log.saver.pred)^2), 
	 	median((log.topredict - log.imputed.el.pred)^2))

	 	
	 all.res2 <- rbind(all.res2, aa)
}

square.loss.median.summary <- data.frame(dataset = c(rep("TimeCourse", nrow(all.res)), rep("CellType", nrow(all.res))), 
           scenario = c(rep("missing2", nrow(all.res)), rep("missing2", nrow(all.res))),
           scimpute = c(all.res[, 1], all.res2[, 1]),
           lnpX = c(all.res[, 2], all.res2[, 2]),
           saver = c(all.res[, 5], all.res2[, 5]),
           imputed.lasso = c(all.res[, 3], all.res2[, 3]),       
           magic = c(all.res[, 4], all.res2[, 4]),
           imputed.elasticnet = c(all.res[, 6], all.res2[, 6]))
                               
           
all.res <- NULL

for(i in 1:10){
	 load(paste0("Chu_cell_type_missing2_copy", i, "_imputation.summary"))
	 log.topredict <- log(topredict+1)
	 log.scimpute.pred <- log(scimpute.pred+1)
	 log.lnpX.pred <- lnpX.pred
	 log.imputed.pred <- log(exp(imputed.pred)-0.1+1)
	 log.imputed.el.pred <- log(exp(imputed.el.pred)-0.1+1)
	 log.magic.pred <- log(magic.pred+1)
	 log.saver.pred <- log(saver.pred+1)
	
	 aa <- c(cor(log.topredict, log.scimpute.pred),
	 	cor(log.topredict[!is.na(log.lnpX.pred)], log.lnpX.pred[!is.na(log.lnpX.pred)]),
	 	cor(log.topredict, log.imputed.pred),
	 	cor(log.topredict, log.magic.pred),
	 	cor(log.topredict, log.saver.pred),
	 	cor(log.topredict, log.imputed.el.pred))
	 	
	 all.res <- rbind(all.res, aa)
}


all.res2 <- NULL

for(i in 1:10){
	 load(paste0("Chu_cell_type_missing2_copy", i, "_imputation.summary"))
	 log.topredict <- log(topredict+1)
	 log.scimpute.pred <- log(scimpute.pred+1)
	 log.lnpX.pred <- lnpX.pred
	 log.imputed.pred <- log(exp(imputed.pred)-0.1+1)
	 log.imputed.el.pred <- log(exp(imputed.el.pred)-0.1+1)
	 log.magic.pred <- log(magic.pred+1)
	 log.saver.pred <- log(saver.pred+1)
	
	 aa <- c(cor(log.topredict, log.scimpute.pred),
	 	cor(log.topredict[!is.na(log.lnpX.pred)], log.lnpX.pred[!is.na(log.lnpX.pred)]),
	 	cor(log.topredict, log.imputed.pred),
	 	cor(log.topredict, log.magic.pred),
	 	cor(log.topredict, log.saver.pred),
	 	cor(log.topredict, log.imputed.el.pred))
	 	
	 all.res2 <- rbind(all.res2, aa)
}

correlation.median.summary <- data.frame(dataset = c(rep("TimeCourse", nrow(all.res)), rep("CellType", nrow(all.res))), 
           scenario = c(rep("missing2", nrow(all.res)), rep("missing2", nrow(all.res))),
           scimpute = c(all.res[, 1], all.res2[, 1]),
           lnpX = c(all.res[, 2], all.res2[, 2]),
           saver = c(all.res[, 5], all.res2[, 5]),
           imputed.lasso = c(all.res[, 3], all.res2[, 3]),       
           magic = c(all.res[, 4], all.res2[, 4]),
           imputed.elasticnet = c(all.res[, 6], all.res2[, 6]))
                   
           
save(correlation.median.summary, square.loss.median.summary, absolute.loss.median.summary, file = "Chu_imputation_with_missing2_summary")



all.res <- NULL

for(i in 1:10){
	 load(paste0("Chu_time_course_missing5_copy", i, "_imputation.summary"))
	 log.topredict <- log(topredict+1)
	 log.scimpute.pred <- log(scimpute.pred+1)
	 log.lnpX.pred <- lnpX.pred
	 log.imputed.pred <- log(exp(imputed.pred)-0.1+1)
	 log.imputed.el.pred <- log(exp(imputed.el.pred)-0.1+1)
	 log.magic.pred <- log(magic.pred+1)
	 log.saver.pred <- log(saver.pred+1)
 
	 aa <- c(median(abs(log.topredict - log.scimpute.pred)),
	 	median(abs(log.topredict - log.lnpX.pred)),
	 	median(abs(log.topredict - log.imputed.pred)),
	 	median(abs(log.topredict - log.magic.pred)),
	 	median(abs(log.topredict - log.saver.pred)), 
	 	median(abs(log.topredict - log.imputed.el.pred)))
	 
	 all.res <- rbind(all.res, aa)
}

all.res2 <- NULL

for(i in 1:10){
	 load(paste0("Chu_cell_type_missing5_copy", i, "_imputation.summary"))
	 log.topredict <- log(topredict+1)
	 log.scimpute.pred <- log(scimpute.pred+1)
	 log.lnpX.pred <- lnpX.pred
	 log.imputed.pred <- log(exp(imputed.pred)-0.1+1)
	 log.imputed.el.pred <- log(exp(imputed.el.pred)-0.1+1)
	 log.magic.pred <- log(magic.pred+1)
	 log.saver.pred <- log(saver.pred+1)
 
	 aa <- c(median(abs(log.topredict - log.scimpute.pred)),
	 	median(abs(log.topredict - log.lnpX.pred)),
	 	median(abs(log.topredict - log.imputed.pred)),
	 	median(abs(log.topredict - log.magic.pred)),
	 	median(abs(log.topredict - log.saver.pred)), 
	 	median(abs(log.topredict - log.imputed.el.pred)))
		 
	 all.res2 <- rbind(all.res2, aa)
}

absolute.loss.median.summary <- data.frame(dataset = c(rep("TimeCourse", nrow(all.res)), rep("CellType", nrow(all.res))), 
           scenario = c(rep("missing5", nrow(all.res)), rep("missing5", nrow(all.res))),
           scimpute = c(all.res[, 1], all.res2[, 1]),
           lnpX = c(all.res[, 2], all.res2[, 2]),
           saver = c(all.res[, 5], all.res2[, 5]),
           imputed.lasso = c(all.res[, 3], all.res2[, 3]),       
           magic = c(all.res[, 4], all.res2[, 4]),
           imputed.elasticnet = c(all.res[, 6], all.res2[, 6]))
           
        
all.res <- NULL

for(i in 1:10){
	 load(paste0("Chu_time_course_missing5_copy", i, "_imputation.summary"))
	 log.topredict <- log(topredict+1)
	 log.scimpute.pred <- log(scimpute.pred+1)
	 log.lnpX.pred <- lnpX.pred
	 log.imputed.pred <- log(exp(imputed.pred)-0.1+1)
	 log.imputed.el.pred <- log(exp(imputed.el.pred)-0.1+1)
	 log.magic.pred <- log(magic.pred+1)
	 log.saver.pred <- log(saver.pred+1)

	 aa <- c(median((log.topredict - log.scimpute.pred)^2),
	 	median((log.topredict - log.lnpX.pred)^2),
	 	median((log.topredict - log.imputed.pred)^2),
	 	median((log.topredict - log.magic.pred)^2),
	 	median((log.topredict - log.saver.pred)^2), 
	 	median((log.topredict - log.imputed.el.pred)^2))
	 
	 all.res <- rbind(all.res, aa)
}

all.res2 <- NULL

for(i in 1:10){
	 load(paste0("Chu_cell_type_missing5_copy", i, "_imputation.summary"))
	 log.topredict <- log(topredict+1)
	 log.scimpute.pred <- log(scimpute.pred+1)
	 log.lnpX.pred <- lnpX.pred
	 log.imputed.pred <- log(exp(imputed.pred)-0.1+1)
	 log.imputed.el.pred <- log(exp(imputed.el.pred)-0.1+1)
	 log.magic.pred <- log(magic.pred+1)
	 log.saver.pred <- log(saver.pred+1)

	 aa <- c(median((log.topredict - log.scimpute.pred)^2),
	 	median((log.topredict - log.lnpX.pred)^2),
	 	median((log.topredict - log.imputed.pred)^2),
	 	median((log.topredict - log.magic.pred)^2),
	 	median((log.topredict - log.saver.pred)^2), 
	 	median((log.topredict - log.imputed.el.pred)^2))

	 	
	 all.res2 <- rbind(all.res2, aa)
}

square.loss.median.summary <- data.frame(dataset = c(rep("TimeCourse", nrow(all.res)), rep("CellType", nrow(all.res))), 
           scenario = c(rep("missing5", nrow(all.res)), rep("missing5", nrow(all.res))),
           scimpute = c(all.res[, 1], all.res2[, 1]),
           lnpX = c(all.res[, 2], all.res2[, 2]),
           saver = c(all.res[, 5], all.res2[, 5]),
           imputed.lasso = c(all.res[, 3], all.res2[, 3]),       
           magic = c(all.res[, 4], all.res2[, 4]),
           imputed.elasticnet = c(all.res[, 6], all.res2[, 6]))
                               
           
all.res <- NULL

for(i in 1:10){
	 load(paste0("Chu_cell_type_missing5_copy", i, "_imputation.summary"))
	 log.topredict <- log(topredict+1)
	 log.scimpute.pred <- log(scimpute.pred+1)
	 log.lnpX.pred <- lnpX.pred
	 log.imputed.pred <- log(exp(imputed.pred)-0.1+1)
	 log.imputed.el.pred <- log(exp(imputed.el.pred)-0.1+1)
	 log.magic.pred <- log(magic.pred+1)
	 log.saver.pred <- log(saver.pred+1)
	
	 aa <- c(cor(log.topredict, log.scimpute.pred),
	 	cor(log.topredict[!is.na(log.lnpX.pred)], log.lnpX.pred[!is.na(log.lnpX.pred)]),
	 	cor(log.topredict, log.imputed.pred),
	 	cor(log.topredict, log.magic.pred),
	 	cor(log.topredict, log.saver.pred),
	 	cor(log.topredict, log.imputed.el.pred))
	 	
	 all.res <- rbind(all.res, aa)
}


all.res2 <- NULL

for(i in 1:10){
	 load(paste0("Chu_cell_type_missing5_copy", i, "_imputation.summary"))
	 log.topredict <- log(topredict+1)
	 log.scimpute.pred <- log(scimpute.pred+1)
	 log.lnpX.pred <- lnpX.pred
	 log.imputed.pred <- log(exp(imputed.pred)-0.1+1)
	 log.imputed.el.pred <- log(exp(imputed.el.pred)-0.1+1)
	 log.magic.pred <- log(magic.pred+1)
	 log.saver.pred <- log(saver.pred+1)
	
	 aa <- c(cor(log.topredict, log.scimpute.pred),
	 	cor(log.topredict[!is.na(log.lnpX.pred)], log.lnpX.pred[!is.na(log.lnpX.pred)]),
	 	cor(log.topredict, log.imputed.pred),
	 	cor(log.topredict, log.magic.pred),
	 	cor(log.topredict, log.saver.pred),
	 	cor(log.topredict, log.imputed.el.pred))
	 	
	 all.res2 <- rbind(all.res2, aa)
}

correlation.median.summary <- data.frame(dataset = c(rep("TimeCourse", nrow(all.res)), rep("CellType", nrow(all.res))), 
           scenario = c(rep("missing5", nrow(all.res)), rep("missing5", nrow(all.res))),
           scimpute = c(all.res[, 1], all.res2[, 1]),
           lnpX = c(all.res[, 2], all.res2[, 2]),
           saver = c(all.res[, 5], all.res2[, 5]),
           imputed.lasso = c(all.res[, 3], all.res2[, 3]),       
           magic = c(all.res[, 4], all.res2[, 4]),
           imputed.elasticnet = c(all.res[, 6], all.res2[, 6]))
                   
           
save(correlation.median.summary, square.loss.median.summary, absolute.loss.median.summary, file = "Chu_imputation_with_missing5_summary")


all.res <- NULL

for(i in 1:10){
	 load(paste0("Chu_time_course_missing10_copy", i, "_imputation.summary"))
	 log.topredict <- log(topredict+1)
	 log.scimpute.pred <- log(scimpute.pred+1)
	 log.lnpX.pred <- lnpX.pred
	 log.imputed.pred <- log(exp(imputed.pred)-0.1+1)
	 log.imputed.el.pred <- log(exp(imputed.el.pred)-0.1+1)
	 log.magic.pred <- log(magic.pred+1)
	 log.saver.pred <- log(saver.pred+1)
 
	 aa <- c(median(abs(log.topredict - log.scimpute.pred)),
	 	median(abs(log.topredict - log.lnpX.pred)),
	 	median(abs(log.topredict - log.imputed.pred)),
	 	median(abs(log.topredict - log.magic.pred)),
	 	median(abs(log.topredict - log.saver.pred)), 
	 	median(abs(log.topredict - log.imputed.el.pred)))
	 
	 all.res <- rbind(all.res, aa)
}

all.res2 <- NULL

for(i in 1:10){
	 load(paste0("Chu_cell_type_missing10_copy", i, "_imputation.summary"))
	 log.topredict <- log(topredict+1)
	 log.scimpute.pred <- log(scimpute.pred+1)
	 log.lnpX.pred <- lnpX.pred
	 log.imputed.pred <- log(exp(imputed.pred)-0.1+1)
	 log.imputed.el.pred <- log(exp(imputed.el.pred)-0.1+1)
	 log.magic.pred <- log(magic.pred+1)
	 log.saver.pred <- log(saver.pred+1)
 
	 aa <- c(median(abs(log.topredict - log.scimpute.pred)),
	 	median(abs(log.topredict - log.lnpX.pred)),
	 	median(abs(log.topredict - log.imputed.pred)),
	 	median(abs(log.topredict - log.magic.pred)),
	 	median(abs(log.topredict - log.saver.pred)), 
	 	median(abs(log.topredict - log.imputed.el.pred)))
		 
	 all.res2 <- rbind(all.res2, aa)
}

absolute.loss.median.summary <- data.frame(dataset = c(rep("TimeCourse", nrow(all.res)), rep("CellType", nrow(all.res))), 
           scenario = c(rep("missing10", nrow(all.res)), rep("missing10", nrow(all.res))),
           scimpute = c(all.res[, 1], all.res2[, 1]),
           lnpX = c(all.res[, 2], all.res2[, 2]),
           saver = c(all.res[, 5], all.res2[, 5]),
           imputed.lasso = c(all.res[, 3], all.res2[, 3]),       
           magic = c(all.res[, 4], all.res2[, 4]),
           imputed.elasticnet = c(all.res[, 6], all.res2[, 6]))
           
        
all.res <- NULL

for(i in 1:10){
	 load(paste0("Chu_time_course_missing10_copy", i, "_imputation.summary"))
	 log.topredict <- log(topredict+1)
	 log.scimpute.pred <- log(scimpute.pred+1)
	 log.lnpX.pred <- lnpX.pred
	 log.imputed.pred <- log(exp(imputed.pred)-0.1+1)
	 log.imputed.el.pred <- log(exp(imputed.el.pred)-0.1+1)
	 log.magic.pred <- log(magic.pred+1)
	 log.saver.pred <- log(saver.pred+1)

	 aa <- c(median((log.topredict - log.scimpute.pred)^2),
	 	median((log.topredict - log.lnpX.pred)^2),
	 	median((log.topredict - log.imputed.pred)^2),
	 	median((log.topredict - log.magic.pred)^2),
	 	median((log.topredict - log.saver.pred)^2), 
	 	median((log.topredict - log.imputed.el.pred)^2))
	 
	 all.res <- rbind(all.res, aa)
}

all.res2 <- NULL

for(i in 1:10){
	 load(paste0("Chu_cell_type_missing10_copy", i, "_imputation.summary"))
	 log.topredict <- log(topredict+1)
	 log.scimpute.pred <- log(scimpute.pred+1)
	 log.lnpX.pred <- lnpX.pred
	 log.imputed.pred <- log(exp(imputed.pred)-0.1+1)
	 log.imputed.el.pred <- log(exp(imputed.el.pred)-0.1+1)
	 log.magic.pred <- log(magic.pred+1)
	 log.saver.pred <- log(saver.pred+1)

	 aa <- c(median((log.topredict - log.scimpute.pred)^2),
	 	median((log.topredict - log.lnpX.pred)^2),
	 	median((log.topredict - log.imputed.pred)^2),
	 	median((log.topredict - log.magic.pred)^2),
	 	median((log.topredict - log.saver.pred)^2), 
	 	median((log.topredict - log.imputed.el.pred)^2))

	 	
	 all.res2 <- rbind(all.res2, aa)
}

square.loss.median.summary <- data.frame(dataset = c(rep("TimeCourse", nrow(all.res)), rep("CellType", nrow(all.res))), 
           scenario = c(rep("missing10", nrow(all.res)), rep("missing10", nrow(all.res))),
           scimpute = c(all.res[, 1], all.res2[, 1]),
           lnpX = c(all.res[, 2], all.res2[, 2]),
           saver = c(all.res[, 5], all.res2[, 5]),
           imputed.lasso = c(all.res[, 3], all.res2[, 3]),       
           magic = c(all.res[, 4], all.res2[, 4]),
           imputed.elasticnet = c(all.res[, 6], all.res2[, 6]))
                               
           
all.res <- NULL

for(i in 1:10){
	 load(paste0("Chu_cell_type_missing10_copy", i, "_imputation.summary"))
	 log.topredict <- log(topredict+1)
	 log.scimpute.pred <- log(scimpute.pred+1)
	 log.lnpX.pred <- lnpX.pred
	 log.imputed.pred <- log(exp(imputed.pred)-0.1+1)
	 log.imputed.el.pred <- log(exp(imputed.el.pred)-0.1+1)
	 log.magic.pred <- log(magic.pred+1)
	 log.saver.pred <- log(saver.pred+1)
	
	 aa <- c(cor(log.topredict, log.scimpute.pred),
	 	cor(log.topredict[!is.na(log.lnpX.pred)], log.lnpX.pred[!is.na(log.lnpX.pred)]),
	 	cor(log.topredict, log.imputed.pred),
	 	cor(log.topredict, log.magic.pred),
	 	cor(log.topredict, log.saver.pred),
	 	cor(log.topredict, log.imputed.el.pred))
	 	
	 all.res <- rbind(all.res, aa)
}


all.res2 <- NULL

for(i in 1:10){
	 load(paste0("Chu_cell_type_missing10_copy", i, "_imputation.summary"))
	 log.topredict <- log(topredict+1)
	 log.scimpute.pred <- log(scimpute.pred+1)
	 log.lnpX.pred <- lnpX.pred
	 log.imputed.pred <- log(exp(imputed.pred)-0.1+1)
	 log.imputed.el.pred <- log(exp(imputed.el.pred)-0.1+1)
	 log.magic.pred <- log(magic.pred+1)
	 log.saver.pred <- log(saver.pred+1)
	
	 aa <- c(cor(log.topredict, log.scimpute.pred),
	 	cor(log.topredict[!is.na(log.lnpX.pred)], log.lnpX.pred[!is.na(log.lnpX.pred)]),
	 	cor(log.topredict, log.imputed.pred),
	 	cor(log.topredict, log.magic.pred),
	 	cor(log.topredict, log.saver.pred),
	 	cor(log.topredict, log.imputed.el.pred))
	 	
	 all.res2 <- rbind(all.res2, aa)
}

correlation.median.summary <- data.frame(dataset = c(rep("TimeCourse", nrow(all.res)), rep("CellType", nrow(all.res))), 
           scenario = c(rep("missing10", nrow(all.res)), rep("missing10", nrow(all.res))),
           scimpute = c(all.res[, 1], all.res2[, 1]),
           lnpX = c(all.res[, 2], all.res2[, 2]),
           saver = c(all.res[, 5], all.res2[, 5]),
           imputed.lasso = c(all.res[, 3], all.res2[, 3]),       
           magic = c(all.res[, 4], all.res2[, 4]),
           imputed.elasticnet = c(all.res[, 6], all.res2[, 6]))
                   
           
save(correlation.median.summary, square.loss.median.summary, absolute.loss.median.summary, file = "Chu_imputation_with_missing10_summary")


