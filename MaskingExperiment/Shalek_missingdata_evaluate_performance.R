all.res <- NULL

for(i in 1:10){
	 load(paste0("Shalek_missing2_copy", i, "_imputation.summary"))
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

absolute.loss.median.summary <- data.frame(dataset = rep("Shalek", nrow(all.res)), 
           scenario = rep("missing2", nrow(all.res)),
           scimpute = all.res[, 1],
           lnpX = all.res[, 2], 
           saver = all.res[, 5],
           imputed.lasso = all.res[, 3],  
           magic = all.res[, 5],
           imputed.elasticnet  = all.res[, 6])           
  
  

all.res <- NULL

for(i in 1:10){
	 load(paste0("Shalek_missing2_copy", i, "_imputation.summary"))
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

square.loss.median.summary <- data.frame(dataset = rep("Shalek", nrow(all.res)), 
           scenario = rep("missing2", nrow(all.res)),
           scimpute = all.res[, 1],
           lnpX = all.res[, 2], 
           saver = all.res[, 5],
           imputed.lasso = all.res[, 3],  
           magic = all.res[, 5],
           imputed.elasticnet  = all.res[, 6])           
  
  
  

all.res <- NULL

for(i in 1:10){
	 load(paste0("Shalek_missing2_copy", i, "_imputation.summary"))
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

correlation.median.summary <- data.frame(dataset = rep("Shalek", nrow(all.res)), 
           scenario = rep("missing2", nrow(all.res)),
           scimpute = all.res[, 1],
           lnpX = all.res[, 2], 
           saver = all.res[, 5],
           imputed.lasso = all.res[, 3],  
           magic = all.res[, 5],
           imputed.elasticnet  = all.res[, 6])           
  
   
save(correlation.median.summary, square.loss.median.summary, absolute.loss.median.summary, file = "Shalek_imputation_with_missing2_summary")



all.res <- NULL

for(i in 1:10){
	 load(paste0("Shalek_missing5_copy", i, "_imputation.summary"))
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

absolute.loss.median.summary <- data.frame(dataset = rep("Shalek", nrow(all.res)), 
           scenario = rep("missing5", nrow(all.res)),
           scimpute = all.res[, 1],
           lnpX = all.res[, 2], 
           saver = all.res[, 5],
           imputed.lasso = all.res[, 3],  
           magic = all.res[, 5],
           imputed.elasticnet  = all.res[, 6])           
  
  

all.res <- NULL

for(i in 1:10){
	 load(paste0("Shalek_missing5_copy", i, "_imputation.summary"))
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

square.loss.median.summary <- data.frame(dataset = rep("Shalek", nrow(all.res)), 
           scenario = rep("missing5", nrow(all.res)),
           scimpute = all.res[, 1],
           lnpX = all.res[, 2], 
           saver = all.res[, 5],
           imputed.lasso = all.res[, 3],  
           magic = all.res[, 5],
           imputed.elasticnet  = all.res[, 6])           
  
  
  

all.res <- NULL

for(i in 1:10){
	 load(paste0("Shalek_missing5_copy", i, "_imputation.summary"))
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

correlation.median.summary <- data.frame(dataset = rep("Shalek", nrow(all.res)), 
           scenario = rep("missing5", nrow(all.res)),
           scimpute = all.res[, 1],
           lnpX = all.res[, 2], 
           saver = all.res[, 5],
           imputed.lasso = all.res[, 3],  
           magic = all.res[, 5],
           imputed.elasticnet  = all.res[, 6])           
  
   
save(correlation.median.summary, square.loss.median.summary, absolute.loss.median.summary, file = "Shalek_imputation_with_missing5_summary")



all.res <- NULL

for(i in 1:10){
	 load(paste0("Shalek_missing10_copy", i, "_imputation.summary"))
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

absolute.loss.median.summary <- data.frame(dataset = rep("Shalek", nrow(all.res)), 
           scenario = rep("missing10", nrow(all.res)),
           scimpute = all.res[, 1],
           lnpX = all.res[, 2], 
           saver = all.res[, 5],
           imputed.lasso = all.res[, 3],  
           magic = all.res[, 5],
           imputed.elasticnet  = all.res[, 6])           
  
  

all.res <- NULL

for(i in 1:10){
	 load(paste0("Shalek_missing10_copy", i, "_imputation.summary"))
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

square.loss.median.summary <- data.frame(dataset = rep("Shalek", nrow(all.res)), 
           scenario = rep("missing10", nrow(all.res)),
           scimpute = all.res[, 1],
           lnpX = all.res[, 2], 
           saver = all.res[, 5],
           imputed.lasso = all.res[, 3],  
           magic = all.res[, 5],
           imputed.elasticnet  = all.res[, 6])           
  
  
  

all.res <- NULL

for(i in 1:10){
	 load(paste0("Shalek_missing10_copy", i, "_imputation.summary"))
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

correlation.median.summary <- data.frame(dataset = rep("Shalek", nrow(all.res)), 
           scenario = rep("missing10", nrow(all.res)),
           scimpute = all.res[, 1],
           lnpX = all.res[, 2], 
           saver = all.res[, 5],
           imputed.lasso = all.res[, 3],  
           magic = all.res[, 5],
           imputed.elasticnet  = all.res[, 6])           
  
   
save(correlation.median.summary, square.loss.median.summary, absolute.loss.median.summary, file = "Shalek_imputation_with_missing10_summary")




