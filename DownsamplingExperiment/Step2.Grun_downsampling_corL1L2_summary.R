



for(i in c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95)){

	scimpute <- read.csv(paste0("Grun/scImpute/filtered_rate", i, "scimpute_count.csv"), as.is = T)
	scimpute.imp <- scimpute[, -1]

	load(paste0("Grun_drImpute_filtered_rate", i))	
	lnpX.imp <- lnpX_imp
	
	load(paste0("Grun_SAVER_filtered_rate", i))
	saver.imp <- saver_imp$estimate
	
	load(paste0("Grun_lasso_filtered_rate", i))
	imputed.imp <- imp_ee$imputed 
	
	magic <- read.csv(paste0("Grun_filtered_rate", i, "_magic.csv"))
	magic.imp <- t(magic[, -1])
	
	load(paste0("Grun_filtered_rate", i))
	imputed.el.imp <- imp_el_ee$imputed	
	 
	save(scimpute.imp, lnpX.imp, saver.imp, imputed.imp, magic.imp, imputed.el.imp,  file = paste0("Grun_filtered_rate", i, "_imputation.summary"))
	
}




for(i in c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95)){

	scimpute <- read.csv(paste0("Grun/scImpute/filtered_uniform_rate", i, "scimpute_count.csv"), as.is = T)
	scimpute.imp <- scimpute[, -1]

	load(paste0("Grun_drImpute_filtered_uniform_rate", i))	
	lnpX.imp <- lnpX_imp
	
	load(paste0("Grun_SAVER_filtered_uniform_rate", i))
	saver.imp <- saver_imp$estimate
	
	load(paste0("Grun_lasso_filtered_uniform_rate", i))
	imputed.imp <- imp_ee$imputed 
	
	magic <- read.csv(paste0("Grun_filtered_uniform_rate", i, "_magic.csv"))
	magic.imp <- t(magic[, -1])
	
	load(paste0("Grun_el_filtered_uniform_el_rate", i))
	imputed.el.imp <- imp_el_ee$imputed
	 
	save(scimpute.imp, lnpX.imp, saver.imp, imputed.imp, magic.imp, imputed.el.imp,  file = paste0("Grun_filtered_uniform_rate", i, "_imputation.summary"))
	
}


load("grun.rda")
xx <- gene.expression
aa <- gsub("\\_[0-9]*", "", colnames(xx))
types <- unique(aa)
gene.sum <- apply(xx, 1, sum)
xx <- xx[gene.sum > 0, ]

all.res <- NULL
all.res2 <- NULL
all.res3 <- NULL
for(i in c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95)){
	load(paste0("Grun_filtered_rate", i))

	sub.xx <- xx[rownames(gene.expression), colnames(gene.expression)]
	
	load(paste0("Grun_downsampling_flag_rate", i))
	
    sub.all.samples <- all.samples[rownames(gene.expression), colnames(gene.expression)]
	all.samples.before.dr <- all.samples.before.dr[rownames(gene.expression), colnames(gene.expression)]
	sub.all.ori <- all.ori[rownames(gene.expression), colnames(gene.expression)]
	
	load(paste0("Grun_filtered_rate", i, "_imputation.summary"))
	log.scimpute.pred <- apply(scimpute.imp, 2, function(x){log(x+1)})
	log.lnpX.pred <- lnpX.imp
	log.imputed.pred <- apply(imputed.imp, 2, function(x){log(exp(x)-0.1+1)})
	log.magic.pred <- apply(magic.imp, 2, function(x){log(x+1)})
	log.saver.pred <- apply(saver.imp, 2, function(x){log(x+1)})
	log.sub.xx <- apply(sub.xx, 2, function(x){log(x+1)})
	log.imputed.el.pred <- apply(imputed.el.imp, 2, function(x){log(exp(x)-0.1+1)})

	nonzero.flag <- sub.all.samples == 0 & sub.all.ori !=0  & all.samples.before.dr !=0

		
	all.res <- rbind(all.res, c(cor(log.sub.xx[nonzero.flag], log.scimpute.pred[nonzero.flag]),
	cor(log.sub.xx[nonzero.flag], log.lnpX.pred[nonzero.flag]),
	cor(log.sub.xx[nonzero.flag], log.saver.pred[nonzero.flag]),
	cor(log.sub.xx[nonzero.flag], log.imputed.pred[nonzero.flag]),
	cor(log.sub.xx[nonzero.flag], log.magic.pred[nonzero.flag]),
	cor(log.sub.xx[nonzero.flag], log.imputed.el.pred[nonzero.flag])))

	all.res2 <- rbind(all.res2, c(median(abs(log.sub.xx[nonzero.flag] - log.scimpute.pred[nonzero.flag])),
	median(abs(log.sub.xx[nonzero.flag] - log.lnpX.pred[nonzero.flag])),
	median(abs(log.sub.xx[nonzero.flag] - log.saver.pred[nonzero.flag])),
	median(abs(log.sub.xx[nonzero.flag] - log.imputed.pred[nonzero.flag])),
	median(abs(log.sub.xx[nonzero.flag] - log.magic.pred[nonzero.flag])),
	median(abs(log.sub.xx[nonzero.flag] - log.imputed.el.pred[nonzero.flag]))))
	
		
	all.res3 <- rbind(all.res3, c(median(abs(log.sub.xx[nonzero.flag] - log.scimpute.pred[nonzero.flag])^2),
	median(abs(log.sub.xx[nonzero.flag] - log.lnpX.pred[nonzero.flag])^2),
	median(abs(log.sub.xx[nonzero.flag] - log.saver.pred[nonzero.flag])^2),
	median(abs(log.sub.xx[nonzero.flag] - log.imputed.pred[nonzero.flag])^2),
	median(abs(log.sub.xx[nonzero.flag] - log.magic.pred[nonzero.flag])^2),
	median(abs(log.sub.xx[nonzero.flag] - log.imputed.el.pred[nonzero.flag])^2)))

}

correlation.summary1 <- data.frame(dataset = rep("Grun", nrow(all.res)), 
           downsampling = c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95),
           scimpute = all.res[, 1],
           lnpX = all.res[, 2], 
           saver = all.res[, 3],
           imputed = all.res[, 4],  
           magic = all.res[, 5],
           imputed.el = all.res[, 6])

L1.summary1 <- data.frame(dataset = rep("Grun", nrow(all.res)), 
           downsampling = c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95),
           scimpute = all.res2[, 1],
           lnpX = all.res2[, 2], 
           saver = all.res2[, 3],
           imputed = all.res2[, 4],  
           magic = all.res2[, 5],
           imputed.el = all.res2[, 6])

L2.summary1 <- data.frame(dataset = rep("Grun", nrow(all.res)), 
           downsampling = c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95),
           scimpute = all.res3[, 1],
           lnpX = all.res3[, 2], 
           saver = all.res3[, 3],
           imputed = all.res3[, 4],  
           magic = all.res3[, 5],
           imputed.el = all.res3[, 6])





load("grun.rda")
xx <- gene.expression
aa <- gsub("\\_[0-9]*", "", colnames(xx))
types <- unique(aa)
gene.sum <- apply(xx, 1, sum)
xx <- xx[gene.sum > 0, ]

all.res <- NULL
all.res2 <- NULL
all.res3 <- NULL
for(i in c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95)){
	load(paste0("Grun_filtered_rate", i))

	sub.xx <- xx[rownames(gene.expression), colnames(gene.expression)]
	
	load(paste0("Grun_downsampling_flag_rate", i))
	
    sub.all.samples <- all.samples[rownames(gene.expression), colnames(gene.expression)]
	all.samples.before.dr <- all.samples.before.dr[rownames(gene.expression), colnames(gene.expression)]
	sub.all.ori <- all.ori[rownames(gene.expression), colnames(gene.expression)]
	
	load(paste0("Grun_filtered_rate", i, "_imputation.summary"))

	log.scimpute.pred <- apply(scimpute.imp, 2, function(x){log(x+1)})
	log.lnpX.pred <- lnpX.imp
	log.imputed.pred <- apply(imputed.imp, 2, function(x){log(exp(x)-0.1+1)})
	log.magic.pred <- apply(magic.imp, 2, function(x){log(x+1)})
	log.saver.pred <- apply(saver.imp, 2, function(x){log(x+1)})
	log.sub.xx <- apply(sub.xx, 2, function(x){log(x+1)})
	log.imputed.el.pred <- apply(imputed.el.imp, 2, function(x){log(exp(x)-0.1+1)})

	zero.flag <- sub.all.samples == 0 & sub.all.ori !=0  & all.samples.before.dr == 0
 	
	
	all.res <- rbind(all.res, c(cor(log.sub.xx[zero.flag], log.scimpute.pred[zero.flag]),
	cor(log.sub.xx[zero.flag], log.lnpX.pred[zero.flag]),
	cor(log.sub.xx[zero.flag], log.saver.pred[zero.flag]),
	cor(log.sub.xx[zero.flag], log.imputed.pred[zero.flag]),
	cor(log.sub.xx[zero.flag], log.magic.pred[zero.flag]),
	cor(log.sub.xx[zero.flag], log.imputed.el.pred[zero.flag])))

	
	all.res2 <- rbind(all.res2, c(median(abs(log.sub.xx[zero.flag] - log.scimpute.pred[zero.flag])),
	median(abs(log.sub.xx[zero.flag] - log.lnpX.pred[zero.flag])),
	median(abs(log.sub.xx[zero.flag] - log.saver.pred[zero.flag])),
	median(abs(log.sub.xx[zero.flag] - log.imputed.pred[zero.flag])),
	median(abs(log.sub.xx[zero.flag] - log.magic.pred[zero.flag])),
	median(abs(log.sub.xx[zero.flag] - log.imputed.el.pred[zero.flag]))))
	
		
	all.res3 <- rbind(all.res3, c(median(abs(log.sub.xx[zero.flag] - log.scimpute.pred[zero.flag])^2),
	median(abs(log.sub.xx[zero.flag] - log.lnpX.pred[zero.flag])^2),
	median(abs(log.sub.xx[zero.flag] - log.saver.pred[zero.flag])^2),
	median(abs(log.sub.xx[zero.flag] - log.imputed.pred[zero.flag])^2),
	median(abs(log.sub.xx[zero.flag] - log.magic.pred[zero.flag])^2),
	median(abs(log.sub.xx[zero.flag] - log.imputed.el.pred[zero.flag])^2)))

}

correlation.summary2 <- data.frame(dataset = rep("Grun", nrow(all.res)), 
           downsampling = c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95),
           scimpute = all.res[, 1],
           lnpX = all.res[, 2], 
           saver = all.res[, 3],
           imputed = all.res[, 4],  
           magic = all.res[, 5],
           imputed.el = all.res[, 6])

L1.summary2 <- data.frame(dataset = rep("Grun", nrow(all.res)), 
           downsampling = c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95),
           scimpute = all.res2[, 1],
           lnpX = all.res2[, 2], 
           saver = all.res2[, 3],
           imputed = all.res2[, 4],  
           magic = all.res2[, 5],
           imputed.el = all.res2[, 6])

L2.summary2 <- data.frame(dataset = rep("Grun", nrow(all.res)), 
           downsampling = c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95),
           scimpute = all.res3[, 1],
           lnpX = all.res3[, 2], 
           saver = all.res3[, 3],
           imputed = all.res3[, 4],  
           magic = all.res3[, 5],
           imputed.el = all.res3[, 6])
  
save(correlation.summary1, L1.summary1, L2.summary1, correlation.summary2, L1.summary2, L2.summary2, file = "Grun_downsampling_comparison_summary_CorL1L2")
           

          


           
 
load("grun.rda")
xx <- gene.expression
aa <- gsub("\\_[0-9]*", "", colnames(xx))
types <- unique(aa)
gene.sum <- apply(xx, 1, sum)
xx <- xx[gene.sum > 0, ]

all.res <- NULL
all.res2 <- NULL
all.res3 <- NULL
for(i in c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95)){
	load(paste0("Grun_filtered_uniform_rate", i))

	sub.xx <- xx[rownames(gene.expression), colnames(gene.expression)]
	
	load(paste0("Grun_uniform_downsampling_flag_rate", i))
	
    sub.all.samples <- all.samples[rownames(gene.expression), colnames(gene.expression)]
	all.samples.before.dr <- all.samples.before.dr[rownames(gene.expression), colnames(gene.expression)]
	sub.all.ori <- all.ori[rownames(gene.expression), colnames(gene.expression)]
	
	load(paste0("Grun_filtered_uniform_rate", i, "_imputation.summary"))
	log.scimpute.pred <- apply(scimpute.imp, 2, function(x){log(x+1)})
	log.lnpX.pred <- lnpX.imp
	log.imputed.pred <- apply(imputed.imp, 2, function(x){log(exp(x)-0.1+1)})
	log.magic.pred <- apply(magic.imp, 2, function(x){log(x+1)})
	log.saver.pred <- apply(saver.imp, 2, function(x){log(x+1)})
	log.sub.xx <- apply(sub.xx, 2, function(x){log(x+1)})
	log.imputed.el.pred <- apply(imputed.el.imp, 2, function(x){log(exp(x)-0.1+1)})

	nonzero.flag <- sub.all.samples == 0 & sub.all.ori !=0  & all.samples.before.dr !=0
	
		
	all.res <- rbind(all.res, c(cor(log.sub.xx[nonzero.flag], log.scimpute.pred[nonzero.flag]),
	cor(log.sub.xx[nonzero.flag], log.lnpX.pred[nonzero.flag]),
	cor(log.sub.xx[nonzero.flag], log.saver.pred[nonzero.flag]),
	cor(log.sub.xx[nonzero.flag], log.imputed.pred[nonzero.flag]),
	cor(log.sub.xx[nonzero.flag], log.magic.pred[nonzero.flag]),
	cor(log.sub.xx[nonzero.flag], log.imputed.el.pred[nonzero.flag])))

	all.res2 <- rbind(all.res2, c(median(abs(log.sub.xx[nonzero.flag] - log.scimpute.pred[nonzero.flag])),
	median(abs(log.sub.xx[nonzero.flag] - log.lnpX.pred[nonzero.flag])),
	median(abs(log.sub.xx[nonzero.flag] - log.saver.pred[nonzero.flag])),
	median(abs(log.sub.xx[nonzero.flag] - log.imputed.pred[nonzero.flag])),
	median(abs(log.sub.xx[nonzero.flag] - log.magic.pred[nonzero.flag])),
	median(abs(log.sub.xx[nonzero.flag] - log.imputed.el.pred[nonzero.flag]))))
	
		
	all.res3 <- rbind(all.res3, c(median(abs(log.sub.xx[nonzero.flag] - log.scimpute.pred[nonzero.flag])^2),
	median(abs(log.sub.xx[nonzero.flag] - log.lnpX.pred[nonzero.flag])^2),
	median(abs(log.sub.xx[nonzero.flag] - log.saver.pred[nonzero.flag])^2),
	median(abs(log.sub.xx[nonzero.flag] - log.imputed.pred[nonzero.flag])^2),
	median(abs(log.sub.xx[nonzero.flag] - log.magic.pred[nonzero.flag])^2),
	median(abs(log.sub.xx[nonzero.flag] - log.imputed.el.pred[nonzero.flag])^2)))

}

correlation.summary1 <- data.frame(dataset = rep("Grun", nrow(all.res)), 
           downsampling = c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95),
           scimpute = all.res[, 1],
           lnpX = all.res[, 2], 
           saver = all.res[, 3],
           imputed = all.res[, 4],  
           magic = all.res[, 5],
           imputed.el = all.res[, 6])

L1.summary1 <- data.frame(dataset = rep("Grun", nrow(all.res)), 
           downsampling = c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95),
           scimpute = all.res2[, 1],
           lnpX = all.res2[, 2], 
           saver = all.res2[, 3],
           imputed = all.res2[, 4],  
           magic = all.res2[, 5],
           imputed.el = all.res2[, 6])

L2.summary1 <- data.frame(dataset = rep("Grun", nrow(all.res)), 
           downsampling = c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95),
           scimpute = all.res3[, 1],
           lnpX = all.res3[, 2], 
           saver = all.res3[, 3],
           imputed = all.res3[, 4],  
           magic = all.res3[, 5],
           imputed.el = all.res3[, 6])



load("grun.rda")
xx <- gene.expression
aa <- gsub("\\_[0-9]*", "", colnames(xx))
types <- unique(aa)
gene.sum <- apply(xx, 1, sum)
xx <- xx[gene.sum > 0, ]

all.res <- NULL
all.res2 <- NULL
all.res3 <- NULL
for(i in c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95)){
	load(paste0("Grun_filtered_uniform_rate", i))

	sub.xx <- xx[rownames(gene.expression), colnames(gene.expression)]
	
	load(paste0("Grun_uniform_downsampling_flag_rate", i))
	
    sub.all.samples <- all.samples[rownames(gene.expression), colnames(gene.expression)]
	all.samples.before.dr <- all.samples.before.dr[rownames(gene.expression), colnames(gene.expression)]
	sub.all.ori <- all.ori[rownames(gene.expression), colnames(gene.expression)]
	
	load(paste0("Grun_filtered_uniform_rate", i, "_imputation.summary"))
	log.scimpute.pred <- apply(scimpute.imp, 2, function(x){log(x+1)})
	log.lnpX.pred <- lnpX.imp
	log.imputed.pred <- apply(imputed.imp, 2, function(x){log(exp(x)-0.1+1)})
	log.magic.pred <- apply(magic.imp, 2, function(x){log(x+1)})
	log.saver.pred <- apply(saver.imp, 2, function(x){log(x+1)})
	log.sub.xx <- apply(sub.xx, 2, function(x){log(x+1)})
	log.imputed.el.pred <- apply(imputed.el.imp, 2, function(x){log(exp(x)-0.1+1)})

	zero.flag <- sub.all.samples == 0 & sub.all.ori !=0  & all.samples.before.dr == 0
	
	all.res <- rbind(all.res, c(cor(log.sub.xx[zero.flag], log.scimpute.pred[zero.flag]),
	cor(log.sub.xx[zero.flag], log.lnpX.pred[zero.flag]),
	cor(log.sub.xx[zero.flag], log.saver.pred[zero.flag]),
	cor(log.sub.xx[zero.flag], log.imputed.pred[zero.flag]),
	cor(log.sub.xx[zero.flag], log.magic.pred[zero.flag]),
	cor(log.sub.xx[zero.flag], log.imputed.el.pred[zero.flag])))

	
	all.res2 <- rbind(all.res2, c(median(abs(log.sub.xx[zero.flag] - log.scimpute.pred[zero.flag])),
	median(abs(log.sub.xx[zero.flag] - log.lnpX.pred[zero.flag])),
	median(abs(log.sub.xx[zero.flag] - log.saver.pred[zero.flag])),
	median(abs(log.sub.xx[zero.flag] - log.imputed.pred[zero.flag])),
	median(abs(log.sub.xx[zero.flag] - log.magic.pred[zero.flag])),
	median(abs(log.sub.xx[zero.flag] - log.imputed.el.pred[zero.flag]))))
	
		
	all.res3 <- rbind(all.res3, c(median(abs(log.sub.xx[zero.flag] - log.scimpute.pred[zero.flag])^2),
	median(abs(log.sub.xx[zero.flag] - log.lnpX.pred[zero.flag])^2),
	median(abs(log.sub.xx[zero.flag] - log.saver.pred[zero.flag])^2),
	median(abs(log.sub.xx[zero.flag] - log.imputed.pred[zero.flag])^2),
	median(abs(log.sub.xx[zero.flag] - log.magic.pred[zero.flag])^2),
	median(abs(log.sub.xx[zero.flag] - log.imputed.el.pred[zero.flag])^2)))

}

correlation.summary2 <- data.frame(dataset = rep("Grun", nrow(all.res)), 
           downsampling = c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95),
           scimpute = all.res[, 1],
           lnpX = all.res[, 2], 
           saver = all.res[, 3],
           imputed = all.res[, 4],  
           magic = all.res[, 5],
           imputed.el = all.res[, 6])

L1.summary2 <- data.frame(dataset = rep("Grun", nrow(all.res)), 
           downsampling = c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95),
           scimpute = all.res2[, 1],
           lnpX = all.res2[, 2], 
           saver = all.res2[, 3],
           imputed = all.res2[, 4],  
           magic = all.res2[, 5],
           imputed.el = all.res2[, 6])

L2.summary2 <- data.frame(dataset = rep("Grun", nrow(all.res)), 
           downsampling = c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95),
           scimpute = all.res3[, 1],
           lnpX = all.res3[, 2], 
           saver = all.res3[, 3],
           imputed = all.res3[, 4],  
           magic = all.res3[, 5],
           imputed.el = all.res3[, 6])
           
save(correlation.summary1, L1.summary1, L2.summary1, correlation.summary2, L1.summary2, L2.summary2, file = "Grun_downsampling_uniform_comparison_summary_CorL1L2")
          
