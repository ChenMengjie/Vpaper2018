

################################################
########. Generate Downsampled data ############
################################################
le <- 0.95 # 0.9, 0.8, 0.7, 0.6, 0.5

log.pi.res <- function(z){
	exp(z)/(1+exp(z))
}

sample.for.dropout <- function(z){

	dropout.list <- NULL
	for(i in 1:length(z)){
		dropout.list <- c(dropout.list, sample(0:1, 1, prob=c(z[i], 1-z[i])))
	}
	return(dropout.list)
}


change.rate <- function(x1, x2, coef){

	log.x1 <- log2(x1 + 0.1)
	log.x2 <- log2(x2 + 0.1)
	log.x1.hat <- coef[1] + log.x1*coef[2]
	log.x2.hat <- coef[1] + log.x2*coef[2]
	rate.x1 <- log.pi.res(log.x1.hat)
	rate.x2 <- log.pi.res(log.x2.hat)
	dropout <- (rate.x2 - rate.x1)/(1 - rate.x1)
	dropout[dropout<0] <- 0.0001
	dropout.label <- sample.for.dropout(dropout)
	return(x2*dropout.label)
}


xx <- read.csv("GSE75748_sc_time_course_ec.csv", as.is = T)
rownames(xx) <- xx[, 1]
xx <- xx[, -1]
aa <- gsub("\\_[^<>]*", "", colnames(xx))
types <- unique(aa)
gene.sum <- apply(xx, 1, sum)
xx <- xx[gene.sum > 0, ]
p <- nrow(xx)
template <- data.frame(id = 1:p, vv = rep(1, p))


all.samples <- NULL
all.samples.before.dr <- NULL
all.samples.names <- NULL
all.ori <- NULL
for(i in 1:length(types)){

	given.type <- types[i]
	selected <- which(aa%in%given.type)
	selected.x <- xx[, selected]
	sizes <- apply(selected.x, 2, sum)
	flag <- log10(sizes) >= 6
	selected.x <- selected.x[, flag]
	logxx <- apply(selected.x, 2, function(x){log2(x+0.1)})
	n <- ncol(selected.x)
	logxx.mean <- apply(logxx, 1, mean)
	logxx.zero <- apply(selected.x, 1, function(y){length(y[y == 0])})
	zero.rate <- round(logxx.zero/n, 2)
	pi_ratio <- log(zero.rate/(1-zero.rate))
	pi_ratio[is.infinite(pi_ratio)] <- NA
	coeff <- summary(lm(pi_ratio~logxx.mean))$coefficients[, 1]
    all.ori <- rbind(all.ori, t(selected.x))
	all.samples.names <- c(all.samples.names, colnames(selected.x))
	for(j in 1:n){

		reads <- sum(selected.x[, j])
		prob <- selected.x[, j]/sum(selected.x[, j])

		new.reads <- round(reads*le)
		newsample <- sample(1:p, size = new.reads, prob = prob, replace = TRUE)
		newsample.tab <- table(newsample)
		newsample.mat <- data.frame(newsample.tab)
		merged <- merge(template, newsample.mat, by.x = "id", by.y = "newsample", all = TRUE)
		yy <- merged[, 3]
		all.samples.before.dr <- cbind(all.samples.before.dr, yy)
		flag2 <- !is.na(yy)
		subyy <- yy[flag2]
		original <- selected.x[flag2, j]
		drop.subyy <- change.rate(original, subyy, coeff)
		yy[flag2] <- drop.subyy

		all.samples <- cbind(all.samples, yy)

	}

}


rownames(all.samples) <- rownames(xx)
colnames(all.samples) <- all.samples.names
all.samples[is.na(all.samples)] <- 0
save(all.samples, file = paste0("Chu_time_course_downsampling_rate", le))



rownames(all.samples.before.dr) <- rownames(xx)
colnames(all.samples.before.dr) <- all.samples.names
all.samples.before.dr[is.na(all.samples.before.dr)] <- 0

all.ori <- t(all.ori)
save(all.samples, all.samples.before.dr, all.ori, file = paste0("Chu_time_course_flag_rate", le))

s1 <- table(all.samples==0 & all.samples.before.dr ==0 & all.ori!=0)
s2 <- table(all.samples==0 & all.samples.before.dr !=0 & all.ori!=0)

save(s1, s2, file = paste0("Chu_time_course_downsampling_zero_summary_rate", le))


############################################
############# prepare input ################ 
############################################ 

for(le in c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95)){

	load(paste0("Chu_time_course_downsampling_rate", le))
	xx <- all.samples
	n <- ncol(xx)
	zero.rate <- apply(xx, 1, function(x){length(x[x == 0])})
	flag2 <- zero.rate <= 0.9*n
	gene.expression <- xx[flag2, ]
	save(gene.expression, file = paste0("Chu_time_course_filtered_rate", le))
	write.csv(gene.expression, paste0("Chu_time_course_filtered_rate", le, ".csv"))
	write.csv(t(gene.expression), paste0("Chu_time_course_filtered_rate", le, "_t.csv"))

}

###############################################  
############# running scImpute ################ 
############################################### 
le <- 0.95 # 0.9, 0.8, 0.7, 0.6, 0.5

library(scImpute)
system(paste0("mkdir scImpute/time_course/filtered_rate", le))
system.time(scimpute(count_path = paste0("Chu_time_course_filtered_rate", le, ".csv"),
  infile = "csv",           # format of input file
  outfile = "csv",          # format of output file
  out_dir = paste0("scImpute/time_course/filtered_rate", le),  # full path to output directory
  drop_thre = 0.5,          # threshold set on dropout probability
  Kcluster = 6,
  ncores = 1))


###############################################  
############# running drImpute ################ 
############################################### 

load(paste0("Chu_time_course_filtered_rate", le))
library(DrImpute)
exdata <- preprocessSC(gene.expression)
sf <- apply(exdata, 2, mean)
npX <- t(t(exdata)/sf)
lnpX <- log(npX+1)
head(lnpX)
system.time(lnpX_imp <- DrImpute(lnpX))
save(lnpX_imp, file = paste0("Chu_drImpute_time_course_filtered_rate", le))


###############################################  
############# running SAVER ################### 
############################################### 

library(SAVER)
system.time(saver_imp <- saver(gene.expression))
save(saver_imp, file = paste0("Chu_SAVER_time_course_filtered_rate", le))


###############################################  
############# running MAGIC ################### 
############################################### 


python MAGIC.py -d Chu_time_course_filtered_rate0.5_t.csv -o Chu_time_course_filtered_rate0.5_magic.csv --mols-per-cell-min 0 --mols-per-cell-max 100000000 csv
python MAGIC.py -d Chu_time_course_filtered_rate0.6_t.csv -o Chu_time_course_filtered_rate0.6_magic.csv --mols-per-cell-min 0 --mols-per-cell-max 100000000 csv
python MAGIC.py -d Chu_time_course_filtered_rate0.7_t.csv -o Chu_time_course_filtered_rate0.7_magic.csv --mols-per-cell-min 0 --mols-per-cell-max 100000000 csv
python MAGIC.py -d Chu_time_course_filtered_rate0.8_t.csv -o Chu_time_course_filtered_rate0.8_magic.csv --mols-per-cell-min 0 --mols-per-cell-max 100000000 csv
python MAGIC.py -d Chu_time_course_filtered_rate0.9_t.csv -o Chu_time_course_filtered_rate0.9_magic.csv --mols-per-cell-min 0 --mols-per-cell-max 100000000 csv
python MAGIC.py -d Chu_time_course_filtered_rate0.95_t.csv -o Chu_time_course_filtered_rate0.95_magic.csv --mols-per-cell-min 0 --mols-per-cell-max 100000000 csv










############################################ 
############# running VIPER ################ 
############################################

######## Lasso ########

le <- 0.95 # 0.9, 0.8, 0.7, 0.6, 0.5

load(paste0("Chu_time_course_filtered_rate", le))
library(VIPER)

system.time(imp_ee <- VIPER(gene.expression, num = 5000, percentage.cutoff = 0.1, minbool = FALSE, alpha = 1))
save(imp_ee, file = paste0("Chu_lasso_time_course_filtered_rate", le))

######## Elastic net ######## 

le <- 0.95 # 0.9, 0.8, 0.7, 0.6, 0.5

load(paste0("Chu_time_course_filtered_rate", le))
library(VIPER)
system.time(imp_el_ee <- VIPER(gene.expression, num = 5000, percentage.cutoff = 0.1, minbool = FALSE, alpha = 0.5))
save(imp_el_ee, file = paste0("Chu_el_time_course_filtered_rate", le))




