###########################################
########## creating missing data ##########
###########################################
load("Shalek_filtered.rda")
for(i in 1:5){
	set.seed(i*500+80224523)
	gene.expression.missing <- gene.expression
	len <- length(gene.expression[gene.expression!=0])
	missing <- sample(c(0, 1), len, replace=T, prob=c(0.02, 0.98))
	gene.expression.missing[gene.expression!=0] <- gene.expression[gene.expression!=0]*missing
    save(gene.expression.missing, file = paste0("Shalek_filtered_missing2_copy", i))
    write.csv(gene.expression.missing, paste0("Shalek_filtered_missing2_copy", i, ".csv"))
	write.csv(t(gene.expression.missing), paste0("Shalek_filtered_missing2_copy_t", i, ".csv"))
}


###############################################  
############# running scImpute ################ 
############################################### 

i <- 1 # 2,3,4,5
library(scImpute)
system.time(scimpute(count_path = paste0("Shalek_filtered_missing2_copy", i, ".csv"),
  infile = "csv",           # format of input file
  outfile = "csv",          # format of output file
  out_dir = paste0("scImpute/Shalek/missing2_copy", i),  # full path to output directory
  drop_thre = 0.5,          # threshold set on dropout probability
  Kcluster = 10,
  ncores = 1))
  


###############################################  
############# running drImpute ################ 
############################################### 

i <- 1 # 2,3,4,5
load(paste0("Shalek_filtered_missing2_copy", i))
library(DrImpute)
exdata <- preprocessSC(gene.expression.missing)
sf <- apply(exdata, 2, mean)
npX <- t(t(exdata)/sf) 
lnpX <- log(npX+1)
head(lnpX)
system.time(lnpX_imp <- DrImpute(lnpX))
save(lnpX_imp, file = paste0("Shalek_filtered_drImpute_missing2_copy", i))

###############################################  
############# running SAVER ################### 
############################################### 

load(paste0("Shalek_filtered_missing2_copy", i))
library(SAVER)
system.time(saver_imp <- saver(gene.expression.missing))
save(saver_imp, file = paste0("Shalek_filtered_SAVER_missing2_copy", i))


###############################################  
############# running MAGIC ################### 
############################################### 
python MAGIC.py -d Shalek_filtered_missing2_copy_t1.csv -o Shalek_filtered_missing2_copy1_magic.csv --mols-per-cell-min 0 --mols-per-cell-max 100000000 csv
python MAGIC.py -d Shalek_filtered_missing2_copy_t2.csv -o Shalek_filtered_missing2_copy2_magic.csv --mols-per-cell-min 0 --mols-per-cell-max 100000000 csv
python MAGIC.py -d Shalek_filtered_missing2_copy_t3.csv -o Shalek_filtered_missing2_copy3_magic.csv --mols-per-cell-min 0 --mols-per-cell-max 100000000 csv
python MAGIC.py -d Shalek_filtered_missing2_copy_t4.csv -o Shalek_filtered_missing2_copy4_magic.csv --mols-per-cell-min 0 --mols-per-cell-max 100000000 csv
python MAGIC.py -d Shalek_filtered_missing2_copy_t5.csv -o Shalek_filtered_missing2_copy5_magic.csv --mols-per-cell-min 0 --mols-per-cell-max 100000000 csv





############################################ 
############# running VIPER ################ 
############################################

######## Lasso ########

load(paste0("Shalek_filtered_missing2_copy", i))
library(VIPER)

system.time(imp_ee <- VIPER(gene.expression, num = 5000, percentage.cutoff = 0.1, minbool = FALSE, alpha = 1))
save(imp_ee, file = paste0("Shalek_filtered_lasso_missing2_copy", i))

######## Elastic net ######## 


load(paste0("Shalek_filtered_missing2_copy", i))
library(VIPER)
system.time(imp_el_ee <- VIPER(gene.expression, num = 5000, percentage.cutoff = 0.1, minbool = FALSE, alpha = 0.5))
save(imp_el_ee, file = paste0("Shalek_filtered_el_missing2_copy", i))



