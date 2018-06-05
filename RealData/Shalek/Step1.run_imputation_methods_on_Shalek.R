
load("Shalek_filtered.rda")
write.csv(gene.expression, "Shalek_filtered.csv")
write.csv(t(gene.expression), "Shalek_filtered_t.csv")

  
###############################################  
############# running drImpute ################ 
############################################### 

load("Shalek_filtered.rda")
library(DrImpute)
exdata <- preprocessSC(gene.expression)
sf <- apply(exdata, 2, mean)
npX <- t(t(exdata)/sf) 
lnpX <- log(npX+1)
head(lnpX)
system.time(lnpX_imp <- DrImpute(lnpX))
save(lnpX_imp, file = "Shalek_drImpute_res")


###############################################  
############# running SAVER ################### 
############################################### 
load("Shalek_filtered.rda")
library(SAVER)
system.time(saver_imp <- saver(gene.expression))
save(saver_imp, file = "Shalek_SAVER_res")

  
###############################################  
############# running MAGIC ################### 
############################################### 

python MAGIC.py -d Shalek_filtered_t.csv -o Shalek_filtered_magic.csv --mols-per-cell-min 0 --mols-per-cell-max 100000000 csv


###############################################  
############# running scImpute ################ 
############################################### 

library(scImpute)
system.time(scimpute(count_path = "Shalek_filtered.csv",
  infile = "csv",           # format of input file
  outfile = "csv",          # format of output file
  out_dir = "scImpute/Shalek_",  # full path to output directory
  drop_thre = 0.5,          # threshold set on dropout probability
  Kcluster = 10,
  ncores = 1))

 
 

############################################ 
############# running VIPER ################ 
############################################

######## Lasso ########

load("Shalek_filtered.rda")

library(VIPER)

system.time(imp_ee <- VIPER(gene.expression, num = 5000, percentage.cutoff = 0.1, minbool = FALSE, alpha = 1))
save(imp_ee, file = "Shalek_lasso_res")


######## Elastic net ######## 

load("Shalek_filtered.rda")
library(VIPER)
system.time(imp_ee <- VIPER(gene.expression, num = 5000, percentage.cutoff = 0.1, minbool = FALSE, alpha = 0.5))
save(imp_ee, file = "Shalek_el_res")






 