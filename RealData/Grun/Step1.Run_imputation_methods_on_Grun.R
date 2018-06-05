
load("grun.rda")
write.csv(gene.expression, "grun.csv")
write.csv(t(gene.expression), "grun_t.csv")



###############################################  
############# running drImpute ################ 
############################################### 

load("grun.rda")
library(DrImpute)
exdata <- preprocessSC(gene.expression)
sf <- apply(exdata, 2, mean)
npX <- t(t(exdata)/sf) 
lnpX <- log(npX+1)
head(lnpX)
system.time(lnpX_imp <- DrImpute(lnpX))
save(lnpX_imp, file = "grun_drImpute_res")


###############################################  
############# running SAVER ################### 
############################################### 
load("grun.rda")
library(SAVER)
system.time(saver_imp <- saver(gene.expression))
save(saver_imp, file = "grun_SAVER_res")

 
  
###############################################  
############# running MAGIC ################### 
############################################### 


python MAGIC.py -d grun_t.csv -o grun_magic.csv --mols-per-cell-min 0 --mols-per-cell-max 100000000 csv


###############################################  
############# running scImpute ################ 
############################################### 
library(scImpute)
system.time(scimpute(count_path = "grun.csv",
  infile = "csv",           # format of input file
  outfile = "csv",          # format of output file
  out_dir = "scImpute/grun_",  # full path to output directory
  drop_thre = 0.5,          # threshold set on dropout probability
  Kcluster = 4,
  ncores = 1))

 
 
 
############################################ 
############# running VIPER ################ 
############################################

######## Lasso ########

load("grun.rda")

library(VIPER)

system.time(imp_ee <- VIPER(gene.expression, num = 5000, percentage.cutoff = 0.1, minbool = FALSE, alpha = 1))
save(imp_ee, file = "grun_lasso_res")


load("grun.rda")

load("Shalek_filtered.rda")
library(VIPER)
system.time(imp_ee <- VIPER(gene.expression, num = 5000, percentage.cutoff = 0.1, minbool = FALSE, alpha = 0.5))
save(imp_ee, file = "grun_el_res")






 