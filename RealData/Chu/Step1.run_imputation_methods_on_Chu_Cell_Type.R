load("GSE75748_sc_cell_type.rda")
write.csv(gene.expression, "GSE75748_sc_cell_type.csv")
write.csv(t(gene.expression), "GSE75748_sc_cell_type_t.csv")


###############################################  
############# running drImpute ################ 
############################################### 
library(DrImpute)
load("GSE75748_sc_cell_type.rda")
exdata <- preprocessSC(gene.expression)
sf <- apply(exdata, 2, mean)
npX <- t(t(exdata) / sf ) 
lnpX <- log(npX+1)
head(lnpX)
system.time(lnpX_imp <- DrImpute(lnpX))
save(lnpX_imp, file = "Chu_impute_ct_drImpute")



###############################################  
############# running SAVER ################### 
############################################### 

library(SAVER)
load("GSE75748_sc_cell_type.rda")
system.time(saver_imp <- saver(gene.expression))
save(saver_imp, file = "Chu_impute_ct_SAVER")



###############################################  
############# running MAGIC ################### 
############################################### 


python MAGIC.py -d GSE75748_sc_cell_type_t.csv -o GSE75748_sc_cell_type_magic.csv --mols-per-cell-min 0 --mols-per-cell-max 100000000 csv

###############################################  
############# running scImpute ################ 
############################################### 

library(scImpute)
scimpute(count_path = "GSE75748_sc_cell_type.csv",
         infile = "csv",           # format of input file
         outfile = "csv",          # format of output file
         out_dir = "Chu/scImpute2_trueK/",           # full path to output directory
         drop_thre = 0.5,          # threshold set on dropout probability
         Kcluster = 7,
         ncores = 1)


############################################ 
############# running VIPER ################ 
############################################

######## Lasso ########

load("GSE75748_sc_cell_type.rda")

library(VIPER)

system.time(imp_ee <- VIPER(gene.expression, num = 5000, percentage.cutoff = 0.1, minbool = FALSE, alpha = 1))
save(imp_ee, file = "Chu_impute_ct_lasso")


######## Elastic net ######## 

load("GSE75748_sc_cell_type.rda")
library(VIPER)
system.time(imp_ee <- VIPER(gene.expression, num = 5000, percentage.cutoff = 0.1, minbool = FALSE, alpha = 0.5))
save(imp_ee, file = "Chu_impute_ct_el")






