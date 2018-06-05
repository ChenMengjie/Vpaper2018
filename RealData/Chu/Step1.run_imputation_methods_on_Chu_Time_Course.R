load("GSE75748_sc_time_course.rda")
write.csv(gene.expression, "GSE75748_sc_time_course.csv")
write.csv(t(gene.expression), "GSE75748_sc_time_course_t.csv")


###############################################  
############# running drImpute ################ 
############################################### 

library(DrImpute)
load("GSE75748_sc_time_course.rda")
exdata <- preprocessSC(gene.expression)
sf <- apply(exdata, 2, mean)
npX <- t(t(exdata) / sf ) 
lnpX <- log(npX+1)
head(lnpX)
system.time(lnpX_imp <- DrImpute(lnpX))
save(lnpX_imp, file = "Chu_impute_tc_drImpute")


###############################################  
############# running SAVER ################### 
############################################### 
library(SAVER)
load("GSE75748_sc_time_course.rda")
system.time(saver_imp <- saver(gene.expression))
save(saver_imp, file = "Chu_impute_tc_SAVER")


###############################################  
############# running MAGIC ################### 
############################################### 



python MAGIC.py -d GSE75748_sc_time_course_t.csv -o GSE75748_sc_time_course_magic.csv --mols-per-cell-min 0 --mols-per-cell-max 100000000 csv



###############################################  
############# running scImpute ################ 
############################################### 

library(scImpute)
scimpute(count_path = "GSE75748_sc_time_course.csv",
  infile = "csv",           # format of input file
  outfile = "csv",          # format of output file
  out_dir = "Chu/scImpute1_trueK/",           # full path to output directory
  drop_thre = 0.5,          # threshold set on dropout probability
  Kcluster = 6,
  ncores = 1)



############################################ 
############# running VIPER ################ 
############################################

######## Lasso ########

load("GSE75748_sc_time_course.rda")

library(VIPER)

system.time(imp_ee <- VIPER(gene.expression, num = 5000, percentage.cutoff = 0.1, minbool = FALSE, alpha = 1))
save(imp_ee, file = "Chu_impute_tc_lasso")


######## Elastic net ######## 

load("GSE75748_sc_time_course.rda")
library(VIPER)
system.time(imp_ee <- VIPER(gene.expression, num = 5000, percentage.cutoff = 0.1, minbool = FALSE, alpha = 0.5))
save(imp_ee, file = "Chu_impute_tc_el")




