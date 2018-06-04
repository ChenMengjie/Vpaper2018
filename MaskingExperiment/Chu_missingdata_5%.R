###########################################
########## creating missing data ##########
###########################################
load("GSE75748_sc_time_course.rda")
for(i in 1:5){
	set.seed(i*500+80224523)
	gene.expression.missing <- gene.expression
	len <- length(gene.expression[gene.expression!=0])
	missing <- sample(c(0, 1), len, replace=T, prob=c(0.05, 0.95))
	gene.expression.missing[gene.expression!=0] <- gene.expression[gene.expression!=0]*missing
    save(gene.expression.missing, file = paste0("Chu_time_course_missing5_copy", i))
    write.csv(gene.expression.missing, paste0("Chu_time_course_missing5_copy", i, ".csv"))
	write.csv(t(gene.expression.missing), paste0("Chu_time_course_missing5_copy_t", i, ".csv"))
}


load("GSE75748_sc_cell_type.rda")
for(i in 1:5){
	set.seed(i*500+80224523)
	gene.expression.missing <- gene.expression
	len <- length(gene.expression[gene.expression!=0])
	missing <- sample(c(0, 1), len, replace=T, prob=c(0.05, 0.95))
	gene.expression.missing[gene.expression!=0] <- gene.expression[gene.expression!=0]*missing
    save(gene.expression.missing, file = paste0("Chu_cell_type_missing5_copy", i))
    write.csv(gene.expression.missing, paste0("Chu_cell_type_missing5_copy", i, ".csv"))
	write.csv(t(gene.expression.missing), paste0("Chu_cell_type_missing5_copy_t", i, ".csv"))
}

###############################################  
############# running scImpute ################ 
############################################### 

i <- 1
library(scImpute)
system.time(scimpute(count_path = paste0("Chu_time_course_missing5_copy", i, ".csv"),
  infile = "csv",           # format of input file
  outfile = "csv",          # format of output file
  out_dir = paste0("scImpute_missing5_copy", i),  # full path to output directory
  drop_thre = 0.5,          # threshold set on dropout probability
  Kcluster = 6,
  ncores = 1))
  


i <- 1
library(scImpute)
system.time(scimpute(count_path = paste0("Chu_cell_type_missing5_copy", i, ".csv"),
  infile = "csv",           # format of input file
  outfile = "csv",          # format of output file
  out_dir = paste0("scImpute_missing5_copy", i),  # full path to output directory
  drop_thre = 0.5,          # threshold set on dropout probability
  Kcluster = 7,
  ncores = 1))
   


###############################################  
############# running drImpute ################ 
############################################### 

load(paste0("Chu_time_course_missing5_copy", i))
library(DrImpute)
exdata <- preprocessSC(gene.expression.missing)
sf <- apply(exdata, 2, mean)
npX <- t(t(exdata)/sf) 
lnpX <- log(npX+1)
head(lnpX)
system.time(lnpX_imp <- DrImpute(lnpX))
save(lnpX_imp, file = paste0("Chu_drImpute_time_course_missing5_copy", i))



load(paste0("Chu_cell_type_missing5_copy", i))

library(DrImpute)
exdata <- preprocessSC(gene.expression.missing)
sf <- apply(exdata, 2, mean)
npX <- t(t(exdata)/sf) 
lnpX <- log(npX+1)
head(lnpX)
system.time(lnpX_imp <- DrImpute(lnpX))
save(lnpX_imp, file = paste0("Chu_drImpute_cell_type_missing5_copy", i))



###############################################  
############# running SAVER ################### 
############################################### 

load(paste0("Chu_time_course_missing5_copy", i))
library(SAVER)
system.time(saver_imp <- saver(gene.expression.missing))
save(saver_imp, file = paste0("Chu_SAVER_time_course_missing5_copy", i))


load(paste0("Chu_cell_type_missing5_copy", i))
library(SAVER)
system.time(saver_imp <- saver(gene.expression.missing))
save(saver_imp, file = paste0("Chu_SAVER_cell_type_missing5_copy", i))




###############################################  
############# running MAGIC ################### 
############################################### 
python MAGIC.py -d Chu_cell_type_missing5_copy_t1.csv -o Chu_cell_type_missing5_copy1_magic.csv --mols-per-cell-min 0 --mols-per-cell-max 100000000 csv
python MAGIC.py -d Chu_cell_type_missing5_copy_t2.csv -o Chu_cell_type_missing5_copy2_magic.csv --mols-per-cell-min 0 --mols-per-cell-max 100000000 csv
python MAGIC.py -d Chu_cell_type_missing5_copy_t3.csv -o Chu_cell_type_missing5_copy3_magic.csv --mols-per-cell-min 0 --mols-per-cell-max 100000000 csv
python MAGIC.py -d Chu_cell_type_missing5_copy_t4.csv -o Chu_cell_type_missing5_copy4_magic.csv --mols-per-cell-min 0 --mols-per-cell-max 100000000 csv
python MAGIC.py -d Chu_cell_type_missing5_copy_t5.csv -o Chu_cell_type_missing5_copy5_magic.csv --mols-per-cell-min 0 --mols-per-cell-max 100000000 csv

python MAGIC.py -d Chu_time_course_missing5_copy_t1.csv -o Chu_time_course_missing5_copy1_magic.csv --mols-per-cell-min 0 --mols-per-cell-max 100000000 csv
python MAGIC.py -d Chu_time_course_missing5_copy_t2.csv -o Chu_time_course_missing5_copy2_magic.csv --mols-per-cell-min 0 --mols-per-cell-max 100000000 csv
python MAGIC.py -d Chu_time_course_missing5_copy_t3.csv -o Chu_time_course_missing5_copy3_magic.csv --mols-per-cell-min 0 --mols-per-cell-max 100000000 csv
python MAGIC.py -d Chu_time_course_missing5_copy_t4.csv -o Chu_time_course_missing5_copy4_magic.csv --mols-per-cell-min 0 --mols-per-cell-max 100000000 csv
python MAGIC.py -d Chu_time_course_missing5_copy_t5.csv -o Chu_time_course_missing5_copy5_magic.csv --mols-per-cell-min 0 --mols-per-cell-max 100000000 csv



############################################ 
############# running VIPER ################ 
############################################

######## Lasso ########

load(paste0("Chu_time_course_missing5_copy", i))
library(VIPER)

system.time(imp_ee <- VIPER(gene.expression, num = 5000, percentage.cutoff = 0.1, minbool = FALSE, alpha = 1))
save(imp_ee, file = paste0("Chu_lasso_time_course_missing5_copy", i))


load(paste0("Chu_cell_type_missing5_copy", i))
library(VIPER)

system.time(imp_ee <- VIPER(gene.expression, num = 5000, percentage.cutoff = 0.1, minbool = FALSE, alpha = 1))
save(imp_ee, file = paste0("Chu_lasso_cell_type_missing5_copy", i))



######## Elastic net ######## 


load(paste0("Chu_time_course_missing5_copy", i))
library(VIPER)
system.time(imp_el_ee <- VIPER(gene.expression, num = 5000, percentage.cutoff = 0.1, minbool = FALSE, alpha = 0.5))
save(imp_el_ee, file = paste0("Chu_el_time_course_missing5_copy", i))


load(paste0("Chu_cell_type_missing5_copy", i))
library(VIPER)
system.time(imp_el_ee <- VIPER(gene.expression, num = 5000, percentage.cutoff = 0.1, minbool = FALSE, alpha = 0.5))
save(imp_el_ee, file = paste0("Chu_el_cell_type_missing5_copy", i))



