###########################################
########## creating two datasets ##########
###########################################
xx <- read.csv("GSE75748_sc_cell_type_ec.csv", as.is = T)
rownames(xx) <- xx[, 1]
xx <- xx[, -1]
aa <- gsub("\\_[^<>]*", "", colnames(xx))
celltype <- unique(aa)
num <- length(unique(aa))

load("GSE75748_sc_cell_type.rda")


for(i in 1:10){
	set.seed(i*500+80224523)
	set1 <- NULL
	set2 <- NULL
	for(j in 1:num){
		flag <- which(aa%in%celltype[j])
		nsamples <- length(flag)
		selected <- sample(flag, floor(nsamples/2))
		set1 <- c(set1, selected)
		set2 <- c(set2, flag[!flag%in%selected])	
	}
    gene.expression1 <- gene.expression[, set1]
	gene.expression2 <- gene.expression[, set2]
	
	save(set1, set2, gene.expression1, gene.expression2, file = paste0("Chu_cell_type_twosets_copy", i))
    write.csv(gene.expression1, paste0("Chu_cell_type_twosets_copy", i, "_1.csv"))
	write.csv(t(gene.expression1), paste0("Chu_cell_type_twosets_copy", i, "_t1.csv"))
	write.csv(gene.expression2, paste0("Chu_cell_type_twosets_copy", i, "_2.csv"))
	write.csv(t(gene.expression2), paste0("Chu_cell_type_twosets_copy", i, "_t2.csv"))
}



###############################################  
############# running scImpute ################ 
############################################### 
i <- 1 #2,3,4,5,6,7,8,9,10

library(scImpute)

system.time(scimpute(count_path = paste0("Chu_cell_type_twosets_copy", i, "_1.csv"),
  infile = "csv",           # format of input file
  outfile = "csv",          # format of output file
  out_dir = paste0("scImpute/cell_type/twosets_copy", i, "_1"),  # full path to output directory
  drop_thre = 0.5,          # threshold set on dropout probability
  Kcluster = 7,
  ncores = 1))
  
system.time(scimpute(count_path = paste0("Chu_cell_type_twosets_copy", i, "_2.csv"),
  infile = "csv",           # format of input file
  outfile = "csv",          # format of output file
  out_dir = paste0("scImpute/cell_type/twosets_copy", i, "_2"),  # full path to output directory
  drop_thre = 0.5,          # threshold set on dropout probability
  Kcluster = 7,
  ncores = 1))
    


###############################################  
############# running drImpute ################ 
############################################### 
i <- 1 #2,3,4,5,6,7,8,9,10

load(paste0("Chu_cell_type_twosets_copy", i))
library(DrImpute)
exdata <- preprocessSC(gene.expression1)
sf <- apply(exdata, 2, mean)
npX <- t(t(exdata)/sf) 
lnpX <- log(npX+1)
system.time(lnpX_imp1 <- DrImpute(lnpX))

exdata <- preprocessSC(gene.expression2)
sf <- apply(exdata, 2, mean)
npX <- t(t(exdata)/sf) 
lnpX <- log(npX+1)
system.time(lnpX_imp2 <- DrImpute(lnpX))
save(lnpX_imp1, lnpX_imp2, file = paste0("Chu_drImpute_cell_type_twosets_tc_copy", i))



###############################################  
############# running SAVER ################### 
############################################### 
i <- 1 #2,3,4,5,6,7,8,9,10
load(paste0("Chu_cell_type_twosets_copy", i))

library(SAVER)
system.time(saver_imp1 <- saver(gene.expression1))
system.time(saver_imp2 <- saver(gene.expression2))
save(saver_imp1, saver_imp2, file = paste0("Chu_SAVER_cell_type_twosets_tc_copy", i))


###############################################  
############# running MAGIC ################### 
############################################### 

python MAGIC.py -d Chu_cell_type_twosets_copy1_t1.csv -o Chu_cell_type_twosets_copy1_t1_magic.csv --mols-per-cell-min 0 --mols-per-cell-max 100000000 csv
python MAGIC.py -d Chu_cell_type_twosets_copy2_t1.csv -o Chu_cell_type_twosets_copy2_t1_magic.csv --mols-per-cell-min 0 --mols-per-cell-max 100000000 csv
python MAGIC.py -d Chu_cell_type_twosets_copy3_t1.csv -o Chu_cell_type_twosets_copy3_t1_magic.csv --mols-per-cell-min 0 --mols-per-cell-max 100000000 csv
python MAGIC.py -d Chu_cell_type_twosets_copy4_t1.csv -o Chu_cell_type_twosets_copy4_t1_magic.csv --mols-per-cell-min 0 --mols-per-cell-max 100000000 csv
python MAGIC.py -d Chu_cell_type_twosets_copy5_t1.csv -o Chu_cell_type_twosets_copy5_t1_magic.csv --mols-per-cell-min 0 --mols-per-cell-max 100000000 csv
python MAGIC.py -d Chu_cell_type_twosets_copy6_t1.csv -o Chu_cell_type_twosets_copy6_t1_magic.csv --mols-per-cell-min 0 --mols-per-cell-max 100000000 csv
python MAGIC.py -d Chu_cell_type_twosets_copy7_t1.csv -o Chu_cell_type_twosets_copy7_t1_magic.csv --mols-per-cell-min 0 --mols-per-cell-max 100000000 csv
python MAGIC.py -d Chu_cell_type_twosets_copy8_t1.csv -o Chu_cell_type_twosets_copy8_t1_magic.csv --mols-per-cell-min 0 --mols-per-cell-max 100000000 csv
python MAGIC.py -d Chu_cell_type_twosets_copy9_t1.csv -o Chu_cell_type_twosets_copy9_t1_magic.csv --mols-per-cell-min 0 --mols-per-cell-max 100000000 csv
python MAGIC.py -d Chu_cell_type_twosets_copy10_t1.csv -o Chu_cell_type_twosets_copy10_t1_magic.csv --mols-per-cell-min 0 --mols-per-cell-max 100000000 csv


python MAGIC.py -d Chu_cell_type_twosets_copy1_t2.csv -o Chu_cell_type_twosets_copy1_t2_magic.csv --mols-per-cell-min 0 --mols-per-cell-max 100000000 csv
python MAGIC.py -d Chu_cell_type_twosets_copy2_t2.csv -o Chu_cell_type_twosets_copy2_t2_magic.csv --mols-per-cell-min 0 --mols-per-cell-max 100000000 csv
python MAGIC.py -d Chu_cell_type_twosets_copy3_t2.csv -o Chu_cell_type_twosets_copy3_t2_magic.csv --mols-per-cell-min 0 --mols-per-cell-max 100000000 csv
python MAGIC.py -d Chu_cell_type_twosets_copy4_t2.csv -o Chu_cell_type_twosets_copy4_t2_magic.csv --mols-per-cell-min 0 --mols-per-cell-max 100000000 csv
python MAGIC.py -d Chu_cell_type_twosets_copy5_t2.csv -o Chu_cell_type_twosets_copy5_t2_magic.csv --mols-per-cell-min 0 --mols-per-cell-max 100000000 csv
python MAGIC.py -d Chu_cell_type_twosets_copy6_t2.csv -o Chu_cell_type_twosets_copy6_t2_magic.csv --mols-per-cell-min 0 --mols-per-cell-max 100000000 csv
python MAGIC.py -d Chu_cell_type_twosets_copy7_t2.csv -o Chu_cell_type_twosets_copy7_t2_magic.csv --mols-per-cell-min 0 --mols-per-cell-max 100000000 csv
python MAGIC.py -d Chu_cell_type_twosets_copy8_t2.csv -o Chu_cell_type_twosets_copy8_t2_magic.csv --mols-per-cell-min 0 --mols-per-cell-max 100000000 csv
python MAGIC.py -d Chu_cell_type_twosets_copy9_t2.csv -o Chu_cell_type_twosets_copy9_t2_magic.csv --mols-per-cell-min 0 --mols-per-cell-max 100000000 csv
python MAGIC.py -d Chu_cell_type_twosets_copy10_t2.csv -o Chu_cell_type_twosets_copy10_t2_magic.csv --mols-per-cell-min 0 --mols-per-cell-max 100000000 csv



############################################ 
############# running VIPER ################ 
############################################

######## Lasso ########

i <- 1 #2,3,4,5,6,7,8,9,10

load(paste0("Chu_cell_type_twosets_copy", i))

library(VIPER)

system.time(imp_ee1 <- VIPER(gene.expression1, num = 5000, percentage.cutoff = 0.1, minbool = FALSE, alpha = 1))
system.time(imp_ee2 <- VIPER(gene.expression2, num = 5000, percentage.cutoff = 0.1, minbool = FALSE, alpha = 1))
save(imp_ee1, imp_ee2, file = paste0("Chu_cell_type_twosets_copy_lasso_copy", i))




######## Elastic net ######## 

i <- 1 #2,3,4,5,6,7,8,9,10

load(paste0("Chu_cell_type_twosets_copy", i)) 
library(VIPER)
system.time(imp_el_ee1 <- VIPER(gene.expression1, num = 5000, percentage.cutoff = 0.1, minbool = FALSE, alpha = 0.5))
system.time(imp_el_ee2 <- VIPER(gene.expression2, num = 5000, percentage.cutoff = 0.1, minbool = FALSE, alpha = 0.5))
save(imp_el_ee1, imp_el_ee2, file = paste0("Chu_cell_type_twosets_copy_el_copy", i))


