
load("Chu_cell_type_downsampling_comparison_summary_CorL1L2")
correlation.summary1[, -c(1:2)] <- apply(correlation.summary1[, -c(1:2)], 2, function(x){round(x, 3)})

aa <- correlation.summary1

load("Chu_time_course_downsampling_comparison_summary_CorL1L2")
correlation.summary1[, -c(1:2)] <- apply(correlation.summary1[, -c(1:2)], 2, function(x){round(x, 3)})

aa2 <- correlation.summary1

load("Grun_downsampling_comparison_summary_CorL1L2")
correlation.summary1[, -c(1:2)] <- apply(correlation.summary1[, -c(1:2)], 2, function(x){round(x, 3)})

aa3 <- correlation.summary1 
  
load("Shalek_downsampling_comparison_summary_CorL1L2")
correlation.summary1[, -c(1:2)] <- apply(correlation.summary1[, -c(1:2)], 2, function(x){round(x, 3)})

aa4 <- correlation.summary1
  
aa.total <- rbind(aa, aa2, aa3, aa4)
correlation.total <- aa.total





load("Chu_cell_type_downsampling_comparison_summary_CorL1L2")
correlation2.summary1[, -c(1:2)] <- apply(correlation2.summary1[, -c(1:2)], 2, function(x){round(x, 3)})

aa <- correlation2.summary1

load("Chu_time_course_downsampling_comparison_summary_CorL1L2")
correlation2.summary1[, -c(1:2)] <- apply(correlation2.summary1[, -c(1:2)], 2, function(x){round(x, 3)})

aa2 <- correlation2.summary1

load("Grun_downsampling_comparison_summary_CorL1L2")
correlation2.summary1[, -c(1:2)] <- apply(correlation2.summary1[, -c(1:2)], 2, function(x){round(x, 3)})

aa3 <- correlation2.summary1 
  
load("Shalek_downsampling_comparison_summary_CorL1L2")
correlation2.summary1[, -c(1:2)] <- apply(correlation2.summary1[, -c(1:2)], 2, function(x){round(x, 3)})

aa4 <- correlation2.summary1
  
aa.total <- rbind(aa, aa2, aa3, aa4)
correlation2.total <- aa.total






load("Chu_cell_type_downsampling_comparison_summary_CorL1L2")
L1.summary1[, -c(1:2)] <- apply(L1.summary1[, -c(1:2)], 2, function(x){round(x, 3)})

aa <- L1.summary1

load("Chu_time_course_downsampling_comparison_summary_CorL1L2")
L1.summary1[, -c(1:2)] <- apply(L1.summary1[, -c(1:2)], 2, function(x){round(x, 3)})

aa2 <- L1.summary1 


load("Grun_downsampling_comparison_summary_CorL1L2")
L1.summary1[, -c(1:2)] <- apply(L1.summary1[, -c(1:2)], 2, function(x){round(x, 3)})

aa3 <- L1.summary1  
  

load("Shalek_downsampling_comparison_summary_CorL1L2")
L1.summary1[, -c(1:2)] <- apply(L1.summary1[, -c(1:2)], 2, function(x){round(x, 3)})

aa4 <- L1.summary1   
  
aa.total <- rbind(aa, aa2, aa3, aa4)
L1.total <- aa.total



load("Chu_cell_type_downsampling_comparison_summary_CorL1L2")
L1.summary2[, -c(1:2)] <- apply(L1.summary2[, -c(1:2)], 2, function(x){round(x, 3)})

aa <- L1.summary2

load("Chu_time_course_downsampling_comparison_summary_CorL1L2")
L1.summary2[, -c(1:2)] <- apply(L1.summary2[, -c(1:2)], 2, function(x){round(x, 3)})

aa2 <- L1.summary2 


load("Grun_downsampling_comparison_summary_CorL1L2")
L1.summary2[, -c(1:2)] <- apply(L1.summary2[, -c(1:2)], 2, function(x){round(x, 3)})

aa3 <- L1.summary2  
  

load("Shalek_downsampling_comparison_summary_CorL1L2")
L1.summary2[, -c(1:2)] <- apply(L1.summary2[, -c(1:2)], 2, function(x){round(x, 3)})

aa4 <- L1.summary2   

aa.total <- rbind(aa, aa2, aa3, aa4)
L1_2.total <- aa.total


save(correlation.total, correlation2.total, L1.total, L1_2.total, file = "downsampling.corL1L2.for.plot")


load("downsampling.corL1L2.for.plot")
res <- NULL
for(i in 1:6){
	tt <- cbind(L1_2.total[, 1:2], correlation.total[, 2+i], L1_2.total[, 2+i], colnames(correlation.total)[2+i])
	colnames(tt) <- c("dataset", "downsampling", "correlation", "coverage", "method")
	res <- rbind(res, tt)
	
}

df <- as.data.frame(res)
bb <- as.character(df[, 1])
bb[bb%in%"Chu_cell_type" ] <- "Cell Type"
bb[bb%in%"Chu_time_course"] <- "Time Course"
df$Data <- bb

aa <- as.character(df[, 5])
aa[aa%in%"lnpX"] <- "DrImpute"
aa[aa%in%"el_ee"] <- "VIPER-Elastic Net"
aa[aa%in%"ee"] <- "VIPER-Lasso"
aa[aa%in%"scimpute"] <- "scImpute"
aa[aa%in%"saver"] <- "SAVER"
aa[aa%in%"magic"] <- "MAGIC"
df$Method <- aa
df$Downsampling <- as.factor(df$downsampling)

library(ggplot2)
theme_set(theme_bw())
p <- ggplot(data = df, 
            aes(y = coverage, 
                x = correlation, 
                colour = Method, shape = Downsampling))  + scale_color_manual(values=c('#66B2FF', '#00CC00', '#FFAAD4', '#B266FF', '#FFD4AA', '#6363EE'))  #+ scale_x_continuous(breaks=c(0.2,0.4,0.6,0.8))
                
p <- p + geom_line(size=2) + geom_point(size=5) + facet_wrap(~ Data, scales = "free", ncol=2) 
    
p <- p + theme(axis.title.x = element_text(size = 20, face = "bold"), 
             axis.title.y = element_text(size = 20, face = "bold", angle = 90), 
             axis.text.x = element_text(size = 15), 
             axis.text.y = element_text(size = 15), 
             strip.text.x = element_text(size = 20, face = "bold"),
             strip.text.y = element_text(size = 20, face = "bold")) 
p <- p + theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        legend.title = element_text(size = 15, face="bold"), 
        legend.text = element_text(size = 15, face="bold"))
p <- p + xlab("Correlation (Drop-out entries)") + ylab("L1 Loss (Down-sampling entries)") 

ggsave(p, file="Figure3A_downsampling.corL1L2_sum.pdf", width = 10, height = 8)




res <- NULL
for(i in 1:6){
	tt <- cbind(L1_2.total[, 1:2], L1.total[, 2+i], L1_2.total[, 2+i], colnames(correlation.total)[2+i])
	colnames(tt) <- c("dataset", "downsampling", "correlation", "coverage", "method")
	res <- rbind(res, tt)
	
}

df <- as.data.frame(res)
bb <- as.character(df[, 1])
bb[bb%in%"Chu_cell_type" ] <- "Cell Type"
bb[bb%in%"Chu_time_course"] <- "Time Course"
df$Data <- bb

aa <- as.character(df[, 5])
aa[aa%in%"lnpX"] <- "DrImpute"
aa[aa%in%"el_ee"] <- "VIPER-Elastic Net"
aa[aa%in%"ee"] <- "VIPER-Lasso"
aa[aa%in%"scimpute"] <- "scImpute"
aa[aa%in%"saver"] <- "SAVER"
aa[aa%in%"magic"] <- "MAGIC"
df$Method <- aa
df$Downsampling <- as.factor(df$downsampling)

library(ggplot2)
theme_set(theme_bw())
p <- ggplot(data = df, 
            aes(y = coverage, 
                x = correlation, 
                colour = Method, shape = Downsampling))  + scale_color_manual(values=c('#66B2FF', '#00CC00', '#FFAAD4', '#B266FF', '#FFD4AA', '#6363EE'))  
                
p <- p + geom_line(size=2) + geom_point(size=5) + facet_wrap(~ Data, scales = "free", ncol=2) 
    
p <- p + theme(axis.title.x = element_text(size = 20, face = "bold"), 
             axis.title.y = element_text(size = 20, face = "bold", angle = 90), 
             axis.text.x = element_text(size = 15), 
             axis.text.y = element_text(size = 15), 
             strip.text.x = element_text(size = 20, face = "bold"),
             strip.text.y = element_text(size = 20, face = "bold")) 
p <- p + theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        legend.title = element_text(size = 15, face="bold"), 
        legend.text = element_text(size = 15, face="bold"))
p <- p + xlab("L1 Loss (Drop-out entries)") + ylab("L1 Loss (Down-sampling entries)") 

ggsave(p, file="FigureS6A_downsampling_L1_.corL1L2_sum.pdf", width = 10, height = 8)

