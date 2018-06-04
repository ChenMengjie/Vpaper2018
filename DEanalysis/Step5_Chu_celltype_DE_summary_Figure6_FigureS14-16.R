
load(paste0("cell_type_pairs_allmethods_topgenes_mean_se_DEsummary_jaccard"))
df <- total.pairs

   	aa <- as.character(df[, 1])
	df$Data <- aa
    df$Per <- df$Mean
    df$Num[df$Num==500] <- 300
    df$Num[df$Num==1000] <- 400


library(ggplot2)
theme_set(theme_bw())
p <- ggplot(data = df[df[, 2]%in%"SCDE", ], 
            aes(x = Num, 
                y = Per, 
                colour = Data)) + scale_x_continuous(breaks=c(100,200,300,400), labels=c("100", "200", "500", "1000"))+ scale_color_manual(values=c( '#CC0000', '#66B2FF', '#00CC00', '#FFAAD4', '#B266FF', '#FFD4AA', '#6363EE'))    
p <- p + geom_line() + geom_point() + facet_wrap(~ pair, scales = "free", ncol=7)  
    
p <- p + theme(axis.title.x = element_text(size = 13, face = "bold"), 
             axis.title.y = element_text(size = 13, face = "bold", angle = 90), 
             strip.text.x = element_text(size = 13, face = "bold"),
             strip.text.y = element_text(size = 13, face = "bold")) 
p <- p + theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        legend.title = element_text(size = 13, face="bold"), 
        legend.text = element_text(size = 11, face="bold"))
p <- p + ylab("Jaccard Index") + xlab("Top genes") 

ggsave(p, file="Figure6.pdf", width = 16, height = 7)





library(ggplot2)
theme_set(theme_bw())
p <- ggplot(data = df[df[, 2]%in%"DESeq2", ], 
            aes(x = Num, 
                y = Per, 
                colour = Data)) + scale_x_continuous(breaks=c(100,200,300,400), labels=c("100", "200", "500", "1000")) + scale_color_manual(values=c( '#CC0000', '#66B2FF', '#00CC00', '#FFAAD4', '#B266FF', '#FFD4AA', '#6363EE')) 
                
p <- p + geom_line() + geom_point() + facet_wrap(~ pair, scales = "free", ncol=7) #+ geom_errorbar(aes(ymin=Per-PerSE, ymax=Per+PerSE), width=.1) 
    
p <- p + theme(axis.title.x = element_text(size = 13, face = "bold"), 
             axis.title.y = element_text(size = 13, face = "bold", angle = 90), 
             strip.text.x = element_text(size = 13, face = "bold"),
             strip.text.y = element_text(size = 13, face = "bold")) 
p <- p + theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        legend.title = element_text(size = 13, face="bold"), 
        legend.text = element_text(size = 11, face="bold"))
p <- p + ylab("Jaccard Index") + xlab("Top genes") 

ggsave(p, file="FigureS14.pdf", width = 16, height = 7)



library(ggplot2)
theme_set(theme_bw())
p <- ggplot(data = df[df[, 2]%in%"edgeR_LRT", ], 
            aes(x = Num, 
                y = Per, 
                colour = Data)) + scale_x_continuous(breaks=c(100,200,300,400), labels=c("100", "200", "500", "1000"))+ scale_color_manual(values=c( '#CC0000', '#66B2FF', '#00CC00', '#FFAAD4', '#B266FF', '#FFD4AA', '#6363EE'))   
p <- p + geom_line() + geom_point() + facet_wrap(~ pair, scales = "free", ncol=7)  
    
p <- p + theme(axis.title.x = element_text(size = 13, face = "bold"), 
             axis.title.y = element_text(size = 13, face = "bold", angle = 90), 
             strip.text.x = element_text(size = 13, face = "bold"),
             strip.text.y = element_text(size = 13, face = "bold")) 
p <- p + theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        legend.title = element_text(size = 13, face="bold"), 
        legend.text = element_text(size = 11, face="bold"))
p <- p + ylab("Jaccard Index") + xlab("Top genes") 

ggsave(p, file="FigureS15.pdf", width = 16, height = 7)




library(ggplot2)
theme_set(theme_bw())
p <- ggplot(data = df[df[, 2]%in%"edgeR_QLF", ], 
            aes(x = Num, 
                y = Per, 
                colour = Data)) + scale_x_continuous(breaks=c(100,200,300,400), labels=c("100", "200", "500", "1000"))+ scale_color_manual(values=c( '#CC0000', '#66B2FF', '#00CC00', '#FFAAD4', '#B266FF', '#FFD4AA', '#6363EE'))   
p <- p + geom_line() + geom_point() + facet_wrap(~ pair, scales = "free", ncol=7)  
    
p <- p + theme(axis.title.x = element_text(size = 13, face = "bold"), 
             axis.title.y = element_text(size = 13, face = "bold", angle = 90), 
             strip.text.x = element_text(size = 13, face = "bold"),
             strip.text.y = element_text(size = 13, face = "bold")) 
p <- p + theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        legend.title = element_text(size = 13, face="bold"), 
        legend.text = element_text(size = 11, face="bold"))
p <- p + ylab("Jaccard Index") + xlab("Top genes") 

ggsave(p, file="FigureS16.pdf", width = 16, height = 7)
