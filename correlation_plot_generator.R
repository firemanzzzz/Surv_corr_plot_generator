library(DESeq2)
library(survminer)
library(survival)
library(broom)
library(readr)
library(ggpubr)

#BRCA as example

BRCA_expression_corre_p <- as.data.frame(t(BRCA_expression))
colnames(BRCA_expression_corre_p) <- BRCA_expression_corre_p[1,]
BRCA_expression_corre_p <- BRCA_expression_corre_p[-1,]

for (i in 1:20530) {
  BRCA_expression_corre_p[,i] <- sapply(BRCA_expression_corre_p[,i],as.numeric)
}

mypval.BRCA.NLRC5 <- vector('numeric', ncol(BRCA_expression_corre_p))
for (i in 1:20530) {
  mypval.BRCA.NLRC5[i] <- cor(BRCA_expression_corre_p$NLRC5, BRCA_expression_corre_p[,i], method = "pearson")
}
mypval.BRCA.NLRC5 <- cbind(mypval.BRCA.NLRC5,1:20530,colnames(BRCA_expression_corre_p))




myplots.BRCA.NLRC5 <- vector('list', ncol(BRCA_expression_corre_p))

for (i in seq_along(BRCA_expression_corre_p)) {
  message(i)
  myplots.BRCA[[i]] <- local({
    i <- i
    p1 <- ggplot(data=BRCA_expression_corre_p, aes(x=NLRC5,y=BRCA_expression_corre_p[,i]))+
      geom_point(color = "DarkBlue")+
      theme_bw()+
      geom_quantile(quantiles = 0.5,size = 2, color = "red") +
      ggtitle(colnames(BRCA_expression_corre_p)[i]) +
      xlab("NLRC5 expression") + ylab(colnames(BRCA_expression_corre_p)[i]) +
      theme(plot.title = element_text(hjust = 0.5))+
      theme(panel.border = element_rect(colour = "black", fill = NA, size= 0.5),    
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            axis.text=element_text(size=16),
            axis.title=element_text(size=18),
            plot.title = element_text(size=20,hjust = 0.5),
            legend.key.size = unit(1, 'cm'),
            legend.title = element_text(size=16),
            legend.text = element_text(size=16))
    p1 <- p1+stat_cor(method="pearson",size=6)
  })
}

BRCA.p.folder <- "I:/TCGA/TCGA_BRCA/correlation_plots/NLRC5/"
BRCA.p.filenames <- vector('character', length = ncol(BRCA_expression_corre_p))

for (i in seq_along(BRCA_expression_corre_p)) {
  BRCA.p.filenames[i] <- paste0(BRCA.p.folder,seq_along(BRCA_expression_corre_p)[i],".svg")
}

i <- 1

for (i in seq_along(BRCA_expression_corre_p)) {
  message(i)
  svg(BRCA.p.filenames[i])
  plot(myplots.BRCA[[i]], col = i)
  dev.off()
}

write.csv(mypval.BRCA.NLRC5, file = "H:/TCGA/TCGA_BRCA/correlation_plots/NLRC5/pval.NLRC5.csv")




#HLA-A
mypval.BRCA.HLAA <- vector('numeric', ncol(BRCA_expression_corre_p))

for (i in 1:20530) {
  mypval.BRCA.HLAA[i] <- cor(BRCA_expression_corre_p$`HLA-A`, BRCA_expression_corre_p[,i], method = "pearson")
}

mypval.BRCA.HLAA <- cbind(mypval.BRCA.HLAA,1:20530,colnames(BRCA_expression_corre_p))

myplots.BRCA.HLAA <- vector('list', ncol(BRCA_expression_corre_p))

for (i in seq_along(BRCA_expression_corre_p)) {
  message(i)
  myplots.BRCA.HLAA[[i]] <- local({
    i <- i
    p1 <- ggplot(data=BRCA_expression_corre_p, aes(x=BRCA_expression_corre_p[,7750],y=BRCA_expression_corre_p[,i]))+
      geom_point(color = "DarkBlue")+
      theme_bw()+
      geom_quantile(quantiles = 0.5,size = 2, color = "red") +
      ggtitle(colnames(BRCA_expression_corre_p)[i]) +
      xlab("HLA-A") + ylab(colnames(BRCA_expression_corre_p)[i]) +
      theme(plot.title = element_text(hjust = 0.5))+
      theme(panel.border = element_rect(colour = "black", fill = NA, size= 0.5),    
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            axis.text=element_text(size=16),
            axis.title=element_text(size=18),
            plot.title = element_text(size=20,hjust = 0.5),
            legend.key.size = unit(1, 'cm'),
            legend.title = element_text(size=16),
            legend.text = element_text(size=16))
    p1 <- p1+stat_cor(method="pearson",size=6)
  })
}

BRCA.p.folder <- "I:/TCGA/TCGA_BRCA/correlation_plots/HLAA/"
BRCA.p.filenames <- vector('character', length = ncol(BRCA_expression_corre_p))

for (i in seq_along(BRCA_expression_corre_p)) {
  BRCA.p.filenames [i] <- paste0(BRCA.p.floder,colnames(BRCA_expression_corre_p)[i],".svg")
}

for (i in seq_along(BRCA_expression_corre_p)) {
    svg(BRCA.p.filenames[i])
    set.seed(i)
    plot(myplots.BRCA.HLAA[[i]], col = i)
    dev.off()
}

write.csv(mypval.BRCA.HLAA, file = "H:/TCGA/TCGA_BRCA/correlation_plots/HLAA/pval.HLAA.csv")
