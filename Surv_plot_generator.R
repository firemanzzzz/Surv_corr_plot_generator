library(survminer)
library(survival)
library(broom)
library(readr)
library(ggplot2)
library(survcomp)

#Load data, all of the data are nicely organized and can be downloaded from TCGA XENA data portal.
#Load gene expression
BRCA_expression <- as.data.frame(
  read_delim("H:/TCGA/TCGA_BRCA/Expression/HiSeqV2", 
             delim = "\t", escape_double = FALSE, 
             trim_ws = TRUE)
)
#Load clinical data
BRCA_survival <- as.data.frame(
  read_delim("H:/TCGA/TCGA_BRCA/Survival/survival_BRCA_survival.txt", 
             delim = "\t", escape_double = FALSE, 
             trim_ws = TRUE)
)
#Load methylaiton data (Not needed here, but can also be used for other purpose)
BRCA_methylation_NLRC5 <- as.data.frame(
  read.delim("H:/TCGA/TCGA_BRCA/Methylation/NLRC5_methylation.txt")
)

#Make the data table easier to use
#Transform axis of expression data table, facilitates z-score calculation
BRCA_expression_corre_p <- as.data.frame(t(BRCA_expression))
#Use gene name as col names
colnames(BRCA_expression_corre_p) <- BRCA_expression_corre_p[1,]
#Remove rebundant name row
BRCA_expression_corre_p <- BRCA_expression_corre_p[-1,]

#Numbers are stored as character, transform into numbers so they can be calculated
for (i in 1:20530) {
  BRCA_expression_corre_p[,i] <- sapply(BRCA_expression_corre_p[,i],as.numeric)
}

#calculate z-score
BRCA_expression_z_score <- as.data.frame(sapply(BRCA_expression_corre_p, function(BRCA_expression_corre_p) 
  (BRCA_expression_corre_p-mean(BRCA_expression_corre_p))/sd(BRCA_expression_corre_p)))
#Label each gene name  
BRCA_expression_z_score$sample <- row.names(BRCA_expression_corre_p)
#Merge expression table with clinical data
BRCA_expression_z_score <- merge(BRCA_expression_z_score, BRCA_survival, by = "sample", all.x = T)
#Check whether it is successfully merged
BRCA_expression_z_score$NLRC5


#Generate costume function: survival curve by higi-high
survival_plotter_screen <- function(x,y,z){
  #Nlrc5 and another gene are grouped by expression level, separately
  y$x.n.group[y$NLRC5 > x] <- 1
  y$x.n.group[y$NLRC5 < -x] <- 2
  y$s.unknown.group[y[,z] > x] <- 0
  y$s.unknown.group[y[,z] < -x] <- 4
  #Nlrc5 and another gene are grouped by expression level, together
  y$x.group[y$x.n.group + y$s.unknown.group == 1] <- 1
  y$x.group[y$x.n.group + y$s.unknown.group == 2] <- NA
  y$x.group[y$x.n.group + y$s.unknown.group == 6] <- 3
  y$x.group[y$x.n.group + y$s.unknown.group == 5] <- NA
  #Make sure that it is data.frame
  y <- as.data.frame(y)
  
  #Plot the figure and store in "plot" object
  plot <- ggsurvplot(
    fit = surv_fit(Surv(OS.time, OS) ~ x.group, data = y), 
    xlab = "Days", 
    ylab = "Overall survival probability",
    legend.labs = c("High-High","Low-Low"))
  
  #generate the survival analysis model
  fit <- surv_fit(Surv(OS.time, OS) ~ x.group, data = y)
  
  #And store in object "pval" and "HR"
  pval <- surv_pvalue(fit)
  HR <- hazard.ratio(y$x.group, surv.time = y$OS.time, 
                     surv.event = y$OS, na.rm = T)
  
  #display p-val and HR in plot
  plot$plot <- plot$plot +
    ggplot2::annotate("text", x = 2000, y = 0.2, label = pval$pval.txt, size = 5) +
    ggplot2::annotate("text", x = 2000, y = 0.1, label = HR$hazard.ratio, size = 5) +
    ggplot2::annotate("text", x = 2000, y = 0.3, label = colnames(y)[z], size = 5)
  return(plot)
}

#Generate costume function: survival curve by low-low
survival_plotter_screen_2 <- function(x,y,z){
  y$x.n.group[y$NLRC5 > x] <- 1
  y$x.n.group[y$NLRC5 < -x] <- 2
  y$s.unknown.group[y[,z] > x] <- 0
  y$s.unknown.group[y[,z] < -x] <- 4
  y$x.group[y$x.n.group + y$s.unknown.group == 1] <- NA
  y$x.group[y$x.n.group + y$s.unknown.group == 2] <- 2
  y$x.group[y$x.n.group + y$s.unknown.group == 6] <- NA
  y$x.group[y$x.n.group + y$s.unknown.group == 5] <- 4
  y <- as.data.frame(y)
  plot <- ggsurvplot(
    fit = surv_fit(Surv(OS.time, OS) ~ x.group, data = y), 
    xlab = "Days", 
    ylab = "Overall survival probability",
    legend.labs = c("High-Low","Low-High"))
  fit <- surv_fit(Surv(OS.time, OS) ~ x.group, data = y)
  pval <- surv_pvalue(fit)
  HR <- hazard.ratio(y$x.group, surv.time = y$OS.time, 
                     surv.event = y$OS, na.rm = T)
  plot$plot <- plot$plot +
    ggplot2::annotate("text", x = 2000, y = 0.2, label = pval$pval.txt, size = 5) +
    ggplot2::annotate("text", x = 2000, y = 0.1, label = HR$hazard.ratio, size = 5) +
    ggplot2::annotate("text", x = 2000, y = 0.3, label = colnames(y)[z], size = 5)
  return(plot)
}

survival_plotter_screen_2(0,BRCA_expression_z_score,11)

#Make a list object to store all the plot
my_surv_plots.BRCA.NLRC5 <- vector('list', 20531)
#Make all the file names for every plot
BRCA.surv.p.high.high.folder <- "H:/TCGA/TCGA_BRCA/survival_plots/high_low_vs_low_high/NLRC5/"
BRCA.surv.p.filenames <- vector('character', length = 20530)
for (i in 2:20531) {
  BRCA.surv.p.filenames[i] <- paste0(BRCA.surv.p.high.high.folder,seq_along(BRCA_expression_corre_p)[i],".svg")
}

#Generate and export all the plot
for (z  in 10154:20531) {
  tryCatch({
    message(z)
    z <- z
    my_surv_plots.BRCA.NLRC5[[z]] <- local({
      p1 <- survival_plotter_screen_2(0,BRCA_expression_z_score,z)
    })
    svg(BRCA.surv.p.filenames[z])
    plot(my_surv_plots.BRCA.NLRC5[[z]]$plot, col = z)
    dev.off()
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}


#Generate table that contain all the p-val and HR
BRCA.surv.p.val.high.high <- vector('character', length = 20530)
BRCA.surv.p.val.high.low <- vector('character', length = 20530)
BRCA.surv.HR.high.high <- vector('character', length = 20530)
BRCA.surv.HR.high.low <- vector('character', length = 20530)

survival_value_screen_1 <- function(x,y,z){
  y$x.n.group[y$NLRC5 > x] <- 1
  y$x.n.group[y$NLRC5 < -x] <- 2
  y$s.unknown.group[y[,z] > x] <- 0
  y$s.unknown.group[y[,z] < -x] <- 4
  y$x.group[y$x.n.group + y$s.unknown.group == 1] <- 1
  y$x.group[y$x.n.group + y$s.unknown.group == 2] <- NA
  y$x.group[y$x.n.group + y$s.unknown.group == 6] <- 3
  y$x.group[y$x.n.group + y$s.unknown.group == 5] <- NA
  y <- as.data.frame(y)
  fit <- surv_fit(Surv(OS.time, OS) ~ x.group, data = y)
  pval <- surv_pvalue(fit)
  HR <- hazard.ratio(y$x.group, surv.time = y$OS.time, 
                     surv.event = y$OS, na.rm = T)
  return(c(pval$pval, HR$hazard.ratio))
}

for (z  in 2:20531) {
  tryCatch({
    message(z)
    z <- z
    temp <- survival_value_screen_1(0,BRCA_expression_z_score,z)
    BRCA.surv.p.val.high.high[[z]] <- temp[1]
    BRCA.surv.HR.high.high[[z]] <- temp[2]
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

survival_value_screen_2 <- function(x,y,z){
  y$x.n.group[y$NLRC5 > x] <- 1
  y$x.n.group[y$NLRC5 < -x] <- 2
  y$s.unknown.group[y[,z] > x] <- 0
  y$s.unknown.group[y[,z] < -x] <- 4
  y$x.group[y$x.n.group + y$s.unknown.group == 1] <- NA
  y$x.group[y$x.n.group + y$s.unknown.group == 2] <- 2
  y$x.group[y$x.n.group + y$s.unknown.group == 6] <- NA
  y$x.group[y$x.n.group + y$s.unknown.group == 5] <- 4
  y <- as.data.frame(y)
  fit <- surv_fit(Surv(OS.time, OS) ~ x.group, data = y)
  pval <- surv_pvalue(fit)
  HR <- hazard.ratio(y$x.group, surv.time = y$OS.time, 
                     surv.event = y$OS, na.rm = T)
    return(c(pval$pval, HR$hazard.ratio))
}

survival_value_screen_2(0,BRCA_expression_z_score,11)

for (z  in 2:20531) {
  tryCatch({
    message(z)
    z <- z
    temp <- survival_value_screen_2(0,BRCA_expression_z_score,z)
    BRCA.surv.p.val.high.low[[z]] <- temp[1]
    BRCA.surv.HR.high.low[[z]] <- temp[2]
    },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

temp <- cbind(BRCA.surv.p.val.high.high, BRCA.surv.HR.high.high, colnames(BRCA_expression_z_score)[1:20531])
write.csv(temp, "H:/TCGA/TCGA_BRCA/survival_plots/high_high_vs_low_low/NLRC5/high_high_p_and_HR.csv")
