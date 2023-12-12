protein_name <- proteome_X$X
protein_name <- substr(protein_name, 1, nchar(protein_name)-2)
write.csv(protein_name, file = "H:/CPTAC/protein_name.csv")

proteome_X_corr <- as.data.frame(t(proteome_X))
View(proteome_X_corr)
colnames(proteome_X_corr) <- proteome_X_corr[1,]
proteome_X_corr <- proteome_X_corr[-1,]

for (i in 1:20742) {
  proteome_X_corr[,i] <- sapply(proteome_X_corr[,i], as.numeric)
}

get.colnumber.subset <- function(gene.name){
  temp.name <- protein_gene_name_conver[which(protein_gene_name_conver$Gene.name == gene.name),3]
  temp.times <-  length(temp.name)
  for (i in 1:temp.times) {
    temp.number <- grep(temp.name[i],colnames(proteome_X_subset))
    print(temp.number)
  }
}

protein_gene_name_conver <- read.delim("H:/CPTAC/mart_export (1).txt")

proteome_X_subset <- as.data.frame(t(proteome_X))
colnames(proteome_X_subset) <- proteome_X_subset[1,]
proteome_X_subset <- proteome_X_subset[-1,]
for (i in 1:20742) {
  proteome_X_subset[,i] <- sapply(proteome_X_subset[,i],as.numeric)
}
proteome_X_subset$Proteome_Sample_ID <- row.names(proteome_X_subset)
proteome_X_subset <- merge(proteome_X_subset, CPTAC_meta, by = "Proteome_Sample_ID", all = T)

unique(CPTAC_meta$cohort)
proteome_X_subset_BRCA <- proteome_X_subset[which(proteome_X_subset$cohort == "BRCA"),]
proteome_X_subset_ccRCC <- proteome_X_subset[which(proteome_X_subset$cohort == "ccRCC"),]
proteome_X_subset_COAD <- proteome_X_subset[which(proteome_X_subset$cohort == "COAD"),]
proteome_X_subset_GBM <- proteome_X_subset[which(proteome_X_subset$cohort == "GBM"),]
proteome_X_subset_LSCC <- proteome_X_subset[which(proteome_X_subset$cohort == "LSCC"),]
proteome_X_subset_LUAD <- proteome_X_subset[which(proteome_X_subset$cohort == "LUAD"),]
proteome_X_subset_HGSC <- proteome_X_subset[which(proteome_X_subset$cohort == "HGSC"),]
proteome_X_subset_PDAC <- proteome_X_subset[which(proteome_X_subset$cohort == "PDAC"),]
proteome_X_subset_UCEC <- proteome_X_subset[which(proteome_X_subset$cohort == "UCEC"),]
proteome_X_subset_HNSCC <- proteome_X_subset[which(proteome_X_subset$cohort == "HNSCC"),]

subset.correl.gen <- function(protein.number){

mypval.temp <- data.frame(matrix(ncol = 31, nrow = 0))
colnames(mypval.temp) <- c("RefSeq.peptide.ID","BRCA.r", "BRCA.p", "BRCA.n",
                            "ccRCC.r", "ccRCC.p", "ccRCC.n",
                            "COAD.r", "COAD.p", "COAD.n",
                            "GBM.r", "GBM.p", "GBM.n",
                            "LSCC.r", "LSCC.p", "LSCC.n",
                            "LUAD.r", "LUAD.p", "LUAD.n",
                            "HGSC.r", "HGSC.p", "HGSC.n",
                            "PDAC.r", "PDAC.p", "PDAC.n",
                            "UCEC.r", "UCEC.p", "UCEC.n",
                            "HNSCC.r", "HNSCC.p", "HNSCC.n")
message(1)
for (i in 1:20742) {
  tryCatch({
    
    mypval.temp[i,1] <- colnames(proteome_X_subset)[i+1]
    temp <- list()
    temp <- cor.test(proteome_X_subset_BRCA[,protein.number], proteome_X_subset_BRCA[,i+1], 
                     na.action=na.omit, method = "pearson")
    mypval.temp[i,2] <- temp$estimate
    mypval.temp[i,3] <- temp$p.value
    mypval.temp[i,4] <- temp$parameter + 1
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}   
message(2)
for (i in 1:20742) {
  tryCatch({
    
    mypval.temp[i,1] <- colnames(proteome_X_subset)[i+1]
    temp <- list()
    temp <- cor.test(proteome_X_subset_ccRCC[,protein.number], proteome_X_subset_ccRCC[,i+1], 
                     na.action=na.omit, method = "pearson")
    mypval.temp[i,5] <- temp$estimate
    mypval.temp[i,6] <- temp$p.value
    mypval.temp[i,7] <- temp$parameter + 1
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}    
message(3)
for (i in 1:20742) {
  tryCatch({
    
    mypval.temp[i,1] <- colnames(proteome_X_subset)[i+1]
    temp <- list()
    temp <- cor.test(proteome_X_subset_COAD[,protein.number], proteome_X_subset_COAD[,i+1], 
                     na.action=na.omit, method = "pearson")
    mypval.temp[i,8] <- temp$estimate
    mypval.temp[i,9] <- temp$p.value
    mypval.temp[i,10] <- temp$parameter + 1
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
} 
message(4)
for (i in 1:20742) {
  tryCatch({
    
    mypval.temp[i,1] <- colnames(proteome_X_subset)[i+1]
    temp <- list()
    temp <- cor.test(proteome_X_subset_GBM[,protein.number], proteome_X_subset_GBM[,i+1], 
                     na.action=na.omit, method = "pearson")
    mypval.temp[i,11] <- temp$estimate
    mypval.temp[i,12] <- temp$p.value
    mypval.temp[i,13] <- temp$parameter + 1
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
message(5)
for (i in 1:20742) {
  tryCatch({
    
    mypval.temp[i,1] <- colnames(proteome_X_subset)[i+1]
    temp <- list()
    temp <- cor.test(proteome_X_subset_LSCC[,protein.number], proteome_X_subset_LSCC[,i+1], 
                     na.action=na.omit, method = "pearson")
    mypval.temp[i,14] <- temp$estimate
    mypval.temp[i,15] <- temp$p.value
    mypval.temp[i,16] <- temp$parameter + 1
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
message(6)
for (i in 1:20742) {
  tryCatch({
    
    mypval.temp[i,1] <- colnames(proteome_X_subset)[i+1]
    temp <- list()
    temp <- cor.test(proteome_X_subset_LUAD[,protein.number], proteome_X_subset_LUAD[,i+1], 
                     na.action=na.omit, method = "pearson")
    mypval.temp[i,17] <- temp$estimate
    mypval.temp[i,18] <- temp$p.value
    mypval.temp[i,19] <- temp$parameter + 1
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}    
message(7) 
for (i in 1:20742) {
  tryCatch({
    
    mypval.temp[i,1] <- colnames(proteome_X_subset)[i+1]
    temp <- list()
    temp <- cor.test(proteome_X_subset_HGSC[,protein.number], proteome_X_subset_HGSC[,i+1], 
                     na.action=na.omit, method = "pearson")
    mypval.temp[i,20] <- temp$estimate
    mypval.temp[i,21] <- temp$p.value
    mypval.temp[i,22] <- temp$parameter + 1
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
} 
message(8)
for (i in 1:20742) {
  tryCatch({
    
    mypval.temp[i,1] <- colnames(proteome_X_subset)[i+1]
    temp <- list()
    temp <- cor.test(proteome_X_subset_PDAC[,protein.number], proteome_X_subset_PDAC[,i+1], 
                     na.action=na.omit, method = "pearson")
    mypval.temp[i,23] <- temp$estimate
    mypval.temp[i,24] <- temp$p.value
    mypval.temp[i,25] <- temp$parameter + 1
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}    
message(9)
for (i in 1:20742) {
  tryCatch({
    
    mypval.temp[i,1] <- colnames(proteome_X_subset)[i+1]
    temp <- list()
    temp <- cor.test(proteome_X_subset_UCEC[,protein.number], proteome_X_subset_UCEC[,i+1], 
                     na.action=na.omit, method = "pearson")
    mypval.temp[i,26] <- temp$estimate
    mypval.temp[i,27] <- temp$p.value
    mypval.temp[i,28] <- temp$parameter + 1
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}  
message(10)
for (i in 1:20742) {
  tryCatch({
    
    mypval.temp[i,1] <- colnames(proteome_X_subset)[i+1]
    temp <- list()
    temp <- cor.test(proteome_X_subset_HNSCC[,protein.number], proteome_X_subset_HNSCC[,i+1], 
                     na.action=na.omit, method = "pearson")
    mypval.temp[i,29] <- temp$estimate
    mypval.temp[i,30] <- temp$p.value
    mypval.temp[i,31] <- temp$parameter + 1
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}   

mypval.temp$RefSeq.peptide.ID <- substr(mypval.temp$RefSeq.peptide.ID, 1, nchar(mypval.temp$RefSeq.peptide.ID)-2)
mypval.temp <- merge(mypval.temp, protein_gene_name_conver, by = "RefSeq.peptide.ID", all.x = T)    
return(mypval.temp)
}

get.colnumber <- function(gene.name){
  temp.name <- protein_gene_name_conver[which(protein_gene_name_conver$Gene.name == gene.name),3]
  temp.times <-  length(temp.name)
  for (i in 1:temp.times) {
    temp.number <- grep(temp.name[i],colnames(proteome_X_subset))
    print(temp.number)
  }
}

get.colnumber.subset("NLRC5")
get.colnumber.subset("SUPT7L")[1]
get.colnumber.subset("FBXO11")
get.colnumber.subset("WDR11")
get.colnumber.subset("NACA")
get.colnumber.subset("NFYA")
get.colnumber.subset("NFYB")
get.colnumber.subset("NFYC")
get.colnumber.subset("HDAC1")
get.colnumber.subset("HDAC2")
get.colnumber.subset("HDAC3")
get.colnumber.subset("HDAC4")
get.colnumber.subset("HDAC5")
get.colnumber.subset("HDAC6")
get.colnumber.subset("HDAC7")
get.colnumber.subset("HDAC8")
get.colnumber.subset("HDAC9")
get.colnumber.subset("HDAC10")
get.colnumber.subset("HDAC11")
get.colnumber.subset("NOD2")

mypval.NACA <- subset.correl.gen(7604)
write.csv(mypval.NACA, file = "H:/CPTAC/CPTAC/mypval.NACA.csv")
mypval.NFYA <- subset.correl.gen(10806)
write.csv(mypval.NFYA, file = "H:/CPTAC/CPTAC/mypval.NFYA.csv")
mypval.NFYB <- subset.correl.gen(9995)
write.csv(mypval.NFYB, file = "H:/CPTAC/CPTAC/mypval.NFYB.csv")
mypval.NFYC <- subset.correl.gen(8062)
write.csv(mypval.NFYC, file = "H:/CPTAC/CPTAC/mypval.NFYC.csv")
mypval.HDAC1 <- subset.correl.gen(3421)
write.csv(mypval.HDAC1, file = "H:/CPTAC/CPTAC/mypval.HDAC1.csv")
mypval.HDAC2 <- subset.correl.gen(3422)
write.csv(mypval.HDAC2, file = "H:/CPTAC/CPTAC/mypval.HDAC2.csv")
mypval.HDAC3 <- subset.correl.gen(6358)
write.csv(mypval.HDAC3, file = "H:/CPTAC/CPTAC/mypval.HDAC3.csv")
mypval.HDAC4 <- subset.correl.gen(3053)
write.csv(mypval.HDAC4, file = "H:/CPTAC/CPTAC/mypval.HDAC4.csv")
mypval.HDAC5 <- subset.correl.gen(3054)
write.csv(mypval.HDAC5, file = "H:/CPTAC/CPTAC/mypval.HDAC5.csv")
mypval.HDAC6 <- subset.correl.gen(2974)
write.csv(mypval.HDAC6, file = "H:/CPTAC/CPTAC/mypval.HDAC6.csv")
mypval.HDAC7 <- subset.correl.gen(3693)
write.csv(mypval.HDAC7, file = "H:/CPTAC/CPTAC/mypval.HDAC7.csv")
mypval.HDAC8 <- subset.correl.gen(9589)
write.csv(mypval.HDAC8, file = "H:/CPTAC/CPTAC/mypval.HDAC8.csv")
mypval.HDAC9 <- subset.correl.gen(7064)
write.csv(mypval.HDAC9, file = "H:/CPTAC/CPTAC/mypval.HDAC9.csv")
mypval.HDAC10 <- subset.correl.gen(9523)
write.csv(mypval.HDAC10, file = "H:/CPTAC/CPTAC/mypval.HDAC10.csv")
mypval.HDAC11 <- subset.correl.gen(7604)
write.csv(mypval.HDAC11, file = "H:/CPTAC/CPTAC/mypval.HDAC11.csv")
mypval.NLRC5 <- subset.correl.gen(10369)
write.csv(mypval.NLRC5, file = "H:/CPTAC/CPTAC/mypval.NLRC5.csv")
mypval.SUPT7L <- subset.correl.gen(9429)
write.csv(mypval.SUPT7L, file = "H:/CPTAC/CPTAC/mypval.SUPT7L.csv")
mypval.WDR11 <- subset.correl.gen(783)
write.csv(mypval.WDR11, file = "H:/CPTAC/CPTAC/mypval.WDR11.csv")
mypval.FBXO11 <- subset.correl.gen(6836)
write.csv(mypval.FBXO11, file = "H:/CPTAC/CPTAC/mypval.FBXO11.csv")