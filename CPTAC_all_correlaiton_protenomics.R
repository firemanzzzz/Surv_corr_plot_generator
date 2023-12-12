# Function to get column number based on gene name
get.colnumber <- function(gene.name){
  # Find the corresponding gene name in protein_gene_name_conver dataframe
  temp.name <- protein_gene_name_conver[which(protein_gene_name_conver$Gene.name == gene.name),3]
  # Count the occurrences of the gene name
  temp.times <- length(temp.name)
  # Loop through each occurrence of the gene name
  for (i in 1:temp.times) {
    # Find column number in proteome_X_corr where gene name matches
    temp.number <- grep(temp.name[i], colnames(proteome_X_corr))
    # Print the column number
    print(temp.number)
  }
}

# Function to calculate correlations for a specific protein
all.correl.gen <- function(protein.number){
  # Create an empty dataframe to store correlation results
  allpval.temp <- data.frame(matrix(ncol = 4, nrow = 0))
  colnames(allpval.temp) <- c("RefSeq.peptide.ID","pearson.r", "p.value", "n.rep")
  # Loop through all columns in proteome_X_corr dataframe
  for (i in 1:20742) {
    # Try to calculate correlation and handle any errors that might occur
    tryCatch({
      # Show progress by printing the iteration number
      message(i)
      # Store column name in the dataframe
      allpval.temp[i,1] <- colnames(proteome_X_corr)[i]
      # Calculate correlation between two protein columns
      temp <- list()
      temp <- cor.test(proteome_X_corr[,protein.number], proteome_X_corr[,i], 
                       na.action=na.omit, method = "pearson")
      # Store correlation coefficient, p-value, and number of observations
      allpval.temp[i,2] <- temp$estimate
      allpval.temp[i,3] <- temp$p.value
      allpval.temp[i,4] <- temp$parameter + 1
    }, error = function(e) {
      cat("ERROR :", conditionMessage(e), "\n")
    })
  }
  # Truncate the last two characters of RefSeq.peptide.ID column
  allpval.temp$RefSeq.peptide.ID <- substr(allpval.temp$RefSeq.peptide.ID, 1, 
                                           nchar(allpval.temp$RefSeq.peptide.ID) - 2)
  # Merge with protein_gene_name_conver dataframe based on RefSeq.peptide.ID
  allpval.temp <- merge(allpval.temp, protein_gene_name_conver, by = "RefSeq.peptide.ID", all.x = TRUE)
  # Return the resulting dataframe with correlation statistics
  return(allpval.temp)
}

# Calling get.colnumber function for specific gene names
get.colnumber("NLRC5")
get.colnumber("NOD2")
get.colnumber("FBXO11")
get.colnumber("WDR11")
get.colnumber("NFYA")
get.colnumber("NFYB")
get.colnumber("NFYC")
get.colnumber("HDAC1")
get.colnumber("HDAC2")
get.colnumber("HDAC3")
get.colnumber("HDAC4")
get.colnumber("HDAC5")
get.colnumber("HDAC6")
get.colnumber("HDAC7")
get.colnumber("HDAC8")
get.colnumber("HDAC9")
get.colnumber("HDAC10")
get.colnumber("HDAC11")
get.colnumber("SUPT7L")
get.colnumber("NACA")

# Calculating correlations for specific proteins and storing in variables
all_corre_NLRC5 <- all.correl.gen(10368)
all_corre_NOD2 <- all.correl.gen(12817)
all_corre_FBXO11 <- all.correl.gen(6835)
all_corre_WDR11 <- all.correl.gen(782)
all_corre_NFYA <- all.correl.gen(10805)
all_corre_NFYB <- all.correl.gen(9994)
all_corre_NFYC <- all.correl.gen(8061)
all_corre_HDAC1 <- all.correl.gen(3420)
all_corre_HDAC2 <- all.correl.gen(3421)
all_corre_HDAC3 <- all.correl.gen(6357)
all_corre_HDAC4 <- all.correl.gen(3052)
all_corre_HDAC5 <- all.correl.gen(3053)
all_corre_HDAC6 <- all.correl.gen(2973)
all_corre_HDAC7 <- all.correl.gen(3692)
all_corre_HDAC8 <- all.correl.gen(9588)
all_corre_HDAC9 <- all.correl.gen(11106)
all_corre_HDAC10 <- all.correl.gen(7063)
all_corre_HDAC11 <- all.correl.gen(9522)
all_corre_SUPT7L <- all.correl.gen(9428)
all_corre_NACA <- all.correl.gen(7603)


# File paths for exporting correlation results to CSV
a <- "H:/CPTAC/CPTAC/"
temp <- c("all_corre_NLRC5.csv", "all_corre_NOD2.csv", "all_corre_FBXO11.csv", "all_corre_WDR11.csv", 
          "all_corre_NFYA.csv", "all_corre_NFYB.csv", "all_corre_NFYC.csv", "all_corre_HDAC1.csv", "all_corre_HDAC2.csv", 
          "all_corre_HDAC3.csv", "all_corre_HDAC4.csv", "all_corre_HDAC5.csv", "all_corre_HDAC6.csv", "all_corre_HDAC7.csv", 
          "all_corre_HDAC8.csv", "all_corre_HDAC9.csv", "all_corre_HDAC10.csv", "all_corre_HDAC11.csv", "all_corre_SUPT7L.csv", "all_corre_NACA.csv")


# List of correlation dataframes
temp1 <- list(all_corre_NLRC5, all_corre_NOD2, all_corre_FBXO11, all_corre_WDR11, 
              all_corre_NFYA, all_corre_NFYB, all_corre_NFYC, all_corre_HDAC1, all_corre_HDAC2, 
              all_corre_HDAC3, all_corre_HDAC4, all_corre_HDAC5, all_corre_HDAC6, all_corre_HDAC7, 
              all_corre_HDAC8, all_corre_HDAC9, all_corre_HDAC10, all_corre_HDAC11, all_corre_SUPT7L, all_corre_NACA)


# Generating full file paths
for (i in 1:20) {
  temp0[i] <- paste(a,temp[i],sep="")
}

# Exporting each dataframe to CSV using respective file paths
for (i in 1:20) {
  write.csv(temp1[[i]], file = temp0[i])
}
