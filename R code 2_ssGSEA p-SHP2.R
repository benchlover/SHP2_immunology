rm(list=ls()) # remove the previous data from environmment

# load the required libraries
library(GSVA)
library(GSEABase)

# import the subsets from upstream data cleaning. 
a = read.table("./mRNA_pSHP2_high.csv", stringsAsFactors = FALSE, header = TRUE, row.names = 1, sep = ",")
b = read.table("./mRNA_pSHP2_low.csv", stringsAsFactors = FALSE, header = TRUE, row.names = 1, sep = ",")

# save the subsets as matrix
a = as.matrix(a)
b = as.matrix(b)

# 1 load the Biocarta pathways
Biocarta <- getGmt("./c2.cp.biocarta.v7.4.symbols.gmt")
# Calculate the ssGSEA scores of the two subsets. 
ssgsea_score3 = gsva(a, Biocarta, method = "ssgsea", ssgsea.norm = TRUE, verbose = TRUE)   # signature 'matrix,list'
ssgsea_score4 = gsva(b, Biocarta, method = "ssgsea", ssgsea.norm = TRUE, verbose = TRUE)
# Output the results as CSV files for downstream analysis. 
write.csv(ssgsea_score3, "ssGSEAhigh2.csv")
write.csv(ssgsea_score4, "ssGSEAlow2.csv")

# 2 load the GO pathways
GO <- getGmt("./c5.go.v7.4.symbols.gmt")
# Calculate the ssGSEA scores of the two subsets. 
ssgsea_score7 = gsva(a, GO, method = "ssgsea", ssgsea.norm = TRUE, verbose = TRUE)   # signature 'matrix,list'
ssgsea_score8 = gsva(b, GO, method = "ssgsea", ssgsea.norm = TRUE, verbose = TRUE)
# Output the results as CSV files for downstream analysis. 
write.csv(ssgsea_score7, "ssGSEAhigh_GO.csv")
write.csv(ssgsea_score8, "ssGSEAlow_GO.csv")

# 3 load the KEGG pathways
KEGG <- getGmt("./c2.cp.kegg.v7.4.symbols.gmt")
# Calculate the ssGSEA scores of the two subsets. 
ssgsea_score9 = gsva(a, KEGG, method = "ssgsea", ssgsea.norm = TRUE, verbose = TRUE)   # signature 'matrix,list'
ssgsea_score10 = gsva(b, KEGG, method = "ssgsea", ssgsea.norm = TRUE, verbose = TRUE)
# Output the results as CSV files for downstream analysis. 
write.csv(ssgsea_score9, "ssGSEAhigh_KEGG.csv")
write.csv(ssgsea_score10, "ssGSEAlow_KEGG.csv")

# 4 load the wikipathways pathways
wiki <- getGmt("./c2.cp.wikipathways.v7.4.symbols.gmt")
# Calculate the ssGSEA scores of the two subsets. 
ssgsea_score11 = gsva(a, wiki, method = "ssgsea", ssgsea.norm = TRUE, verbose = TRUE)   # signature 'matrix,list'
ssgsea_score12 = gsva(b, wiki, method = "ssgsea", ssgsea.norm = TRUE, verbose = TRUE)
# Output the results as CSV files for downstream analysis. 
write.csv(ssgsea_score11, "ssGSEAhigh_wiki.csv")
write.csv(ssgsea_score12, "ssGSEAlow_wiki.csv")
