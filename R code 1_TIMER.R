rm(list=ls())
# introduce the required package: immunedeconv
# citation for immnuedeconv: https://icbi-lab.github.io/immunedeconv/
require(immunedeconv)
# import the mRNA expression data of patients (n=207 following merging with RPPA and mRNA)
raw<-read.table("exp.txt", header=TRUE, row.names=1)
# excuete TIMER deconvolution, indicators show that the patients are breast cancer patients
res<-deconvolute(raw, "timer", indications=rep("BRCA", ncol(raw)))
# export the results to csv files named result_TIMER
write.csv(res,"result_TIMER.csv")
