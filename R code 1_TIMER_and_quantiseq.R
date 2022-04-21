rm(list=ls())
# introduce the required package: immunedeconv
# citation for immnuedeconv: https://icbi-lab.github.io/immunedeconv/
require(immunedeconv)
# import the mRNA expression data of patients (n=207 following merging with RPPA and mRNA)
raw<-read.table("exp.txt", header=TRUE, row.names=1)
# excuete TIMER deconvolution, indicators show that the patients are breast cancer patients, arrays=TRUE show that the raw data are microarray data
res<-deconvolute(raw, "timer", indications=rep("BRCA", ncol(raw)), arrays=TRUE)
# excuete quantiseq deconvolution, this method works for both RNAseq data and microarray data, arrays=TRUE show that the raw data are microarray data
res2<-deconvolute(raw, "quantiseq", arrays=TRUE)
# export the results to csv files named result_TIMER and result_quantiseq
write.csv(res,"result_TIMER.csv")
write.csv(res2,"result_quantiseq.csv")
