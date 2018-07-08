rm(list=ls())
library(ggplot2)
library(ROCR)
library(flux)

args = commandArgs(T)

prefix = args[1]
outputdir = args[2]

summaryfile = paste(outputdir,"/prauc_val.txt",sep = "")
pltfile = paste(outputdir,"/",prefix,"pr_table.txt",sep = "")
pruac_file=paste(outputdir,"/",prefix,"prauc.txt",sep = "")
plot_pr_table<-read.table(pruac_file, sep='\t')
colnames(plot_pr_table)=c("rank","label")

pred<-prediction(0- plot_pr_table[,1], plot_pr_table[,2])
perf <- performance(pred,"prec","rec")
A=data.frame(perf@x.values[[1]], perf@y.values[[1]])
A[is.na(A)] <- 1
aucs<-auc(A[,1][1:sum(A[,1] <= 1)], A[,2][1:sum(A[,1] <= 1)])
aucs<-round(aucs,4)

cat(prefix, paste(aucs,"\n",sep=''), sep='\t',file=summaryfile,append=TRUE)
write.table(A, file=pltfile,sep='\t',quote=F,row.names=F,col.names=F)
