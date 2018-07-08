library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(ggthemes)
library(reshape2)
library(dplyr)
library(ROCR)
library(flux)
rm(list=ls())

setwd("/Users/apple/Desktop/project/dnase_project/thesis/figures/Latest_figure/Contrast_All_PR_curve")

mytemp<-theme(panel.background=element_rect(fill="white",colour=NA), panel.grid.major =element_blank(), panel.grid.minor = element_blank(),axis.line = element_line(colour = "black",size=0.5),axis.text.x = element_text(face = "bold"),axis.text.y =element_text(face = "bold"),
              axis.title.y=element_text(face="bold"),axis.title.x=element_text(face="bold"),legend.text=element_text(face="bold"),legend.title=element_text(face="bold"))

plot_tab<-read.table("./data/pr_all.txt",sep='\t')
pr_few<-read.table("pr_plot.txt",sep='\t')
tissue1=unique(pr_few[,5])[1:5]
for (tis in tissue1){
  pr_sum=c()
  pr_table<-subset(pr_few,pr_few[,4]=="DFilter DNase max signal difference"&pr_few[,5]==tis)
  pr_sum<-rbind(pr_sum,pr_table[pr_table[,1]!=1,c(1:4)])
  pr_table<-subset(pr_few,pr_few[,4]=="Hotspot2 DNase max signal difference"&pr_few[,5]==tis)
  pr_sum<-rbind(pr_sum,pr_table[pr_table[,1]!=1,c(1:4)])
  pr_table<-subset(pr_few,pr_few[,4]=="HOMER H3K27ac max Z-score difference"&pr_few[,5]==tis)
  pr_sum<-rbind(pr_sum,pr_table[pr_table[,1]!=1,c(1:4)])
  colnames(pr_sum)<-c("x_value","y_value","AUC","mark")
  temp<-subset(plot_tab,plot_tab[,3]=="DFilter"&plot_tab[,5]==tis)
  pred<-prediction(0-temp[,1],temp[,2])
  perf <- performance(pred,"prec","rec")
  A=data.frame(perf@x.values[[1]], perf@y.values[[1]])
  A[is.na(A)] <- 1
  aucs<-auc(A[,1][1:sum(A[,1] <= 1)], A[,2][1:sum(A[,1] <= 1)])
  aucs<-round(aucs,2)
  pr_table<-data.frame(x_value=A[,1],y_value=A[,2],AUC=aucs,mark="DFilter DNase weighted difference tissue, differential signal")
  pr_sum<-rbind(pr_sum,pr_table[-nrow(pr_table),])
  outtext<-paste(unique(pr_table$mark),tis,aucs,sep='\t')

  temp<-subset(plot_tab,plot_tab[,3]=="HOMER"&plot_tab[,5]==tis)
  pred<-prediction(0-temp[,1],temp[,2])
  perf <- performance(pred,"prec","rec")
  A=data.frame(perf@x.values[[1]], perf@y.values[[1]])
  A[is.na(A)] <- 1
  aucs<-auc(A[,1][1:sum(A[,1] <= 1)], A[,2][1:sum(A[,1] <= 1)])
  aucs<-round(aucs,2)
  pr_table<-data.frame(x_value=A[,1],y_value=A[,2],AUC=aucs,mark="HOMER H3K27ac weighted difference tissue, differential Z-score")
  pr_sum<-rbind(pr_sum,pr_table[-nrow(pr_table),])
  outtext<-paste(unique(pr_table$mark),tis,aucs,sep='\t')
  #write(x=outtext, file="./table2.txt",append = TRUE)
  
  temp<-subset(plot_tab,plot_tab[,3]=="Hotspot2"&plot_tab[,5]==tis)
  pred<-prediction(0-temp[,1],temp[,2])
  perf <- performance(pred,"prec","rec")
  A=data.frame(perf@x.values[[1]], perf@y.values[[1]])
  A[is.na(A)] <- 1
  aucs<-auc(A[,1][1:sum(A[,1] <= 1)], A[,2][1:sum(A[,1] <= 1)])
  aucs<-round(aucs,2)
  pr_table<-data.frame(x_value=A[,1],y_value=A[,2],AUC=aucs,mark="Hotspot2 DNase weighted difference tissue, differential signal")
  pr_sum<-rbind(pr_sum,pr_table[-nrow(pr_table),])
  outtext<-paste(unique(pr_table$mark),tis,aucs,sep='\t')
  #write(x=outtext, file="./table2.txt",append = TRUE)
  pr_end<-pr_table[nrow(pr_table),]
  
  pr_sum$mark<-factor(pr_sum$mark,levels=c("DFilter DNase weighted difference tissue, differential signal","Hotspot2 DNase weighted difference tissue, differential signal", "HOMER H3K27ac weighted difference tissue, differential Z-score","DFilter DNase max signal difference","Hotspot2 DNase max signal difference"  ,"HOMER H3K27ac max Z-score difference" ))
  p<-ggplot(pr_sum,aes(x=x_value,y=y_value,col=mark,group=mark,linetype=mark))+geom_line()+scale_linetype_manual(values = c("solid","solid","solid","dotted","dotted","dotted"))+geom_point(aes(x=1,y=pr_end[1,2]),col="black",size=0.6)
  #mytitle<-paste("PR curve of ",tis," DNase peaks")
  q<-p+xlab("recall")+ylab("precision")+mytemp+scale_color_manual(values =c('#d7191c','#045a8d','#31a354','#d7191c','#045a8d','#31a354'))+theme(legend.title=element_blank(),legend.position = "none")+xlim(c(0,1))+ylim(c(0,1))  
  #outfile=paste(tis,"_improve_prcurve.pdf",sep='')
  outfile=paste("1212_",tis,"_improve_prcurve.pdf",sep='')
  ggsave(q,file=outfile,width=3.5,height = 3)
}
 



tissue2<-c("forebrain","heart")

for (tis in tissue2){
  pr_sum=c()
  pr_table<-subset(pr_few,pr_few[,4]=="HOMER H3K27ac max Z-score difference"&pr_few[,5]==tis)
  print(pr_table[1,3])
  pr_sum<-rbind(pr_sum,pr_table[-nrow(pr_table),c(1:4)])
  colnames(pr_sum)<-c("x_value","y_value","AUC","mark")
  temp<-subset(plot_tab,plot_tab[,3]=="HOMER"&plot_tab[,5]==tis)
  pred<-prediction(0-temp[,1],temp[,2])
  perf <- performance(pred,"prec","rec")
  A=data.frame(perf@x.values[[1]], perf@y.values[[1]])
  A[is.na(A)] <- 1
  aucs<-auc(A[,1][1:sum(A[,1] <= 1)], A[,2][1:sum(A[,1] <= 1)])
  aucs<-round(aucs,2)
  print(aucs)
  pr_table<-data.frame(x_value=A[,1],y_value=A[,2],AUC=aucs,mark="HOMER H3K27ac weighted Z-score difference")
  pr_sum<-rbind(pr_sum,pr_table[-nrow(pr_table),])
  outtext<-paste(unique(pr_table$mark),tis,aucs,sep='\t')
  #write(x=outtext, file="./table2.txt",append = TRUE)
  pr_sum$mark<-factor(pr_sum$mark,levels=c("HOMER H3K27ac weighted Z-score difference","HOMER H3K27ac max Z-score difference"))
  p<-ggplot(pr_sum,aes(x=x_value,y=y_value,group=mark,linetype=mark))+geom_line(col='#31a354')+scale_linetype_manual(values = c("solid","dotted"))

  q<-p+xlab("recall")+ylab("precision")+mytemp+theme(legend.title=element_blank())+theme(legend.position="none")+geom_point(aes(x=1,y=pr_end[1,2]),col="black",size=0.6)
  outfile=paste("1212_",tis,"_improve_prcurve.pdf",sep='')
  ggsave(q,file=outfile,width=3.5,height = 3)
}

