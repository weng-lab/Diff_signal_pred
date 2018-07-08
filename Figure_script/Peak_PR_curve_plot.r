library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(ggthemes)
library(reshape2)
library(dplyr)
library(ROCR)
library(flux)
rm(list=ls())
boxcol<-function(n){
  colall<-c('#d7191c','#31a354','#0571b0')
  return(colall[c(1:n)])
}
linecol<-function(n){
  colall<-c('#d7191c','#31a354','#756bb1','#0571b0','#d95f0e','#bdbdbd')
  return(colall[c(1:n)])
}

mytemp<-theme(panel.background=element_rect(fill="white",colour=NA), panel.grid.major =element_blank(), panel.grid.minor = element_blank(),axis.line = element_line(colour = "black",size=0.5),axis.text.x = element_text(face = "bold"),axis.text.y =element_text(face = "bold"),
              axis.title.y=element_text(face="bold"),axis.title.x=element_text(face="bold"),legend.text=element_text(face="bold"),legend.title=element_text(face="bold"))

###  H3K27ac PR curve

dnase_tab<-read.table("pr_all.txt",sep='\t')

colnames(dnase_tab)<-c("rank","label","soft","tissue")

#softs<-c("DFilter", "HOMER", "MACS2", "MUSIC.broad", "MUSIC.punctuate", "F-seq", "RSEG", "FindER", "BCP", "MOSAiCS")
# only use top 4 software
softs<-c("DFilter", "HOMER", "MACS2","MUSIC punctate")
tissues<-levels(dnase_tab$tissue)
pr_sum<-c()
top200_list<-c()
pr_end<-c()
for (tis in tissues){
  for (software in softs){
    temp<-subset(dnase_tab,tissue==tis&soft==software)
    top200_temp<-subset(temp,temp[,1]<=200)
    top200_tab<-data.frame(x_value=sum(top200_temp[,2])/sum(temp[,2]), y_value=sum(top200_temp[,2])/ nrow(top200_temp),soft=software,mark="H3K27ac",tissue=tis)
    top200_list<-rbind(top200_list,top200_tab)
    if (nrow(temp)>0){
      pred<-prediction(0-temp[,1],temp[,2])
      perf <- performance(pred,"prec","rec")
      A=data.frame(perf@x.values[[1]], perf@y.values[[1]])
      A[is.na(A)] <- 1
      aucs<-auc(A[,1][1:sum(A[,1] <= 1)], A[,2][1:sum(A[,1] <= 1)])
      aucs<-round(aucs,2)
      #print(aucs)
      pr_table<-data.frame(x_value=A[,1],y_value=A[,2],AUC=aucs,soft=software,mark="H3K27ac",tissue=tis)
      #auc_tmp<-c(software,"DNase",tis,aucs)
      #prauc_table<-rbind(prauc_table,auc_tmp)
      pr_sum<-rbind(pr_sum,pr_table[-nrow(pr_table),])
      
    }
  }
  pr_end<-rbind(pr_end,pr_table[nrow(pr_table),])
}
labels<-paste(pr_sum$soft,pr_sum$AUC,sep=":")
pr_sum<-cbind(pr_sum,labels)


for (tis in tissues){
  i=1
  plot_table<-subset(pr_sum, tissue==tis)
  top200_p<-subset(top200_list,tissue==tis)
  p<-ggplot(plot_table,aes(x=x_value,y=y_value,col=soft,group=soft))+geom_line()+geom_point(data=top200_p,aes(x=x_value,y=y_value,col=soft),size=3,shape=1)+scale_shape(solid = FALSE)
  #mytitle<-paste("PR curve of ",tis," DNase peaks")
  q<-p+xlab("recall")+ylab("precision")+mytemp+scale_color_manual("Algorithm",values =c(linecol(3),"#e6550d")) +theme(legend.title=element_blank())+geom_point(aes(x=1,y=pr_end[i,2]),col="black",size=0.6)
  outfile=paste("1212_H3K27ac_",tis,"_prcurve.pdf",sep='')
  
  ggsave(q,file=outfile,width=5,height = 3)
  i=i+1
}


#label
for (tis in tissues[c(2,3)]){
  outname<-paste(tis,"_label.pdf",sep='')
  plot_table<-subset(pr_sum, tissue==tis)
  top200_p<-subset(top200_list,tissue==tis)
  p<-ggplot(plot_table,aes(x=x_value,y=y_value,col=labels,group=labels))+geom_line()#+geom_point(data=top200_p,aes(x=x_value,y=y_value,col=soft),size=3,shape=1)+scale_shape(solid = FALSE)
  #mytitle<-paste("PR curve of ",tis," DNase peaks")
  q<-p+xlab("recall")+ylab("precision")+mytemp+scale_color_manual("Algorithm",values =c(linecol(3),"#e6550d")) +theme(legend.title=element_blank())
  outfile=paste("H3K27ac_",tis,"_label.pdf",sep='')
  ggsave(q,file=outfile,width=5.5,height = 3)
}





