library(ggplot2)
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

mytemp<-theme(panel.background=element_rect(fill="white",colour="black"), panel.grid.major =element_blank(), panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),axis.text.x = element_text(face = "bold"),axis.text.y =element_text(face = "bold"),
              axis.title.y=element_text(face="bold"),axis.title.x=element_text(face="bold"),legend.text=element_text(face="bold"),legend.title=element_text(face="bold"))


args = commandArgs(T)
pltdir = args[1]
outdir = args[2]

allf = dir(pltdir)
alltab = c()

for (myf in allf){
    intmp =read.table(paste(pltdir,"/",myf,sep=''), sep='\t')
    A = intmp
    A[is.na(A)] <- 1
    aucs<-auc(A[,1][1:sum(A[,1] <= 1)], A[,2][1:sum(A[,1] <= 1)])
    label<-paste(myf,round(aucs,2),sep=":")

    intab = cbind(intmp, myf, label)
    alltab<-rbind(alltab, intab)
}
colnames(alltab)<-c("x_value","y_value","sample","label")

p<-ggplot(alltab,aes(x=x_value,y=y_value,col=label,group=label))+geom_line()

q<-p+xlab("recall")+ylab("precision")+mytemp+scale_color_manual(values =c(linecol(2)))+theme(legend.title=element_blank())
outfile=paste(outdir,"/contrast_prcurve.pdf",sep='')
ggsave(q,file=outfile,width=5.5,height = 3)

