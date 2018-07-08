library(ggplot2)

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


softs<-c("DFilter", "Hotspot2")
tissue<-c("cranioface","cranioface","hindbrain","hindbrain","limb","limb","midbrain",
  "midbrain","neural_tube","neural_tube","forebrain","forebrain","heart","heart")
# raw top2000
for (soft in softs){
  Ptable2<-c()
  Ptable3<-c()
  xlabs<-paste("Rank by ",soft,sep='')
  ylabs<-("Rank by differential DNase signal")
  for (tis in unique(tissue)[1:5]){
    raw_name<-paste("./raw/",soft,"/",tis,".bed",sep='')
    raw_peak<-read.table(raw_name,sep='\t')
    raw_rank<-raw_peak[,4]
    weighted_name<-paste("./weighted/",soft,"/",tis,".bed",sep='')
    weighted_peak<-read.table(weighted_name,sep='\t')
    weighted_rank<-weighted_peak[,4]
    ptable<-data.frame(xval=raw_rank, yval=weighted_rank,Tissue=tis)
    ptable2<-subset(ptable, xval<=2000)
    ptable3<-subset(ptable, yval<=2000)
    print(soft)

    if (tis==tis){
      # raw peaks
      xlabs<-paste("Rank by ",soft,sep='')
      ylabs<-("Rank by differential DNase signal")
      p<-ggplot(ptable2,aes(x=xval,y=yval))+geom_point(size=0.2,col=rainbow(7)[1])+mytemp+xlab(xlabs)+ylab(ylabs)
      outname<-paste("./top2k_raw_rank/",soft,".",tis,".scatter.pdf",sep='')
      #ggsave(p,filename = outname,height = 4,width = 4.5)
      print(nrow(subset(ptable2,ptable2$yval>20000)))
      p<-ggplot(ptable2,aes(x=yval))+geom_density(size=1,alpha=0,col=rainbow(7)[1])+mytemp+xlab(xlabs)+ylab("Density")+geom_hline(yintercept = 0,color="white",size=1.2)+xlim(0,max(ptable2$yval*1.05))
      outname<-paste("./top2k_raw_rank/",soft,".",tis,".density.pdf",sep='')
      #ggsave(p,filename = outname,height = 4,width = 6)
      # rerank peaks
      p<-ggplot(ptable3,aes(x=yval,y=xval))+geom_point(size=0.2,col=rainbow(7)[1])+mytemp+ylab(xlabs)+xlab(ylabs)
      outname<-paste("./top2k_diff_rank/",soft,".",tis,".scatter.pdf",sep='')
      #ggsave(p,filename = outname,height = 4,width = 4.5)
     
      p<-ggplot(ptable3,aes(x=xval))+geom_density(size=1,alpha=0,col=rainbow(7)[1])+mytemp+xlab(xlabs)+ylab("Density")+geom_hline(yintercept = 0,color="white",size=1.2)+xlim(0,max(ptable3$xval*1.05))
      outname<-paste("./top2k_diff_rank/",soft,".",tis,".density.pdf",sep='')
      #ggsave(p,filename = outname,height = 4,width = 6)
    }
    Ptable2<-rbind(Ptable2, ptable2)
    Ptable3<-rbind(Ptable3, ptable3)
  }
  
  xlabs<-paste(xlabs)
  p<-ggplot(Ptable2,aes(x=yval, col=Tissue))+geom_density(size=1,alpha=0)+mytemp+xlab(ylabs)+ylab("Density")+geom_hline(yintercept = 0,color="white",size=1.2)+xlim(0,max(Ptable2$yval*1.05))+scale_color_manual("Tissue",values=rainbow(7)[1:5])
  outname<-paste("./top2k_raw_rank/",soft,"all.density.pdf",sep='')
  #ggsave(p,filename = outname,height = 4,width = 6)
  # rerank peaks
  p<-ggplot(Ptable3,aes(x=xval, col=Tissue))+geom_density(size=1,alpha=0)+mytemp+xlab(xlabs)+ylab("Density")+geom_hline(yintercept = 0,color="white",size=1.2)+xlim(0,max(Ptable3$xval*1.05))+scale_color_manual("Tissue",values=rainbow(7)[1:5])
  outname<-paste("./top2k_diff_rank/",soft,".all.density.pdf",sep='')
  #ggsave(p,filename = outname,height = 4,width = 6)
}



for (soft in "HOMER"){
  Ptable2<-c()
  Ptable3<-c()
  xlabs<-paste("Rank by ",soft,sep='')
  ylabs<-("Rank by differential H3K27ac signal")
  for (tis in unique(tissue)){
    raw_name<-paste("./raw/",soft,"/",tis,".bed",sep='')
    raw_peak<-read.table(raw_name,sep='\t')
    raw_rank<-raw_peak[,4]
    weighted_name<-paste("./weighted/",soft,"/",tis,".bed",sep='')
    weighted_peak<-read.table(weighted_name,sep='\t')
    weighted_rank<-weighted_peak[,4]
    ptable<-data.frame(xval=raw_rank, yval=weighted_rank,Tissue=tis)
    ptable2<-subset(ptable, xval<=2000)
    ptable3<-subset(ptable, yval<=2000)
    if (tis in tissue[1]){
      # raw peaks
      xlabs<-paste("Rank by ",soft,sep='')
      ylabs<-("Rank by differential DNase signal")
      p<-ggplot(ptable2,aes(x=xval,y=yval))+geom_point(size=0.2,col=rainbow(7)[1])+mytemp+xlab(xlabs)+ylab(ylabs)
      outname<-paste("./top2k_raw_rank/",soft,".",tis,".scatter.pdf",sep='')
      ggsave(p,filename = outname,height = 4,width = 4.5)
      
      p<-ggplot(ptable2,aes(x=yval))+geom_density(size=1,alpha=0,col=rainbow(7)[1])+mytemp+xlab(xlabs)+ylab("Density")+geom_hline(yintercept = 0,color="white",size=1.2)+xlim(0,max(ptable2$yval*1.05))
      outname<-paste("./top2k_raw_rank/",soft,".",tis,".density.pdf",sep='')
      ggsave(p,filename = outname,height = 4,width = 6)
      # rerank peaks
      p<-ggplot(ptable3,aes(x=yval,y=xval))+geom_point(size=0.2,col=rainbow(7)[1])+mytemp+ylab(xlabs)+xlab(ylabs)
      outname<-paste("./top2k_diff_rank/",soft,".",tis,".scatter.pdf",sep='')
      ggsave(p,filename = outname,height = 4,width = 4.5)
      
      p<-ggplot(ptable3,aes(x=xval))+geom_density(size=1,alpha=0,col=rainbow(7)[1])+mytemp+xlab(xlabs)+ylab("Density")+geom_hline(yintercept = 0,color="white",size=1.2)+xlim(0,max(ptable3$xval*1.05))
      outname<-paste("./top2k_diff_rank/",soft,".",tis,".density.pdf",sep='')
      ggsave(p,filename = outname,height = 4,width = 6)
    }
    Ptable2<-rbind(Ptable2, ptable2)
    Ptable3<-rbind(Ptable3, ptable3)
  }
  
  xlabs<-paste(xlabs)
  p<-ggplot(Ptable2,aes(x=yval, col=Tissue))+geom_density(size=1,alpha=0)+mytemp+xlab(ylabs)+ylab("Density")+geom_hline(yintercept = 0,color="white",size=1.2)+xlim(0,max(Ptable2$yval*1.05))+scale_color_manual("Tissue",values=rainbow(7))
  outname<-paste("./top2k_raw_rank/",soft,"all.density.pdf",sep='')
  #ggsave(p,filename = outname,height = 4,width = 6)
  # rerank peaks
  p<-ggplot(Ptable3,aes(x=xval, col=Tissue))+geom_density(size=1,alpha=0)+mytemp+xlab(xlabs)+ylab("Density")+geom_hline(yintercept = 0,color="white",size=1.2)+xlim(0,max(Ptable3$xval*1.05))+scale_color_manual("Tissue",values=rainbow(7))
  outname<-paste("./top2k_diff_rank/",soft,".all.density.pdf",sep='')
  #ggsave(p,filename = outname,height = 4,width = 6)
  #p<-ggplot(Ptable3,aes(x=xval,y=yval, col=Tissue))+geom_line()+mytemp+xlab(xlabs)+ylab("Density")+geom_hline(yintercept = 0,color="white",size=1.2)+xlim(0,max(Ptable3$xval*1.05))+scale_color_manual("Tissue",values=rainbow(7))
  #ggsave(p,filename = "label.pdf",height=4, width = 2)
}






