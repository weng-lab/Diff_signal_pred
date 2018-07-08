library(ggplot2)

boxcol<-function(n){
  colall<-c('#1a9850','#984ea3','#fc8d59')
  return(colall[c(1:n)])
}
mytemp<-theme(panel.background=element_rect(fill="white",colour=NA), panel.grid.major =element_blank(), panel.grid.minor = element_blank(),axis.line = element_line(colour = "black",size=0.6),axis.text.x = element_text(face = "bold"),axis.text.y =element_text(face = "bold"),
              axis.title.y=element_text(face="bold"),axis.title.x=element_text(face="bold"),legend.text=element_text(face="bold"),legend.title=element_text(face="bold"))

mytab<-read.table("./summary.txt",sep='\t')
colnames(mytab)<-c("Motif","Seq","Pval","Count","Soft","Tissue")

for (tis in c("cranioface","limb","neural_tube","hindbrain","midbrain","forebrain","heart")){
  for (soft in c("HOMER","Hotspot2","DFilter")){
    sub1<-subset(mytab, Soft==soft & Tissue==tis)
    if (nrow(sub1)>0){
      sub2<-subset(mytab, Soft==paste(soft,"_weighted",sep='') & Tissue==tis)
      subp<-merge(sub1,sub2,by=c("Tissue","Motif","Seq"))
      subplots<-subp[,c(1,4,5,7,8)]
      
      colnames(subplots)<-c("Tissue","RawP","RawC","WeightP","WeightC")
      subplots$RawP[subplots$RawP>100] <- 100
      subplots$WeightP[subplots$WeightP>100] <- 100
      subplots2<-subset(subplots,subplots$WeightP>2 | subplots$RawP>2)
      p<-ggplot(subplots2)+geom_point(aes(x=RawP,y=WeightP))+mytemp+xlab(paste(soft," origin peak -log10 pvalue",sep=''))+ylab(paste(soft," differential signal peak -log10 pvalue",sep=''))+geom_abline()+xlim(c(0,100))+ylim(c(0,100))
      mytext <- paste("italic(p)"," == ", format(t.test(subplots2$WeightP,subplots2$RawP,paried=T)$p.value, scientific=T, digits = 3),sep='')
      q<-p+annotate(geom="text",label=mytext,x = 50,y=30,fontface="bold",size=4.5,parse=T)
      print(paste(soft,tis))
      print(t.test(subplots2$WeightP,subplots2$RawP,paried=T)$p.value)
      #ggsave(q,filename = paste("./figures/pval/",soft,"_",tis,"_Pval.pdf",sep=''), useDingbats=FALSE, height = 5,width=6)
      # p<-ggplot(subplots2)+geom_point(aes(x=RawC,y=WeightC))+mytemp+xlab(paste("Number of ",soft," origin peak TF bings",sep=''))+ylab(paste("Number of ",soft," differential signal peak TF bindings",sep=''))+geom_abline()
      #ggsave(p,filename = paste("./figures/count/",soft,"_",tis,"_Count.pdf",sep = ""))

    }
  }
}
