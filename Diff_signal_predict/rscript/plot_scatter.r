rm(list=ls())
library(ggplot2)

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


# 输入dist文件和PRAUC文件

args = commandArgs(T)
inputdir = args[1]
outputfile = args[2]

dist_file = paste(inputdir,"/distance.txt",sep='')
prauc_file = paste(inputdir,"/prauc_val.txt",sep='')

dists <- read.table(dist_file,sep='\t')
praucs <-read.table(prauc_file,sep='\t')
praucs[,2] <- praucs[,2]-as.numeric(praucs[1,2])
my_tab<-data.frame(cbind(praucs,dists))
colnames(my_tab)<-c("samples","prauc","distance")
p<-ggplot(my_tab,aes(x=distance, y=prauc))+geom_point()+mytemp+xlab("distance")+ylab("prauc diff")
cor_square <-round(cor(my_tab$prauc,my_tab$distance)^2, 2)
r_square_text<-paste('R^{2}==',cor_square,sep='')
q<- p+annotate(geom='text',x=-Inf,y=Inf,vjust=1.3,hjust=-0.1,label=r_square_text,colour='black',size=5,fontface='bold',parse=TRUE)+theme(legend.position='none')
ggsave(q,filename=outputfile,height=5,width=6)