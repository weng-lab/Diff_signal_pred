
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(ggthemes)

rm(list=ls())
boxcol<-function(n){
  colall<-c('#1a9850','#fc8d59')
  return(colall[c(1:n)])
}
mytemp<-theme(panel.background=element_rect(fill="white",colour=NA), panel.grid.major =element_blank(), panel.grid.minor = element_blank(),axis.line = element_line(colour = "black",size=0.5),axis.text.x = element_text(colour = "black"),axis.text.y =element_text(colour = "black"),
              axis.title.y=element_text(colour = "black"))
widetemp<-theme(panel.background=element_rect(fill="white",colour=NA), panel.grid.major =element_blank(), panel.grid.minor = element_blank(),axis.line = element_line(colour = "black",size=0.5),axis.text.x = element_text(angle = 30,size = 12,hjust = 1,colour = "black"),axis.text.y =element_text(size = 12,colour = "black"),
                axis.title.y=element_text(size=12,colour = "black"),legend.text=element_text(size=12))


rep_cons<-read.table("replicate_consistency.txt",sep='\t')
size_cons<-read.table('consistency_summary.txt',sep='\t')

rep_cons<-subset(rep_cons,rep_cons[,2]!="MACS2.broad"&rep_cons[,2]!="MACS2.gapped")
rep_cons<-subset(rep_cons,rep_cons[,2]!="F-seq"|rep_cons[,3]!='DNase')
size_cons<-subset(size_cons,size_cons[,3]!="MACS2.broad"&size_cons[,3]!="MACS2.gapped")
size_cons<-subset(size_cons,size_cons[,3]!="F-seq"|size_cons[,2]!='DNase')
rep_table<-data.frame(value=as.numeric(rep_cons[,1]),soft=rep_cons[,2],mark=rep_cons[,3],types=rep_cons[,4])
marks<-levels(rep_table$mark)
softs<-levels(rep_table$soft)
size_table<-data.frame(value=as.numeric(size_cons[,1]),soft=size_cons[,3],mark=size_cons[,2],types=size_cons[,4])
size_table<-subset(size_table, types !="Rep0 vs Rep0.pr1,pr2")
size_table$types=factor(size_table$types, levels<-c("Rep0 vs Rep1,Rep2","Rep0 vs Rep.pr1,pr2"))
rep_table$types=factor(rep_table$types,levels = c("Top 1000","Top 2000","Top 10000","Top 20000","Top 40000","All Peaks"))


mar="DNase"
plot_table<-subset(rep_table, mark==mar)
p<-ggplot(plot_table,aes(x=soft, y=value,fill=types))+xlab(NULL)+ylim(c(0,1))+geom_boxplot(width=0.7,position = position_dodge(0.7),alpha=0.6,outlier.size = 0.5)+ylab(NULL)#+ylab("Consistency between 2 bio-replicates")
q<-p+mytemp+theme(legend.title=element_blank())+scale_fill_manual(values = rainbow(6))+theme(legend.position = "none")

plot_table2<-subset(size_table, mark==mar)
p2<-ggplot(plot_table2,aes(x=soft, y=value,fill=types))+xlab(NULL)+ylab(NULL)+ylim(c(0,1))+geom_boxplot(width=0.6,position = position_dodge(0.7),alpha=1,outlier.size = 0.5)#+ylab("Consistency between Rep0\n and Reps and pseudo-Reps")
q2<-p2+mytemp+scale_fill_manual(values = boxcol(2))+theme(legend.title=element_blank())+theme(legend.position = "none")


ggsave(q, file="figure1a.pdf", height = 3, width=4,useDingbats=FALSE )
ggsave(q2, file="figure1b.pdf", height = 3, width=4,useDingbats=FALSE)

for (mar in marks[-1]){
  plot_table<-subset(rep_table, mark==mar)
  p<-ggplot(plot_table,aes(x=soft, y=value,fill=types))+xlab(NULL)+ylab("Consistency between 2 bio-replicates")+ylim(c(0,1))+geom_boxplot(width=0.7,position = position_dodge(0.7),alpha=0.6,outlier.size = 0.5)
  q<-p+widetemp+scale_fill_manual(values = c(rainbow(6)[1:4],rainbow(6)[6]))+theme(legend.title=element_blank())
  out_title<-paste("S3_",mar,"_cons_1.pdf",sep='')
  ggsave(q, file=out_title, height = 4, width=11,useDingbats=FALSE)
  plot_table2<-subset(size_table, mark==mar)
  p2<-ggplot(plot_table2,aes(x=soft, y=value,fill=types))+xlab(NULL)+ylab("Consistency between Rep0\n and Reps and pseudo-Reps")+ylim(c(0,1))+geom_boxplot(width=0.6,position = position_dodge(0.7),alpha=1,outlier.size = 0.5)
  q2<-p2+widetemp+scale_fill_manual(values = boxcol(2))+theme(legend.title=element_blank())
  
  out_title<-paste("S3_",mar,"_cons_2.pdf",sep='')
  ggsave(q2, file=out_title, height = 4, width=12,useDingbats=FALSE)
}

mar="H3K27ac"

plot_table<-subset(rep_table, mark==mar)
p<-ggplot(plot_table,aes(x=soft, y=value,fill=types))+xlab(NULL)+ylim(c(0,1))+geom_boxplot(width=0.7,position = position_dodge(0.7),alpha=0.6,outlier.size = 0.5)+ylab(NULL)#+ylab("Consistency between 2 bio-replicates")
q<-p+widetemp+scale_fill_manual(values = c(rainbow(6)[1:4],rainbow(6)[6]))+theme(legend.title=element_blank())+theme(legend.position = "none")
out_title<-paste("S3_",mar,"_cons_1.pdf",sep='')
ggsave(q, file=out_title, height = 4, width=9)
plot_table2<-subset(size_table, mark==mar)
p2<-ggplot(plot_table2,aes(x=soft, y=value,fill=types))+xlab(NULL)+ylim(c(0,1))+geom_boxplot(width=0.7,position = position_dodge(0.7),alpha=1,outlier.size = 0.5)+ylab(NULL)#+ylab("Consistency between Rep0\n and Reps and pseudo-Reps")
q2<-p2+widetemp+scale_fill_manual(values = boxcol(3))+theme(legend.title=element_blank())+theme(legend.position = "none")
out_title<-paste("S3_",mar,"_cons_2.pdf",sep='')
ggsave(q2, file=out_title, height = 4, width=9)

