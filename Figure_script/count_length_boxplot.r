library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(ggthemes)

rm(list=ls())
#definite genral ggplot template
boxcol<-function(n){
  colall<-c('#d7191c','#31a354','#0571b0')
  return(colall[c(1:n)])
}

mytemp<-theme(panel.background=element_rect(fill="white",colour=NA), panel.grid.major =element_blank(), panel.grid.minor = element_blank(),axis.line = element_line(colour = "black",size=0.5),axis.text.x = element_text(face = "bold"),axis.text.y =element_text(face = "bold"),
              axis.title.y=element_text(face="bold"),legend.text=element_text(face="bold"),legend.title=element_text(face="bold"))
widetemp<-theme(panel.background=element_rect(fill="white",colour=NA), panel.grid.major =element_blank(), panel.grid.minor = element_blank(),axis.line = element_line(colour = "black",size=0.5),axis.text.x = element_text(face = "bold",angle = 30,size = 12,hjust = 1),axis.text.y =element_text(face = "bold",size = 12),
                axis.title.y=element_text(face="bold",size=12),legend.text=element_text(face="bold",size=12),legend.title=element_text(face="bold"))

counts<-read.table("thesis_count_summary.txt",sep='\t')
counts<-data.frame(peaks = counts[,1],softs= counts[,2], mark=counts[,3],tissue=counts[,4],reps=counts[,5])
lengs<- read.table("thesis_length_summary.txt",sep='\t')
lengs<-data.frame(peaks = lengs[,1]+1, softs = lengs[,2], mark=lengs[,3],tissue = lengs[,4], reps= lengs[,5])




marks<-levels(counts$mark)
mar='DNase'
plot_table<-subset(counts, mark==mar&softs!="john"&reps!="Rep0 Peaks")
p<-ggplot(plot_table, aes(x=softs, y=peaks,fill=reps))+ylab("Number of peaks")+geom_boxplot(width=0.7,position = position_dodge(0.7),alpha=0.6,outlier.size = 0.1)+xlab(NULL)
#plot_title<-paste("Number of peaks in ",mar," peak calling algorithms",sep='')
q<-p+mytemp+scale_fill_manual("Sample",values = boxcol(3))+theme(legend.title=element_blank())+theme(legend.position='none')


plot_table2<-subset(lengs, mark==mar&softs!="john"&reps!="Rep0 Peaks")
p2<-ggplot(plot_table2, aes(x=softs, y=peaks, fill=reps))+ylab("log10(peak length)")+scale_y_log10()+xlab(NULL)+theme(legend.position='none')
q2<-p2+geom_boxplot(width=0.7,position = position_dodge(0.7),alpha=0.6,outlier.size = 0.5)+mytemp+scale_fill_manual("Sample",values = boxcol(3))+theme(legend.title=element_blank())

q3<-grid.arrange(q,q2,ncol=2)
#out_title<-paste(mar,"_counts.pdf",sep='')
ggsave(q, file="figure2a.pdf", height = 3, width=3.73)
ggsave(q2, file="figure2b.pdf", height = 3, width=3.5)

# plot for header
mar="H3K27ac"
plot_table<-subset(counts, mark==mar&softs!="john"&reps!="Rep0 Peaks")

p<-ggplot(plot_table, aes(x=softs, y=peaks,fill=reps))+ylab("Number of peaks")+geom_boxplot(width=0.7,position = position_dodge(0.7),alpha=0.6,outlier.size = 0.5)+xlab(NULL)
#plot_title<-paste("Number of peaks in ",mar," peak calling algorithms",sep='')
q<-p+widetemp+scale_fill_manual("Sample",values = boxcol(3))+theme(legend.title=element_blank())
q2<-q+theme(legend.position="top")
q2
ggsave(q2, file="box_header.pdf", height = 4, width=11)

for (mar in marks[-1]){
  plot_table<-subset(counts, mark==mar&softs!="john"&reps!="Rep0 Peaks")
  p<-ggplot(plot_table, aes(x=softs, y=peaks,fill=reps))+ylab("Number of peaks")+geom_boxplot(width=0.7,position = position_dodge(0.7),alpha=0.6,outlier.size = 0.5)+xlab(NULL)
  #plot_title<-paste("Number of peaks in ",mar," peak calling algorithms",sep='')
  q<-p+widetemp+scale_fill_manual("Sample",values = boxcol(3))+theme(legend.title=element_blank())
  
  plot_table2<-subset(lengs, mark==mar&softs!="john"&reps!="Rep0 Peaks")
  p2<-ggplot(plot_table2, aes(x=softs, y=peaks, fill=reps))+ylab("log10(peak length)")+scale_y_log10()+xlab(NULL)
  q2<-p2+geom_boxplot(width=0.7,position = position_dodge(0.7),alpha=0.6,outlier.size = 0.5)+widetemp+scale_fill_manual("Sample",values = boxcol(3))+theme(legend.title=element_blank())
  q3<-grid.arrange(q,q2,nrow=2)
  out_title<-paste("S1_",mar,".pdf",sep='')
  ggsave(q3, file=out_title, height = 6, width=11)
}  
  
# output count table for each mark
for (mar in marks){
  temptable<-subset(counts, mark==mar)
  newlabel<-unique(paste(temptable$softs,temptable$reps,sep=' '))
  tempmat<-matrix(NA, nrow = length(newlabel), ncol =length(unique(temptable$tissue)))
  tempmat<-data.frame(tempmat,row.names =newlabel)
  colnames(tempmat)<-unique(temptable$tissue)
  for (software in unique(temptable$softs)){
    for (tis in unique(temptable$tissue)){
      for (rep in unique(temptable$reps)){
        tempmat[paste(software,rep,sep=' '),tis]<-subset(temptable, softs==software&tissue==tis&reps==rep)$peaks
      }
    }
  }
  outname<-paste(mar,"_count.txt",sep='')
  write.table(tempmat,file=outname, sep='\t',quote=F,row.names = T, col.names = T)
}

# output length quantile table
myscale<-c("Min Length","1st quantile Length","Medium Length","3rd quantile Length","Max Length")
myscale2<-c(0,0.25,0.5,0.75,1)
newlabel<-c()
temprep<-c("Rep1 Peaks","Rep2 Peaks","Final Peaks")
for (mys in myscale){
  for (myrep in temprep){
    newlabel<-c(newlabel,paste(myrep, mys,sep=' '))
  }
}
for (mar in marks){
  temptable<-subset(lengs, mark==mar)
  
  tempmat<-matrix(NA, nrow=length(unique(temptable$softs)),ncol=length(newlabel))
  tempmat<-data.frame(tempmat, row.names = unique(temptable$softs))
  colnames(tempmat)<-newlabel
  for (soft in unique((temptable$softs))){
    for (rep in unique(temptable$reps)){
      tempdata<-subset(temptable,softs==soft&rep==reps)$peaks
      for (pos in c(1:5)){
        tempmat[soft, paste(rep, myscale[pos],sep=' ')]=floor(quantile(tempdata,myscale2[pos]))
      }
    }
  }
  outname<-paste(mar,"_length.txt",sep='')
  write.table(tempmat,file=outname, sep='\t',quote=F,row.names = T, col.names = T)
}

# replot H3K27ac without Hotspot and F-seq
plot_table<-subset(counts, mark=="H3K27ac"&softs!="john"&softs!="F-seq"&softs!="HotSpot"&reps!="Rep0 Peaks")
p<-ggplot(plot_table, aes(x=softs, y=peaks,fill=reps))+ylab("Number of peaks")+geom_boxplot(width=0.7,position = position_dodge(0.7),alpha=0.6,outlier.size = 0.5)+xlab(NULL)
#plot_title<-paste("Number of peaks in ",mar," peak calling algorithms",sep='')
q<-p+widetemp+scale_fill_manual("Sample",values = boxcol(3))+theme(legend.title=element_blank())
q<-p+widetemp+scale_fill_manual("Sample",values = boxcol(3))+theme(legend.title=element_blank())
ggsave(q, file ="figure3_new_d.pdf",height = 4, width = 10)

plot_table2<-subset(lengs, mark=="H3K27ac"&softs!="john"&softs!="F-seq"&softs!="HotSpot")
p2<-ggplot(plot_table2, aes(x=softs, y=peaks, fill=reps))+ylab("log10(peak length)")+scale_y_log10()+xlab(NULL)
q2<-p2+geom_boxplot(width=0.7,position = position_dodge(0.7),alpha=0.6,outlier.size = 0.5)+widetemp+scale_fill_manual("Sample",values = boxcol(3))+theme(legend.title=element_blank())
ggsave(q2, file ="figure3_new_e.pdf",height = 4, width = 10)



