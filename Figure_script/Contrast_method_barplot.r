rm(list = ls())
library(ggplot2)
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


peak_table<-read.table("./data/prauc_table.txt",sep='\t')
prauc_table<-read.table("./data/full_pr_table.txt",sep='\t')
p1<-subset(peak_table, peak_table[,1]=="DFilter"|peak_table[,1]=="HotSpot2")
peaks<-subset(p1,p1[,2]=="H3K27ac"|p1[,2]=="DNase")
mypeak<-data.frame(mark=paste(peaks[,1],peaks[,2],sep="_"),tissue=peaks[,3],tissue2=peaks[,3],prauc=peaks[,4])
colnames(prauc_table)<-c("mark","tissue","tissue2","prauc")

#figure1 compare fcmethod and peak
fig1_tab<-subset(prauc_table,mark=="DFilter_DNase"|mark=="HotSpot2_DNase")
fig1_tab2<-merge(fig1_tab, mypeak,by=c("mark","tissue"))
prauc<-fig1_tab2[,4]-fig1_tab2[,6]
fig1_tab3<-cbind(fig1_tab2[,c(1:3)],prauc)
#write.table(fig1_tab3,file="./data/figure1.txt",sep='\t',quote = F,row.names = F,col.names = F)
#remove in one group
fig1_tab4<-read.table("./data/figure1.txt",sep='\t')
colnames(fig1_tab4)<-c("mark","tissue","tissue2","prauc")
fig1_tab5<-subset(fig1_tab4,mark=="DFilter_DNase")




fig1_new1<-read.table("./data/figure1_3.txt",sep='\t')
colnames(fig1_new1)<-c("mark","tissue","tissue2","prauc")

# figure2 zscore and fc
fig2_tab<-subset(prauc_table,mark=="DFilter_DNase_zscore"|mark=="HotSpot2_DNase_zscore"|mark=="DFilter_DNase_minus"|mark=="HotSpot2_DNase_minus")
#write.table(fig2_tab, file = "./data/figure2.txt",sep = "\t",quote = F,row.names = F,col.names = F)
for (tis in unique(fig1_tab$tissue)){
#tis="limb"
  fig2_tab2<-subset(fig2_tab, tissue==tis)
  fig2_tab3<-subset(fig1_tab,tissue==tis&tissue2!="average")
  fig2_tab4<-rbind(fig2_tab3,fig2_tab2)
  outname=paste("./figure2_",tis,".pdf",sep='')
  ptable<-as.data.frame(as.matrix(fig2_tab4))
  #ptable$mark<-factor(ptable$mark,levels=c("DFilter_DNase","HotSpot2_DNase","DFilter_DNase_zscore","HotSpot2_DNase_zscore", "DFilter_DNase_minus", "HotSpot2_DNase_minus"))
  levels(ptable$mark)<-c("DFilter DNase differential log signal","Hotspot2 DNase differential log signal","DFilter DNase differential Z-score","Hotspot2 DNase differential Z-score", "DFilter DNase differential signal", "Hotspot2 DNase differential signal")
  ptable$prauc<-as.data.frame(as.numeric(as.matrix(ptable$prauc)))
  p<-ggplot(ptable,aes(x=tissue2,y=prauc,group=mark,fill=mark))+geom_bar(stat = "identity",position = "dodge")
  q<-p+mytemp+scale_fill_manual(values = c('#fcbba1','#c6dbef','#fb6a4a','#6baed6','#a50f15','#08519c'))+xlab("Contrast Tissue")+ylab("PR-AUC")+theme(legend.title=element_blank())+ylim(c(0,0.55))
  #ggsave(q,filename = outname,height = 4,width = 5)
  ggsave(q,filename = "legend.pdf",height = 4,width = 5)
}

 








## weight or most distant
comp_table<-read.table("./data/weight_pr_summary.txt",sep='\t')
colnames(comp_table)<-c("mark","tissue1","tissue2","prauc")
plot_table<-subset(comp_table,mark=="H3K27ac_zscore")
p<-ggplot(plot_table,aes(x=tissue1,y=prauc,group=tissue2,fill=tissue2))+geom_bar(stat = "identity",position = "dodge")
q<-p+mytemp+xlab(NULL)+ylab("PR-AUC")+theme(legend.title=element_blank())+ylim(c(0,0.55))+scale_fill_manual(values = c("#e66101","#7b3294"))
ggsave(q,filename = "weight_compare.pdf",height = 4,width = 7)
