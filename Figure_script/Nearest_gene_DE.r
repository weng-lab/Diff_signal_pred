rm(list=ls())
library(ggplot2)
library(DESeq2)
library(plyr)
library("TxDb.Mmusculus.UCSC.mm10.ensGene")


txdb<-TxDb.Mmusculus.UCSC.mm10.ensGene
mm10_gene<-genes(txdb)
setwd("./data/mouse_fpkm/")
file_list=dir("./")
tissue<-c("cranioface","cranioface","forebrain","forebrain","heart","heart",
          "hindbrain","hindbrain","limb","limb","liver","liver","midbrain",
          "midbrain","neural_tube","neural_tube")

temp<-read.table(file_list[1],sep='\t')
colnames(temp)<-c("gene",strsplit(file_list[1],".tsv")[[1]])
expr_table<-temp
for (myfile in file_list[-1]){
  temp<-read.table(myfile,sep='\t')
  colnames(temp)<-c("gene",strsplit(myfile,".tsv")[[1]])
  expr_table<-merge(expr_table,temp,by="gene")
}

expr_table<-data.frame(expr_table)
expr_gene<-subset(expr_table,apply(expr_table[,-1],1,max)>1)

a=data.frame(mm10_gene)

high<-as.vector(expr_gene[,1])
high2<-c()
for (ii in high){
  j<-strsplit(ii,"\\.")[[1]][1]
  high2<-c(high2,j)
}
out_gene<-a[a$gene_id %in% high2, ]

write.table(out_gene,file="../expressed_gene.txt", row.names=F,col.names = F,quote = F,sep='\t')



## DE analysis
rm(list=ls())
setwd("/Users/apple/Desktop/project/dnase_project/thesis/figures/Latest_figure/Nearest_gene_analysis/data/mouse_expr/")
tissue<-c("cranioface","cranioface","forebrain","forebrain","heart","heart",
          "hindbrain","hindbrain","limb","limb","liver","liver","midbrain",
          "midbrain","neural_tube","neural_tube")

file_list=dir("./")
tissue<-c("cranioface","cranioface","forebrain","forebrain","heart","heart",
          "hindbrain","hindbrain","limb","limb","liver","liver","midbrain",
          "midbrain","neural_tube","neural_tube")

temp<-read.table(file_list[1],sep='\t')
colnames(temp)<-c("gene",strsplit(file_list[1],".txt")[[1]])
expr_table<-temp
for (myfile in file_list[-1]){
  temp<-read.table(myfile,sep='\t')
  colnames(temp)<-c("gene",strsplit(myfile,".txt")[[1]])
  expr_table<-merge(expr_table,temp,by="gene")
}
gene_table<-data.frame(expr_table[,-1])
rownames(gene_table)<-expr_table[,1]
condition = rep("control",10)
brain_list=c(3,4,7,8,13,14,15,16)
for (i in seq(1,16,2)){
  if (i %in% brain_list){
    total_list=c(i,i+1,c(1:16)[-brain_list])
  }else{
    total_list=c(i,i+1,c(1:16)[brain_list])
  }
  gene_table2<-gene_table[,total_list]
  condition2=condition
  condition2[c(1,2)]=c("treat","treat")
  coldata<- data.frame(row.names = colnames(gene_table2), condition2)
  dds <- DESeqDataSetFromMatrix(countData = gene_table2, colData = coldata, design = ~condition2)
  dds <- DESeq(dds)
  res <- results(dds,contrast = c("condition2","treat", "control"))
  res <- res[order(res$pvalue),]
  outdata<-data.frame(res)
  outrow<-rownames(outdata)
  high2<-c()
  for (ii in outrow){
    j<-strsplit(ii,"\\.")[[1]][1]
    high2<-c(high2,j)
  }
  outdata2<-data.frame(gene=high2,fc=res$log2FoldChange)
  outdata2<-subset(outdata2,fc!="NA")
  tis=tissue[i]
  outname<-paste("../",tis,"_fc.txt",sep='')
  write.table(outdata2,file = outname,sep='\t',quote=F,row.names = F,col.names = F)
}




### get nearest gene expression distribution
rm(list = ls())
setwd("../closest_gene/")

ptype<-c("DFilter_DNase_peak","DFilter_weighted_signal","Hotspot2_DNase_peak","Hotspot2_weighted_signal"
         ,"HOMER_H3K27ac_peak","HOMER_weighted_zscore")
tissue<-c("cranioface","cranioface","hindbrain","hindbrain","limb","limb","midbrain",
          "midbrain","neural_tube","neural_tube","forebrain","forebrain","heart","heart")
out_fc<-c()
for (tis in unique(tissue)[1:5]){
  inname<-paste("../",tis,"_fc.txt",sep='')
  fc_tab<-read.table(inname,sep='\t')
  colnames(fc_tab)<-c("gene","fc")
  for (peaks in ptype){
    inname2<-paste(peaks,"_",tis,"_gene.txt",sep='')
    gen_tab<-read.table(inname2,sep='\t')
    colnames(gen_tab)<-"gene"
    temp<-merge(gen_tab,fc_tab,by="gene")
    temp2<-data.frame(fc=temp[,2],method=peaks, tissue=tis)
    out_fc<-rbind(out_fc,temp2)
  }
}
for (tis in unique(tissue)[6:7]){
  inname<-paste("../",tis,"_fc.txt",sep='')
  fc_tab<-read.table(inname,sep='\t')
  colnames(fc_tab)<-c("gene","fc")
  for (peaks in ptype[5:6]){
    inname2<-paste(peaks,"_",tis,"_gene.txt",sep='')
    gen_tab<-read.table(inname2,sep='\t')
    colnames(gen_tab)<-"gene"
    temp<-merge(gen_tab,fc_tab,by="gene")
    temp2<-data.frame(fc=temp[,2],method=peaks, tissue=tis)
    out_fc<-rbind(out_fc,temp2)
  }
}

write.table(out_fc,file="../../figures/fc_all_table.txt",sep='\t',quote=F,row.names = F,col.names = F)

# plot and t test analysis
Tissue=c("face","face","hindbrain","hindbrain","limb","limb","midbrain",
         "midbrain","neural tube","neural tube","forebrain","forebrain","heart","heart")
setwd("/Users/apple/Desktop/project/dnase_project/thesis/figures/Latest_figure/Nearest_gene_analysis")
outfc<-read.table("./figures/fc_all_table.txt",sep='\t')
colnames(outfc)<-c("fc","method","tissue")
mytemp<-theme(panel.background=element_rect(fill="white",colour=NA), panel.grid.major =element_blank(), panel.grid.minor = element_blank(),axis.line = element_line(colour = "black",size=0.5),axis.text.x = element_text(face = "bold"),axis.text.y =element_text(face = "bold"),
              axis.title.y=element_text(face="bold"),axis.title.x=element_text(face="bold"),legend.text=element_text(face="bold"),legend.title=element_text(face="bold"),panel.border=element_blank())

out_fc2<-subset(outfc,method=="DFilter DNase peak"|method=="DFilter weighted signal")
levels(out_fc2$method)<-c("DFilter DNase peaks, ranked by DFilter","DFilter DNase peaks, ranked by differential signal against weighted tissues","","","","")
meds<-ddply(out_fc2,.(tissue,method),summarise,med=median(fc))
p<-ggplot(out_fc2,aes(x=tissue,y=fc,col=method))+xlab(NULL)+ylab("Nearest gene Log2 fold changes")+ylim(c(-5,5))+geom_boxplot(width=0.7,position = position_dodge(0.8),alpha=1,outlier.size = 0.3)
q<-p+mytemp+theme(legend.title=element_blank())+geom_text(data=meds,aes(x=tissue,y=med,label=round(med,2),group=method),size=3,hjust=1.15,vjust=-3,position = position_dodge(0.8))+theme(legend.position = "none")
#ggsave(q,filename = "./figures/DFilter_fc_diff_label.pdf",height=3.5, width = 6)

ggsave(q,filename = "./figures/DFilter_fc_diff.pdf",height=3.5, width = 5)

for (tis in unique(Tissue)[1:5]){
  a=subset(out_fc2,method=="DFilter DNase peak"&tissue==tis)$fc
  b=subset(out_fc2,method=="DFilter DNase weighted distance tissue, differential signal"&tissue==tis)$fc
  outtext<-paste("DFilter",tis,t.test(x=b,y=a)$p.value,t.test(x=b,y=a)$stat)
  #write.table(outtext,file="./figures/pvalues.txt",sep = '\t',quote=F,row.names = F,col.names = F,append = T)
}



out_fc2<-subset(outfc,method=="Hotspot2 DNase peak"|method=="Hotspot2 weighted signal")
levels(out_fc2$method)<-c("","","","","Hotspot2 DNase peak, ranked by Hotspot2","Hotspot2 DNase peaks, ranked by differential signal against weighted tissues")
meds<-ddply(out_fc2,.(tissue,method),summarise,med=median(fc))
p<-ggplot(out_fc2,aes(x=tissue,y=fc,col=method))+xlab(NULL)+ylab("Nearest gene Log2 fold changes")+ylim(c(-5,5))+geom_boxplot(width=0.7,position = position_dodge(0.8),alpha=1,outlier.size = 0.3)
q<-p+mytemp+theme(legend.title=element_blank())+geom_text(data=meds,aes(x=tissue,y=med,label=round(med,2),group=method),size=3,hjust=1.15,vjust=-3,position = position_dodge(0.8))+theme(legend.position = "none")
ggsave(q,filename = "./figures/Hotspot2_fc_diff.pdf",height=3.5, width = 5)
#ggsave(q,filename = "./figures/Hotspot2_fc_diff_label.pdf",height=3.5, width = 6.5)
for (tis in unique(Tissue)[1:5]){
  a=subset(out_fc2,method=="Hotspot2 DNase peak"&tissue==tis)$fc
  print(length(a[a>1])/length(a))
  b=subset(out_fc2,method=="Hotspot2 DNase weighted distance tissue, differential signal"&tissue==tis)$fc
  #print(length(b[b>1])/length(b))
  outtext<-paste("HotSpot2",tis,t.test(x=b,y=a)$p.value,t.test(x=b,y=a)$stat)
  #write.table(outtext,file="./figures/pvalues.txt",sep = '\t',quote=F,row.names = F,col.names = F,append = T)
}

out_fc2<-subset(outfc,method=="HOMER H3K27ac peak"|method=="HOMER weighted zscore")
levels(out_fc2$method)<-c("","","HOMER H3K27ac peak, ranked by HOMER","HOMER H3K27ac peaks, ranked by differential Z-scores against weighted tissues","","")
meds<-ddply(out_fc2,.(tissue,method),summarise,med=median(fc))
p<-ggplot(out_fc2,aes(x=tissue,y=fc,col=method))+xlab(NULL)+ylab("Nearest gene Log2 fold changes")+ylim(c(-5,5))+geom_boxplot(width=0.7,position = position_dodge(0.8),alpha=1,outlier.size = 0.3)
q<-p+mytemp+theme(legend.title=element_blank())+geom_text(data=meds,aes(x=tissue,y=med,label=round(med,2),group=method),size=3,hjust=1.15,vjust=-4,position = position_dodge(0.8))+theme(legend.position = "none")
ggsave(q,filename = "./figures/HOMER_fc_diff.pdf",height=3.5, width = 7)
#ggsave(q,filename = "./figures/HOMER_fc_diff_label.pdf",height=3.5, width = 6)
for (tis in unique(Tissue)[1:7]){
  a=subset(out_fc2,method=="HOMER H3K27ac peak"&tissue==tis)$fc
  print(length(a[a>1])/length(a))
  b=subset(out_fc2,method=="HOMER H3K27ac weighted distance tissue, differential Z-score"&tissue==tis)$fc
  #print(length(b[b>1])/length(b))
  outtext<-paste("HOMER",tis,t.test(x=b,y=a)$p.value,t.test(x=b,y=a)$stat)
  #write.table(outtext,file="./figures/pvalues.txt",sep = '\t',quote=F,row.names = F,col.names = F,append = T)
}






