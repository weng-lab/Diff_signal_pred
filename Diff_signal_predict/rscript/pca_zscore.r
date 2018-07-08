library(DESeq2)

args = commandArgs(T)
expr_file = args[1]
peak_file = args[2]
outdir = args[3]

test<-read.table(expr_file,sep='\t',header=TRUE)
condition=colnames(test)
countdata<-floor(test*10)+1
rownames(countdata)<-paste("peak_",c(1:nrow(countdata)),sep='')


coldata<- data.frame(row.names = colnames(countdata), condition)
dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~condition)
pca_dds <- estimateSizeFactors(dds)
se <- SummarizedExperiment(log2(counts(pca_dds, normalized=TRUE)),colData=colData(pca_dds))

# mydat<-plotPCA( DESeqTransform(se),returnData =TRUE)
# mydist<-as.matrix(dist(mydat[,c(1:6)]))
mydat<-prcomp(t(counts(pca_dds, normalized=TRUE)),center=T,scale.=T)
mydist <- as.matrix(dist(mydat$x[,c(1:10)]))

outd <- paste(outdir,"/distance.txt",sep='')
write.table(mydist[,1],file=outd,sep='\t',row.names=F,col.names=F,quote=F)

mycount<-scale(log2(countdata),center=T,scale.=T)

myfc<-mycount[,1]
test2<-read.table(peak_file,sep='\t')
outtab<-cbind(test2[,c(1:3)],myfc)
colnames(outtab)<-c("chrom","start","end","val")
outtab<-outtab[order(outtab$val, decreasing=TRUE),]
outf<-paste(outdir,"/pca_peak/tissue","0","peak.txt",sep='')
for (i in c(2:ncol(countdata))){
    myfc<-mycount[,1]-mycount[,i]
    test2<-read.table(peak_file,sep='\t')
    outtab<-cbind(test2[,c(1:3)],myfc)
    colnames(outtab)<-c("chrom","start","end","val")
    outtab<-outtab[order(outtab$val, decreasing=TRUE),]
    outf<-paste(outdir,"/pca_peak/tissue",i,"peak.txt",sep='')
    # outf的第四列为重新rank后的score。
    write.table(outtab,file=outf, row.names=F, col.names=F, quote=F, sep='\t')
}
