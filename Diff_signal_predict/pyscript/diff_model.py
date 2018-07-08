#!/usr/bin/env python
# coding:utf-8
import numpy as np
import os
from scipy.stats.mstats import zscore

def sh(args):
    # f = open("cmd.txt","w")
    # print >> f, args
    return os.system(args)


class Diff_model(object):

    def __init__(self,expr_file, peak_file,enhancer="",outpre="./"):
        self.tissue = "input"
        self.enhancer = enhancer
        self.validate = ("cranioface_11.5_day","limb_11.5_day","forebrain_11.5_day","midbrain_11.5_day",
                         "hindbrain_11.5_day","neural_tube_11.5_day","liver_11.5_day","heart_11.5_day")
        self.expr_file = expr_file
        self.peak_file = peak_file
        self.outprefix = outpre
        print "Start differential signal recalibration."

    def search_db(self):
        myf = [i.rstrip().split("\t") for i in open(self.expr_file)]
        my_expr = myf[1:]
        my_tis = myf[0]
        my_expr_trans = map(list, zip(*my_expr))
        my_expr_db = my_expr_trans[1:]
        my_expr_input = my_expr_trans[0]

        for index, item in enumerate(my_expr_db):
            if item == my_expr_input :
                print "The input tissue is {}.".format(my_tis[index+1])
                self.tissue = my_tis[index+1]
                break
        if self.tissue == 'input':
            print "The tissue is not found in our backend."

    def weight_fc(self):
        sh("Rscript $DIFF_PRED/rscript/weighted_fc.r {0} {1}".format(self.expr_file,self.peak_file))
        sh("sort -k 4gr,4gr {0}.weightfc.txt >{0}.weightfc.sorted.txt".format(self.peak_file))
        if self.tissue in self.validate:
            self.enhancer = "$DIFF_PRED/lib/enhancer_db/{0}_enhancer.txt".format(self.tissue)
        if self.enhancer:
            self._count_pca_fc()
            self._prcurve_plot(self.peak_file, "{0}.weightfc.sorted.txt".format(self.peak_file),
                               self.enhancer)

    def weight_zscore(self):
        sh("Rscript $DIFF_PRED/rscript/weighted_zscore.r {0} {1}".format(self.expr_file, self.peak_file))
        sh("sort -k 4gr,4gr {0}.weightzscore.txt >{0}.weightzscore.sorted.txt".format(self.peak_file))
        if self.tissue in self.validate:
            self.enhancer = "$DIFF_PRED/lib/enhancer_db/{0}_enhancer.txt".format(self.tissue)
        if self.enhancer:
            self._count_pca_zscore()
            self._prcurve_plot(self.peak_file, "{0}.weightzscore.sorted.txt".format(self.peak_file),
                               self.enhancer)

    # sub main function
    def _prcurve_plot(self, rawpeak, adjustpeak, enhancers):
        # select only top 5
        # self._count_pruac(rawpeak, enhancers, "raw",self.outprefix)
        self._count_pruac(adjustpeak, enhancers, "adjusted",self.outprefix)
        sh("mkdir -p {0}/prauc_tmp/plotdir; cp {0}/pca_tmp/pca_pr/inputpr_table.txt {0}/prauc_tmp/plotdir/raw ;\
           cp {0}/adjustedpr_table.txt {0}/prauc_tmp/plotdir/adjusted".format(self.outprefix))
        self._plt_prcurve("{0}/prauc_tmp/plotdir/".format(self.outprefix),self.outprefix)

    def _count_pruac(self, peak, enhancers, prefix, outdir):
        sh("bash $DIFF_PRED/shscript/prauc_step1.sh {0} {1} {2} {3}".format(peak, enhancers, prefix, outdir))
        # all files in prauc_tmp folder
        # print "Prefix is {}".format(prefix)
        sh("Rscript $DIFF_PRED/rscript/prauc_step2.r {0} {1}".format(prefix, outdir))

    def _count_pca_fc(self):
        _outdir = "{0}/pca_tmp/".format(self.outprefix)
        sh("mkdir -p {0}".format(_outdir))
        sh("mkdir -p {0}/pca_peak/".format(_outdir))
        sh("mkdir -p {0}/pca_pr/".format(_outdir))
        if os.path.isfile("{0}/pca_pr/prauc_val.txt".format(_outdir)):
            sh("rm {0}/pca_pr/prauc_val.txt".format(_outdir))
        sh("Rscript $DIFF_PRED/rscript/pca_fc.r {0} {1} {2}"\
           .format(self.expr_file, self.peak_file, _outdir))
        peakf = os.listdir("{0}/pca_peak/".format(_outdir))
        self._count_pruac("{0}/pca_peak/tissue0peak.txt".format(_outdir),
                          self.enhancer, "input",
                          "{0}/pca_pr/".format(_outdir))

        for index, item in enumerate(peakf):
            self._count_pruac("{0}/pca_peak/tissue{1}peak.txt".format(_outdir,index+2),
                              self.enhancer,"tis{0}".format(str(index+2)),
                              "{0}/pca_pr/".format(_outdir))
        sh("mkdir -p {0}/plotdir/ ;cp {0}/distance.txt {0}/plotdir/distance.txt; \
            cp {0}/pca_pr/prauc_val.txt {0}/plotdir/prauc_val.txt".format(_outdir))
        # output the figure to main directory
        self._plt_scatter("{0}/plotdir/".format(_outdir), self.outprefix)

    def _count_pca_zscore(self):
        _outdir = "{0}/pca_tmp/".format(self.outprefix)
        sh("mkdir -p {0}".format(_outdir))
        sh("mkdir -p {0}/pca_peak/".format(_outdir))
        sh("mkdir -p {0}/pca_pr/".format(_outdir))
        if os.path.isfile("{0}/pca_pr/prauc_val.txt".format(_outdir)):
            sh("rm {0}/pca_pr/prauc_val.txt".format(_outdir))
        sh("Rscript $DIFF_PRED/rscript/pca_zscore.r {0} {1} {2}"\
           .format(self.expr_file, self.peak_file, _outdir))

        peakf = os.listdir("{0}/pca_peak/".format(_outdir))
        self._count_pruac("{0}/pca_peak/tissue0peak.txt".format(_outdir),
                         self.enhancer, "input",
                         "{0}/pca_pr/".format(_outdir))

        for index, item in enumerate(peakf):
            self._count_pruac("{0}/pca_peak/tissue{1}peak.txt".format(_outdir,index+2),
                              self.enhancer,"tis{0}".format(str(index+2)),
                              "{0}/pca_pr/".format(_outdir))
        sh("mkdir -p {0}/plotdir/ ;cp {0}/distance.txt {0}/plotdir/distance.txt; \
           cp {0}/pca_pr/prauc_val.txt {0}/plotdir/prauc_val.txt".format(_outdir))
        # output the figure to main directory
        self._plt_scatter("{0}/plotdir/".format(_outdir),self.outprefix)


    def _plt_scatter(self,pltdir, outdir):
        sh("Rscript $DIFF_PRED/rscript/plot_scatter.r {0} {1}/improvement.pdf".format(pltdir, outdir))

    # sub function
    def _plt_prcurve(self, pltdir, outdir):
        sh("Rscript $DIFF_PRED/rscript/plot_pr_curve.r {0} {1}".format(pltdir, outdir))

















# python code for pca
def zeroMean(dataMat):
    meanVal=np.mean(dataMat,axis=0)
    newData=dataMat-meanVal
    return newData,meanVal

def percentage2n(eigVals,percentage):
    sortArray=np.sort(eigVals)   #升序
    sortArray=sortArray[-1::-1]  #逆转，即降序
    arraySum=sum(sortArray)
    tmpSum=0
    num=0
    for i in sortArray:
        tmpSum+=i
        num+=1
        if tmpSum>=arraySum*percentage:
            return num


def pca(dataMat,percentage=0.80):
    newData,meanVal=zeroMean(dataMat)
    covMat=np.cov(newData,rowvar=0)    #求协方差矩阵,return ndarray；若rowvar非0，一列代表一个样本，为0，一行代表一个样本
    eigVals,eigVects=np.linalg.eig(np.mat(covMat))#求特征值和特征向量,特征向量是按列放的，即一列代表一个特征向量
    n=percentage2n(eigVals,percentage)                 #要达到percent的方差百分比，需要前n个特征向量
    eigValIndice=np.argsort(eigVals)            #对特征值从小到大排序
    n_eigValIndice=eigValIndice[-1:-(n+1):-1]   #最大的n个特征值的下标
    n_eigVect=eigVects[:,n_eigValIndice]        #最大的n个特征值对应的特征向量
    lowDDataMat=newData*n_eigVect               #低维特征空间的数据
    reconMat=(lowDDataMat*n_eigVect.T)+meanVal  #重构数据
    return lowDDataMat,reconMat
