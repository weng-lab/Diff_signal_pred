#/usr/bin/env python
import os, gzip, sys
from multiprocessing import Process, Pool
from bx.bbi.bigwig_file import BigWigFile
sh = os.system

MARKS = ("dnase", "H3K4me1","H3K4me2","H3K4me3","H3K27ac","H3K9ac")
TISSUE = ("cranioface", "hindbrain", "limb","midbrain", "neural_tube","forebrain","liver","heart")
SOFTWARE = ("macs2","music","fseq","homer","hotspot","rseg","finder","bcp","mosaics","dfilter")


def refine_with_summit(_soft,_mark,_tissue):
    _temp_peak = [i.rstrip().split('\t') for i in open("/Data/adam/dnase/sort_bed/{0}.{1}.{2}.bed"\
                                                       .format(_soft,_mark,_tissue))]

    _temp_bw = open("/Data/adam/dnase/bigwig/{0}.{1}.rep0.bw".format(_mark,_tissue))
    _temp_enrich = open("/Data/adam/dnase/enrich_all_bed/{0}.{1}.{2}.bed".format(_soft,_mark,_tissue),'w')
    _bw = BigWigFile(file=_temp_bw)

    for line in _temp_peak:
        vals = _bw.get(line[0],int(line[1]),int(line[2]))
        vals =tuple(vals)
        if len(vals)>0:
            maxs = 0
            for _key in vals:
                if float(_key[2])>maxs:
                    maxs = float(_key[2])
                    summit = _key[:2]
            summit_p=int((float(summit[0])+float(summit[1]))/2)
            if summit_p-1000>0:
                print >> _temp_enrich, "{0}\t{1}\t{2}".format(line[0],str(summit_p-1000),str(summit_p+999))
            else:
                print >> _temp_enrich, "{0}\t{1}\t{2}".format(line[0],1,2000)
    _temp_enrich.close()

    awk_args='{printf "%s\\t%s\\n", $0,NR}'
    sh("sort -k 1,1 -k 2g,2g /Data/adam/dnase/enrich_all_bed/{0}.{1}.{2}.bed| bedtools merge -i stdin\
     | awk '{3}'>/Data/adam/dnase/enrich_all_merge_bed/{0}.{1}.{2}.bed".format(_soft,_mark,_tissue,awk_args))
    enhancer_dir ="/Data/adam/dnase/enhancer/tissue_enhancer/"
    sh("bash /Data/adam/dnase/enhancer/roc_pr.sh {0}{1}_enhancer.txt {0}negative_enhancer.txt \
    {2} {3} {1}".format(enhancer_dir, _tissue,_soft,_mark))
    raw_pr = [i.rstrip().split('\t') for i in open("/Data/adam/dnase/roc_pr_value/{0}.{1}.{2}.bed"\
                                                   .format(_soft,_mark,_tissue))]
    pr_refine = []
    temp_positive = [0,0]
    temp_negative = [0,0]
    for _line in raw_pr:
        if _line[3]==0 and _line[4]==0:
            temp_positive = _line[1:3]
            out_line = _line[:]
            out_line[3:5] = temp_negative
            precision = float(out_line[1])/(float(out_line[1])+float(out_line[3]))
            out_line.append(str(precision))
            pr_refine.append(out_line)
        elif _line[1]==0 and _line[2]==0:
            temp_negative = _line[3:5]
            out_line = _line[:]
            out_line[1:3] = temp_positive
            precision = float(out_line[1])/(float(out_line[1])+float(out_line[3]))
            out_line.append(str(precision))
            pr_refine.append(out_line)
        else:
            print "error in {0}.{1}.{2}, {3}".format(_soft,_mark,_tissue,_line)
    pr_refine2 = ['\t'.join(i) for i in pr_refine]
    with open("/Data/adam/dnase/roc_pr_final/{0}.{1}.{2}.bed".format(_soft,_mark,_tissue),'w') as f:
        for item in pr_refine2:
            print >>f, item


def Get_summit(_file):
    _soft, _mark, _tissue = _file.split('/')[1:4]
    _temp = [i.rstrip() for i in gzip.open(_file)]
    if len(_temp)<2000:
        print "skip {0}".format(_file)
    else:
        if _soft=="bcp":
            sh("zcat {0} |sort -k 5gr,5gr | cut -f1-3 >/Data/adam/dnase/sort_bed/{1}.{2}.{3}.bed"\
               .format(_file,_soft,_mark,_tissue))
            refine_with_summit(_soft,_mark,_tissue)
        elif _soft=="dfilter":
            sh("zcat {0} |sort -k 7gr,7gr | cut -f1-3 >/Data/adam/dnase/sort_bed/{1}.{2}.{3}.bed"\
               .format(_file,_soft,_mark,_tissue))
            refine_with_summit(_soft,_mark,_tissue)
        elif _soft=="finder":
            sh("zcat {0} |sort -k 4gr,4gr | cut -f1-3 >/Data/adam/dnase/sort_bed/{1}.{2}.{3}.bed"\
               .format(_file,_soft,_mark,_tissue))
            refine_with_summit(_soft,_mark,_tissue)
        elif _soft=="fseq":
            sh("zcat {0} |sort -k 7gr,7gr | cut -f1-3 >/Data/adam/dnase/sort_bed/{1}.{2}.{3}.bed"\
               .format(_file,_soft,_mark,_tissue))
            refine_with_summit(_soft,_mark,_tissue)
        elif _soft=="homer":
            sh("zcat {0} |sort -k 5gr,5gr | cut -f1-3 >/Data/adam/dnase/sort_bed/{1}.{2}.{3}.bed"\
               .format(_file,_soft,_mark,_tissue))
            refine_with_summit(_soft,_mark,_tissue)
        elif _soft=="hotspot":
            sh("zcat {0} |sort -k 5g,5g -k 4gr,4gr | cut -f1-3 >/Data/adam/dnase/sort_bed/{1}.{2}.{3}.bed"\
               .format(_file,_soft,_mark,_tissue))
            refine_with_summit(_soft,_mark,_tissue)
        elif _soft=="macs2":
            if "broad" in _file:
                sh("zcat {0} |sort -k 8gr,8gr | cut -f1-3 >/Data/adam/dnase/sort_bed/{1}.broad.{2}.{3}.bed"\
                   .format(_file,_soft,_mark,_tissue))
                refine_with_summit("{0}.broad".format(_soft),_mark,_tissue)
            elif "narrow" in _file:
                sh("zcat {0} |sort -k 8gr,8gr| cut -f1-3 >/Data/adam/dnase/sort_bed/{1}.narrow.{2}.{3}.bed"\
                   .format(_file,_soft,_mark,_tissue))
                refine_with_summit("{0}.narrow".format(_soft),_mark,_tissue)
            elif "gapped" in _file:
                sh("zcat {0} |sort -k 14gr,14gr| cut -f1-3 >/Data/adam/dnase/sort_bed/{1}.gapped.{2}.{3}.bed"\
                   .format(_file,_soft,_mark,_tissue))
                refine_with_summit("{0}.gapped".format(_soft),_mark,_tissue)
            else:
                pass
        elif _soft=="music":
            if "broad" in _file:
                sh("zcat {0} |sort -k 8g,8g| cut -f1-3 >/Data/adam/dnase/sort_bed/{1}.broad.{2}.{3}.bed"\
                   .format(_file,_soft,_mark,_tissue))
                refine_with_summit("{0}.broad".format(_soft),_mark,_tissue)
            elif "punctuate" in _file:
                sh("zcat {0} |sort -k 8g,8g | cut -f1-3 >/Data/adam/dnase/sort_bed/{1}.punctuate.{2}.{3}.bed"\
                   .format(_file,_soft,_mark,_tissue))
                refine_with_summit("{0}.punctuate".format(_soft),_mark,_tissue)
            else:
                pass
        elif _soft=="mosaics":
            sh("zcat {0} |sort -k 7gr,7gr| cut -f1-3 >/Data/adam/dnase/sort_bed/{1}.{2}.{3}.bed"\
               .format(_file,_soft,_mark,_tissue))
            refine_with_summit(_soft,_mark,_tissue)
        elif _soft=="rseg":
            sh("zcat {0} |sort -k 6gr,6gr| cut -f1-3 >/Data/adam/dnase/sort_bed/{1}.{2}.{3}.bed"\
               .format(_file,_soft,_mark,_tissue))
            refine_with_summit(_soft,_mark,_tissue)
        else:
            print "CANNOT FIND SOFTWARE"
            sys.exit(1)

if __name__=='__main__':
    os.chdir("/home/adam/dnase/encode/data_out/")
    if not os.path.exists("/Data/adam/dnase/roc_pr_value/"):
        sh("mkdir /Data/adam/dnase/roc_pr_value/")
    if not os.path.exists("/Data/adam/dnase/roc_pr_final/"):
        sh("mkdir /Data/adam/dnase/roc_pr_final/")
    if not os.path.exists("/Data/adam/dnase/sort_bed/"):
        sh("mkdir /Data/adam/dnase/sort_bed/")
    if not os.path.exists("/Data/adam/dnase/enrich_all_bed/"):
        sh("mkdir /Data/adam/dnase/enrich_all_bed/")
    if not os.path.exists("/Data/adam/dnase/enrich_all_merge_bed/"):
        sh("mkdir /Data/adam/dnase/enrich_all_merge_bed/")
    pool = Pool(processes=5)
    result =[]
    bedlist = [i.rstrip() for i in open('/Data/adam/dnase/top_bed/final_peak_list.txt')]
    print "bedlist",bedlist[0]
    for item in bedlist:
        result.append(pool.apply_async(Get_summit,(item,)))

    pool.close()
    pool.join()
