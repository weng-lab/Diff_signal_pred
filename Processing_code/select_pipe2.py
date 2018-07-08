#/usr/bin/env python
import os, gzip, sys
from multiprocessing import Process, Pool
from bw_summit import get_summit

dfilter_dir = "/Data/Peakcalling/data_out/dfilter/dnase"
hotspot_dir = "/Data/Peakcalling/data_out/hotspot/dnase"
homer_dir = "/Data/Peakcalling/data_out/homer/dnase"
macs2_dir = "/Data/Peakcalling/data_out/macs2/dnase"
john_dir = "/Data/Peakcalling/ENCODE_PEAK/dnase"


out_dir = "/Data/Peakcalling/pr_calling/select_peak"
TISSUE = ("cranioface", "hindbrain", "limb","midbrain", "neural_tube")
HISTONE = ("H3K27ac","H3K4me1","H3K4me2","H3K4me3","H3K9ac")

def sh(cmds):
    os.system(cmds)


def get_sorted_bed():
    for tissue in TISSUE:
        for reps in ["final","rep1","rep2"]:
            final_file = "{0}/{1}/dnase.{1}.{2}.bed.gz".format(hotspot_dir, tissue,reps)
            tmp_file="{2}/tmp_peak/hotspot.dnase.{0}.{1}.bed"\
                .format(tissue, reps,out_dir)
            sh('zcat {0} | sort -k 5g,5g -k 4gr,4gr | cut -f1-3 | grep -v "chrM" >{1}'.format(final_file, tmp_file))
            # dfilter
            final_file = "{0}/{1}/dnase.{1}.{2}.bed.gz".format(dfilter_dir, tissue,reps)
            tmp_file="{2}/tmp_peak/dfilter.dnase.{0}.{1}.bed"\
                .format(tissue, reps,out_dir)
            sh('zcat {0} |sort -k 6gr,6gr -k 4gr,4gr | cut -f1-3 | grep -v "chrM" >{1}'.format(final_file, tmp_file))
            #
            # homer
            final_file = "{0}/{1}/dnase.{1}.{2}.bed.gz".format(homer_dir, tissue,reps)
            tmp_file="{2}/tmp_peak/homer.dnase.{0}.{1}.bed"\
                .format(tissue, reps,out_dir)
            sh('zcat {0} | sort -k 5gr,5gr | cut -f1-3 | grep -v "chrM" >{1}'.format(final_file, tmp_file))
            # # macs2
            final_file = "{0}/{1}/dnase.{1}.{2}.narrowPeak.gz".format(macs2_dir, tissue,reps)
            tmp_file="{2}/tmp_peak/macs2.dnase.{0}.{1}.bed"\
                .format(tissue, reps,out_dir)
            sh('zcat {0} | sort -k 8gr,8gr | cut -f1-3 | grep -v "chrM" >{1}'.format(final_file, tmp_file))


    print "sort finished"


def get_sorted_histone():
    # appending module
    # for tissue in TISSUE:
    for tissue in TISSUE+("liver","heart","forebrain"):
        for histone in HISTONE:

            for reps in ["final","rep1","rep2","rep0"]:

                # hotspot
                final_file = "/Data/Peakcalling/data_out/hotspot/{2}/{0}/{2}.{0}.{1}.bed.gz".format(tissue,reps,histone)
                tmp_file="{2}/tmp_peak/hotspot.{3}.{0}.{1}.bed"\
                    .format(tissue, reps,out_dir,histone)
                sh('zcat {0} | sort -k 5g,5g -k 4gr,4gr | cut -f1-3 | grep -v "chrM" >{1}'.format(final_file, tmp_file))
                # dfilter
                final_file = "/Data/Peakcalling/data_out/dfilter/{2}/{0}/{2}.{0}.{1}.bed.gz".format(tissue,reps,histone)
                tmp_file="{2}/tmp_peak/dfilter.{3}.{0}.{1}.bed"\
                    .format(tissue, reps,out_dir,histone)
                sh('zcat {0} |sort -k 6gr,6gr -k 4gr,4gr | cut -f1-3 | grep -v "chrM" >{1}'.format(final_file, tmp_file))

                # homer
                final_file = "/Data/Peakcalling/data_out/homer/{2}/{0}/{2}.{0}.{1}.bed.gz".format(tissue,reps,histone)
                tmp_file="{2}/tmp_peak/homer.{3}.{0}.{1}.bed"\
                    .format(tissue, reps,out_dir,histone)
                sh('zcat {0} | sort -k 5gr,5gr | cut -f1-3 | grep -v "chrM" >{1}'.format(final_file, tmp_file))

                # macs2
                final_file = "/Data/Peakcalling/data_out/macs2/{2}/{0}/{2}.{0}.{1}.narrowPeak.gz".format(tissue,reps,histone)
                tmp_file="{2}/tmp_peak/macs2.{3}.{0}.{1}.bed"\
                    .format(tissue, reps,out_dir,histone)
                sh('zcat {0} | sort -k 8gr,8gr | cut -f1-3 | grep -v "chrM" >{1}'.format(final_file, tmp_file))
                # final_file = "{0}/{1}/{3}.{1}.{2}.broadPeak.gz".format(macs2_dir, tissue,reps,histone)
                # tmp_file="{2}/tmp_peak/macs2.broad.{3}.{0}.{1}.bed"\
                #     .format(tissue, reps,out_dir,histone)
                # sh('zcat {0} | sort -k 8gr,8gr | cut -f1-3 | grep -v "chrM" >{1}'.format(final_file, tmp_file))
                # final_file = "{0}/{1}/{3}.{1}.{2}.gappedPeak.gz".format(macs2_dir, tissue,reps,histone)
                # tmp_file="{2}/tmp_peak/macs2.gapped.{3}.{0}.{1}.bed"\
                #     .format(tissue, reps,out_dir,histone)
                # sh('zcat {0} | sort -k 14gr,14gr | cut -f1-3 | grep -v "chrM" >{1}'.format(final_file, tmp_file))
                #music
                final_file = "/Data/Peakcalling/data_out/music/{2}/{0}/{2}.{0}.{1}.broad.bed.gz".format(tissue,reps,histone)
                tmp_file="{2}/tmp_peak/music.broad.{3}.{0}.{1}.bed"\
                    .format(tissue, reps,out_dir,histone)
                sh('zcat {0} | sort -k 8g,8g | cut -f1-3 | grep -v "chrM" >{1}'.format(final_file, tmp_file))
                final_file = "/Data/Peakcalling/data_out/music/{2}/{0}/{2}.{0}.{1}.punctuate.bed.gz".format(tissue,reps,histone)
                tmp_file="{2}/tmp_peak/music.punctuate.{3}.{0}.{1}.bed"\
                    .format(tissue, reps,out_dir,histone)
                sh('zcat {0} | sort -k 8g,8g | cut -f1-3 | grep -v "chrM" >{1}'.format(final_file, tmp_file))
                # BCP
                final_file = "/Data/Peakcalling/data_out/bcp/{2}/{0}/{2}.{0}.{1}_results_HM.bed.gz".format(tissue,reps,histone)
                tmp_file="{2}/tmp_peak/bcp.{3}.{0}.{1}.bed"\
                    .format(tissue, reps,out_dir,histone)
                sh('zcat {0} | sort -k 5gr,5gr | cut -f1-3 | grep -v "chrM" >{1}'.format(final_file, tmp_file))
                # FindER
                final_file = "/Data/Peakcalling/data_out/finder/{2}/{0}/{2}.{0}.{1}.bed.gz".format(tissue,reps,histone)
                tmp_file="{2}/tmp_peak/finder.{3}.{0}.{1}.bed"\
                    .format(tissue, reps,out_dir,histone)
                sh('zcat {0} | sort -k 4gr,4gr | cut -f1-3 | grep -v "chrM" >{1}'.format(final_file, tmp_file))
                # Fseq
                final_file = "/Data/Peakcalling/data_out/fseq/{2}/{0}/{2}.{0}.{1}.narrowPeak.gz".format(tissue,reps,histone)
                tmp_file="{2}/tmp_peak/fseq.{3}.{0}.{1}.bed"\
                    .format(tissue, reps,out_dir,histone)
                sh('zcat {0} | sort -k 7gr,7gr | cut -f1-3 | grep -v "chrM" >{1}'.format(final_file, tmp_file))
                # mosaics
                final_file = "/Data/Peakcalling/data_out/mosaics/{2}/{0}/{2}.{0}.{1}.bed.gz".format(tissue,reps,histone)
                tmp_file="{2}/tmp_peak/mosaics.{3}.{0}.{1}.bed"\
                    .format(tissue, reps,out_dir,histone)
                sh('zcat {0} | sort -k 7gr,7gr | cut -f1-3 | grep -v "chrM" >{1}'.format(final_file, tmp_file))
                # rseg
                final_file = "/Data/Peakcalling/data_out/rseg/{2}/{0}/{2}.{0}.{1}.bed.gz".format(tissue,reps,histone)
                tmp_file="{2}/tmp_peak/rseg.{3}.{0}.{1}.bed"\
                    .format(tissue, reps,out_dir,histone)
                sh('zcat {0} | sort -k 6gr,6gr | cut -f1-3 | grep -v "chrM" >{1}'.format(final_file, tmp_file))



def process_histone_pr():
    for tissue in TISSUE:
        for histone in HISTONE:
            for _soft in ["hotspot","homer","macs2","dfilter","fseq","bcp",
                          "music.broad","music.punctuate","mosaics","rseg","finder"]:
                pr_name ="{0}/2000bp/adjust_bed/{1}.{3}.{2}.final.pr.bed".format(out_dir, _soft, tissue, histone)
                if os.path.exists(pr_name):
                    pr_file =[i.rstrip().split('\t') for i in open(pr_name)]
                    pr_file2 = [i+[_soft,tissue,histone] for i in pr_file]
                    pr_file2 = ["\t".join(item) for item in pr_file2]
                    with open("{0}/histone_pr_summary.txt".format(out_dir),"a") as f:
                        for _line in pr_file2:
                            print >>f, _line
                pr_name2 ="{0}/2000bp/adjust_bed/{1}.{3}.{2}.final.combine.pr.bed".format(out_dir, _soft, tissue, histone)
                if os.path.exists(pr_name2):
                    pr_file =[i.rstrip().split('\t') for i in open(pr_name2)]
                    pr_file2 = [i+[_soft,tissue,"{0}.combine".format(histone)] for i in pr_file]
                    pr_file2 = ["\t".join(item) for item in pr_file2]
                    with open("{0}/histone_pr_summary.txt".format(out_dir),"a") as f:
                        for _line in pr_file2:
                            print >>f, _line


def process_promoter_pr():
    if os.path.exists("{0}/promoter_pr_summary.txt".format(out_dir)):
        sh("rm {0}/promoter_pr_summary.txt".format(out_dir))

    for tissue in TISSUE+("liver","heart","forebrain"):
        for histone in HISTONE:
            for _soft in ["hotspot","homer","macs2","dfilter","fseq","bcp",
                          "music.broad","music.punctuate","mosaics","rseg","finder"]:
                for reps in ["final","rep1","rep2"]:
                    pr_name ="{0}/promoter_pr/{1}.{3}.{2}.{4}.pr.bed".format(out_dir, _soft, tissue, histone,reps)
                    if os.path.exists(pr_name):
                        pr_file =[i.rstrip().split('\t') for i in open(pr_name)]
                        pr_file2 = [i+[_soft,tissue,histone,reps] for i in pr_file]
                        pr_file2 = ["\t".join(item) for item in pr_file2]
                        with open("{0}/promoter_pr_summary.txt".format(out_dir),"a") as f:
                            for _line in pr_file2:
                                print >>f, _line
    for tissue in TISSUE:
        for _soft in ["hotspot","homer","macs2","dfilter","hotspot2"]:
            # for reps in ["final","rep1","rep2"]:
            for reps in ["final"]:
                pr_name = "{0}/promoter_pr/{1}.dnase.{2}.{3}.pr.bed".format(out_dir, _soft, tissue, reps)
                if os.path.exists(pr_name):
                    pr_file = [i.rstrip().split('\t') for i in open(pr_name)]
                    pr_file2 = [i + [_soft, tissue, "dnase",reps] for i in pr_file]
                    pr_file2 = ["\t".join(item) for item in pr_file2]
                    with open("{0}/promoter_pr_summary.txt".format(out_dir), "a") as f:
                        for _line in pr_file2:
                            print >> f, _line


def process_pr():
    if os.path.exists("{0}/dnase_pr_summary.txt".format(out_dir)):
        sh("rm {0}/dnase_pr_summary.txt".format(out_dir))
    for tissue in TISSUE:
        for reps in ["final","rep1","rep2"]:
            for _soft in ["dfilter","hotspot","homer","john","macs2"]:
                for peaklen in [300,1000]:
                    pr_name ="{0}/{4}bp/adjust_bed/{1}.dnase.{2}.{3}.pr.bed".format(out_dir, _soft, tissue, reps,peaklen)
                    if os.path.exists(pr_name):
                        pr_file =[i.rstrip().split('\t') for i in open(pr_name)]
                        pr_file2 = [i+[_soft,tissue,reps,str(peaklen)] for i in pr_file]
                        pr_file2 = ["\t".join(item) for item in pr_file2]
                        with open("{0}/dnase_pr_summary.txt".format(out_dir),"a") as f:
                            for _line in pr_file2:
                                print >>f, _line


def process_adjust_dnase_pr():
    if os.path.exists("{0}/dnase_pr_adjust.txt".format(out_dir)):
        sh("rm {0}/dnase_pr_adjust.txt".format(out_dir))
    for tissue in TISSUE:
        pr_name ="{0}/300bp/adjust_bed/dfilter.dnase.{1}.final.combine.pr.bed".format(out_dir,tissue)
        if os.path.exists(pr_name):
            # print pr_name,"prname"
            pr_file =[i.rstrip().split('\t') for i in open(pr_name)]
            pr_file2 = [i+[tissue] for i in pr_file]
            pr_file2 = ["\t".join(item) for item in pr_file2]
            with open("{0}/dnase_pr_adjust.txt".format(out_dir),"a") as f:
                for _line in pr_file2:
                    print >>f, _line


def process_adjust_john_pr():
    if os.path.exists("{0}/john_pr_adjust.txt".format(out_dir)):
        sh("rm {0}/john_pr_adjust.txt".format(out_dir))
    for tissue in TISSUE:
        for rep in ["rep1","rep2"]:
            pr_name ="{0}/500bp/adjust_bed/john.dnase.{1}.{2}.combine.pr.bed".format(out_dir,tissue,rep)
            if os.path.exists(pr_name):
                pr_file =[i.rstrip().split('\t') for i in open(pr_name)]
                pr_file2 = [i+["dnase.combine",tissue,rep] for i in pr_file]
                pr_file2 = ["\t".join(item) for item in pr_file2]
                with open("{0}/john_pr_adjust.txt".format(out_dir),"a") as f:
                    for _line in pr_file2:
                        print >>f, _line
            pr_name2 ="{0}/500bp/adjust_bed/john.dnase.{1}.{2}.pr.bed".format(out_dir,tissue,rep)
            if os.path.exists(pr_name2):
                pr_file =[i.rstrip().split('\t') for i in open(pr_name2)]
                pr_file2 = [i+["dnase",tissue,rep] for i in pr_file]
                pr_file2 = ["\t".join(item) for item in pr_file2]
                with open("{0}/john_pr_adjust.txt".format(out_dir),"a") as f:
                    for _line in pr_file2:
                        print >>f, _line


def process_test_dfilter_pr():
    if os.path.exists("{0}/dfilter_test.txt".format(out_dir)):
        sh("rm {0}/dfilter_test.txt".format(out_dir))
    for tissue in TISSUE:
        pr_name ="{0}/test_dfilter2/dfilter.dnase.{1}.final.pr.bed".format(out_dir,tissue)
        if os.path.exists(pr_name):
            pr_file =[i.rstrip().split('\t') for i in open(pr_name)]
            pr_file2 = [i+[tissue] for i in pr_file]
            pr_file2 = ["\t".join(item) for item in pr_file2]
            with open("{0}/dfilter_test.txt".format(out_dir),"a") as f:
                for _line in pr_file2:
                    print >>f, _line


def process_dnase_histone_pr():
    if os.path.exists("{0}/histone_dnase_pr.txt".format(out_dir)):
        sh("rm {0}/histone_dnase_pr.txt".format(out_dir))
    for _soft in ["hotspot","homer","macs2","dfilter","fseq","bcp",
                              "music.broad","music.punctuate","mosaics","rseg","finder"]:
        for tissue in TISSUE:
            pr_name ="{0}/2000bp/histone_dnase_pr/{2}.H3K27ac_dnase.{1}.final.dnase_hist.pr.bed".format(out_dir,tissue,_soft)
            if os.path.exists(pr_name):
                pr_file =[i.rstrip().split('\t') for i in open(pr_name)]
                pr_file2 = [i+[_soft,"H3K27ac",tissue] for i in pr_file]
                pr_file2 = ["\t".join(item) for item in pr_file2]
                with open("{0}/histone_dnase_pr.txt".format(out_dir),"a") as f:
                    for _line in pr_file2:
                        print >>f, _line


def process_macs2_dnase_pr():
    if os.path.exists("{0}/dfilter_adjust_macs_pr.txt".format(out_dir)):
        sh("rm {0}/dfilter_adjust_macs_pr.txt".format(out_dir))
    for tissue in TISSUE:
        for soft in ["dfilter_origin","dfilter_macs2"]:
            pr_name ="{0}/origin_peak/adjust_bed/{1}.dnase.{2}.final.pr.bed".format(out_dir,soft,tissue)
            if os.path.exists(pr_name):
                pr_file =[i.rstrip().split('\t') for i in open(pr_name)]
                pr_file2 = [i+[soft,tissue] for i in pr_file]
                pr_file2 = ["\t".join(item) for item in pr_file2]
                with open("{0}/dfilter_adjust_macs_pr.txt".format(out_dir),"a") as f:
                    for _line in pr_file2:
                        print >>f, _line


def get_enrichment():
    if os.path.exists("/Data/Peakcalling/pr_calling/select_peak/dnase_enrichment.txt"):
        sh("rm /Data/Peakcalling/pr_calling/select_peak/dnase_enrichment.txt")
    sh("mkdir /Data/Peakcalling/pr_calling/select_peak/tmp/")
    pool =Pool(processes=15)
    result = []
    for tissue in TISSUE:
        for width in ["150","300","500","1000"]:
            for _soft in ["hotspot","homer","dfilter","macs2"]:
                peaks = "/Data/Peakcalling/pr_calling/select_peak/{2}bp/adjust_bed/{0}.dnase.{1}.final.adjust.bed"\
                .format(_soft, tissue,width)
                alignments = "/Data/Peakcalling/tagAlign_dir/dnase.{0}.rep0.tagAlign.gz".format(tissue)
                outfile = "/Data/Peakcalling/pr_calling/select_peak/dnase_enrichment.txt"
                result.append(pool.apply_async(sh,("bash ./bin/enrichment.sh {0} {1} {2} {3} {4}"\
                                               .format(peaks,alignments, _soft, tissue,width),)))


    pool.close()
    pool.join()
    sh("cat /Data/Peakcalling/pr_calling/select_peak/tmp/*.txt > /Data/Peakcalling/pr_calling/select_peak/dnase_enrichment.txt")
    sh("rm -rf /Data/Peakcalling/pr_calling/select_peak/tmp/")


def get_peak_statistic():
    os.chdir("/Data/Peakcalling/pr_calling/select_peak/tmp_peak")
    file_list=os.listdir("./")
    if os.path.exists("/Data/Peakcalling/pr_calling/select_peak/thesis_count_summary.txt"):
        sh("rm /Data/Peakcalling/pr_calling/select_peak/thesis_count_summary.txt")
    if os.path.exists("/Data/Peakcalling/pr_calling/select_peak/thesis_length_summary.txt"):
        sh("rm /Data/Peakcalling/pr_calling/select_peak/thesis_length_summary.txt")

    for _file in file_list:
        if not _file.startswith("."):
            if "punctuate" in _file or "broad" in _file:
                mark, tissue, rep = _file.split(".")[2:5]
                soft="{0}_{1}".format(_file.split(".")[0],_file.split(".")[1])
            else:
                soft, mark, tissue, rep = _file.split(".")[:4]
            sh("bash /Data/Peakcalling/pr_calling/scripts/bin/thesis_count_length.sh {0} {1} {2} {3} {4}"\
               .format(_file,soft, mark, tissue, rep))




def process_mouse_pr():

    if not os.path.exists("/Data/Peakcalling/pr_calling/select_peak/tmp_peak"):
        os.system("mkdir /Data/Peakcalling/pr_calling/select_peak/tmp_peak")
    if not os.path.exists("/Data/Peakcalling/pr_calling/select_peak/tmp_peak2"):
        os.system("mkdir /Data/Peakcalling/pr_calling/select_peak/tmp_peak2")
    step1 process tagAlign to rank sort bed
    get_sorted_bed()
    get_sorted_histone()
    
    # step 2: refine summit and get bed file for pr
    pool =Pool(processes=10)
    result = []
    for tissue in TISSUE:
        for reps in ["final","rep1","rep2"]:
            for _soft in ["hotspot","homer","macs2"]:

                peakfile = "{0}/tmp_peak/{1}.dnase.{2}.{3}.bed".format(out_dir, _soft, tissue, reps)
                bwfile = "/Data/adam/dnase/bigwig/dnase.{0}.rep0.bw".format(tissue)
                if reps == "final":
                    bwfile = "/Data/adam/dnase/bigwig/dnase.{0}.rep0.bw".format(tissue)
                else:
                    bwfile = "/Data/adam/dnase/bigwig/dnase.{0}.{1}.bw".format(tissue, reps)
                outprefix = "{0}/tmp_peak2/{1}.dnase.{2}.{3}".format(out_dir, _soft, tissue, reps)
                if os.path.exists(peakfile):
                    result.append(pool.apply_async(get_summit, (peakfile,bwfile,outprefix, )))
    
    pool.close()
    pool.join()
    ###
    generate dfilter bigwig summit in "test" directory other than tmp_peak2
    ##
    if not os.path.exists("/Data/Peakcalling/pr_calling/select_peak/test_dfilter"):
        os.system("mkdir /Data/Peakcalling/pr_calling/select_peak/test_dfilter")
    if not os.path.exists("/Data/Peakcalling/pr_calling/select_peak/test_dfilter2"):
        os.system("mkdir /Data/Peakcalling/pr_calling/select_peak/test_dfilter2")
    for tissue in TISSUE:
    
        peakfile = "{0}/tmp_peak/dfilter.dnase.{1}.final.bed".format(out_dir,tissue)
    
        bwfile = "/Data/adam/dnase/bigwig/dnase.{0}.rep0.bw".format(tissue)
    
        outprefix = "{0}/test_dfilter/dfilter.dnase.{1}.final".format(out_dir, tissue)
        if os.path.exists(peakfile):
            get_summit(peakfile,bwfile,outprefix)
    for tissue in TISSUE:

        peakfile = "{0}/tmp_peak/dfilter.dnase.{1}.rep0.bed".format(out_dir,tissue)

        bwfile = "/Data/adam/dnase/bigwig/dnase.{0}.rep0.bw".format(tissue)

        outprefix = "{0}/test_dfilter/dfilter.dnase.{1}.rep0".format(out_dir, tissue)
        if os.path.exists(peakfile):
            get_summit(peakfile,bwfile,outprefix)

    # get summit of histone
    pool =Pool(processes=10)
    result = []

    for tissue in TISSUE+("liver","heart","forebrain"):
        for reps in ["final","rep1","rep2"]:
            for histone in HISTONE:
                for _soft in ["hotspot","homer","macs2","dfilter","fseq","bcp","music.broad","music.punctuate","mosaics","rseg","finder"]:
                    peakfile = "{0}/tmp_peak/{1}.{4}.{2}.{3}.bed".format(out_dir, _soft, tissue, reps,histone)
                    if reps == "final":
                        bwfile = "/Data/adam/dnase/bigwig/{1}.{0}.rep0.bw".format(tissue,histone)
                    else:
                        bwfile = "/Data/adam/dnase/bigwig/{2}.{0}.{1}.bw".format(tissue, reps,histone)
                    outprefix = "{0}/tmp_peak2/{1}.{4}.{2}.{3}".format(out_dir, _soft, tissue, reps,histone)
                    if os.path.exists(peakfile):
                        result.append(pool.apply_async(get_summit, (peakfile,bwfile,outprefix, )))
    
    pool.close()
    pool.join()
    
    # step 3: extend and merge and pr
    # get dfilter tmp2(summit file) for uniform processing
    for tissue in TISSUE:
        for reps in  ["rep1","rep2"]:

            inputfile = "{0}/{1}/dnase.{1}.{2}.bed.gz".format(dfilter_dir, tissue,reps)
            outputfile = "{0}/tmp_peak2/dfilter.dnase.{1}.{2}.bed".format(out_dir,tissue,reps)
            if os.path.exists(inputfile):
                sh('zcat {0} |sort -k 6gr,6gr -k 4gr,4gr | cut -f1,5 | grep -v "chrM" >{1}'.format(inputfile, outputfile))

    pool =Pool(processes=10)
    result = []
    for tissue in TISSUE:
        for reps in ["final","rep1","rep2"]:
            for _soft in ["dfilter","hotspot","homer","macs2"]:
                inputfile = "{0}/tmp_peak2/{1}.dnase.{2}.{3}.bed".format(out_dir, _soft, tissue, reps)
                if os.path.exists(inputfile):

                    for peaklen in [300,2000]:
                        outprefix = "{0}/{4}bp/adjust_bed/{1}.dnase.{2}.{3}".format(out_dir, _soft, tissue, reps,str(peaklen))
                        result.append(pool.apply_async(sh, ("bash ./bin/refine_call_pr.sh {0} {1} {2} {3}".format(inputfile, outprefix,tissue, peaklen/2),)))
    
    
    pool.close()
    pool.join()

    # process histone final files
    pool = Pool(processes=10)
    result = []
    for tissue in TISSUE+("liver","heart","forebrain"):
        for reps in ["rep0","final", "rep1", "rep2"]:
            for histone in HISTONE:
                for _soft in ["hotspot","homer","macs2","dfilter","fseq","bcp",
                              "music.broad","music.punctuate","mosaics","rseg","finder"]:
                    inputfile = "{0}/tmp_peak2/{1}.{4}.{2}.{3}.bed"\
                        .format(out_dir, _soft, tissue, reps,histone)
                    if os.path.exists(inputfile):
    
                        outprefix = "{0}/2000bp/adjust_bed/{1}.{4}.{2}.{3}"\
                            .format(out_dir, _soft, tissue, reps,histone)
                        result.append(pool.apply_async(sh, (
                        "bash ./bin/refine_call_pr.sh {0} {1} {2} {3}".format(inputfile, outprefix, tissue,
                                                                                  1000),)))
    
    pool.close()
    pool.join()


    process_macs2_dnase_pr()
    adjust john stam results by 500bp dnase and 2kb h3k27ac
    pool =Pool(processes=1)
    result = []
    for tissue in TISSUE[1:]:
        # for reps in ["rep1","rep2"]:
        for reps in ["rep1"]:
            inputfile = "{0}/tmp_peak2/john.dnase.{1}.{2}.bed".format(out_dir, tissue, reps)
            if os.path.exists(inputfile):
                histone_bw = "/Data/adam/dnase/bigwig/H3K27ac.{0}.rep0.bw".format(tissue)
                dnase_bw="/Data/adam/dnase/bigwig/dnase.{0}.{1}.bw".format(tissue, reps)
                outprefix = "{0}/500bp/adjust_bed/john.dnase.{1}.{2}".format(out_dir, tissue, reps)
                result.append(pool.apply_async(sh, ("bash ./bin/refine_john_pr.sh {0} {1} {2} {3} {4}"\
                                                    .format(inputfile, outprefix,tissue,histone_bw,dnase_bw),)))
    
    pool.close()
    pool.join()


    for tissue in TISSUE:
        inputfile = "{0}/test_dfilter/dfilter.dnase.{1}.final.bed".format(out_dir,tissue)
        if os.path.exists(inputfile):
            outprefix = "{0}/test_dfilter2/dfilter.dnase.{1}.final".format(out_dir,tissue)
            sh("bash ./bin/refine_call_pr.sh {0} {1} {2} 150".format(inputfile, outprefix,tissue))
    for tissue in TISSUE:
        inputfile = "{0}/test_dfilter/dfilter.dnase.{1}.rep0.bed".format(out_dir,tissue)
        if os.path.exists(inputfile):
            outprefix = "{0}/test_dfilter2/dfilter.dnase.{1}.rep0".format(out_dir,tissue)
            sh("bash ./bin/refine_call_pr.sh {0} {1} {2} 150".format(inputfile, outprefix, tissue))


    if not os.path.isdir("{0}/promoter_pr/".format(out_dir)):
    sh("mkdir {0}/promoter_pr/".format(out_dir))
    
    pool =Pool(processes=10)
    result = []
    for tissue in TISSUE+("liver","heart","forebrain"):
        # for reps in ["final", "rep1", "rep2"]:
        for reps in ["final"]:
            for histone in HISTONE:
                for _soft in ["hotspot","homer","macs2","dfilter","fseq","bcp",
                              "music.broad","music.punctuate","mosaics","rseg","finder"]:
                # for _soft in ["fseq"]:
                    outprefix = "{0}/promoter_pr/{2}.{3}.{1}.{4}".format(out_dir, tissue,_soft,histone,reps)
                    inputfile ="/home/adam/prauc/data/{2}.{3}.{1}.{4}.adjust.bed".format(out_dir,tissue,_soft,histone,reps)
                    if os.path.exists(inputfile):
                        pass
                        result.append(pool.apply_async(sh, ("bash ./bin/promoter_pr.sh {0} {1} {2} histone"\
                                                            .format(inputfile, outprefix,tissue),)))
    pool.close()
    pool.join()
    pool = Pool(processes=10)
    result = []
    for tissue in TISSUE:
        for reps in ["final", "rep1", "rep2"]:
            for _soft in ["hotspot","homer","dfilter","macs2","hotspot2"]:
                outprefix = "{0}/promoter_pr/{2}.dnase.{1}.{3}".format(out_dir, tissue, _soft, reps)
                inputfile = "/home/adam/prauc/data/{2}.dnase.{1}.{3}.adjust.bed".format(out_dir, tissue, _soft,reps)
                if os.path.exists(inputfile):

                    result.append(pool.apply_async(sh, ("bash ./bin/promoter_pr.sh {0} {1} {2}".format(inputfile, outprefix, tissue),)))
    pool.close()
    pool.join()
    process_promoter_pr()


if __name__ == "__main__":

    process_mouse_pr()
    get_enrichment()
    get_peak_statistic()
