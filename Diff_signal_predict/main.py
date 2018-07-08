#!/usr/bin/env python

import argparse
import logging
import os
import sys
from multiprocessing import Pool
from pyscript.diff_model import Diff_model
from pyscript.peak_pro import peak_process




def sh(args):
    return os.system(args)


def prepare_parser():
    description = "%(prog)s -- rerank resized peaks by diff signal"
    usage = """Usage: %(prog)s <-p peak> <--bigwig/--bam bw/bam file> <-m mode> [-e enhancer file -r -o outdir]
            Example: %(prog)s -p limb.e11.5.bed --bigwig limb.e11.5.bw -m 1 -r -o test -d mydb_path
            """
    parser = argparse.ArgumentParser("Diffpeak", description=description, usage=usage)
    parser.add_argument('-p','--peak',help="Origin peak file.",dest="peak",required=True)
    # group = parser.add_mutually_exclusive_group(required=True)
    # group.add_argument('--bigwig',help="Origin bigwig file", action="store",dest="bigwig")
    # group.add_argument('--bam',help="Origin bam file",dest="bam",action="store")
    parser.add_argument('--bigwig',help="Origin bigwig file", action="store",dest="bigwig",
                        required=True)
    parser.add_argument('-m','--mode',help="Please select 1, 2 or 3. 1:ATAC-seq; 2:H3K27ac ChIP-seq; 3:DNase-seq",
                        dest="mode",action="store", type=int, required=True)
    parser.add_argument('-e','--enhancer',help="User defined experiment validated enhancer set",
                        dest="enhancer",action="store")
    # parser.add_argument('-s','--species',help="Choose the species, mouse only now",default="mouse",
    #                     dest="species")
    parser.add_argument('-d','--db',help="Absolute path of database",dest="db",default="db")
    parser.add_argument('-o','--outpre',help="Output prefix",dest="outpre",default="./")
    parser.add_argument('-r','--resize',help="Resize to same peak width. All peaks resize to\
                        300bp for DNase-seq and ATAC-seq, 2000bp for H3K27ac ChIP-seq,\
                        corresponding to signal summit of bigwig(--bigwig). Recommend for\
                        enhancer prediction. Default: False.", dest="resize",
                        default = False,action="store_true")
    return parser


def arg_validate(parser):
    parser = parser.parse_args()
    if int(parser.mode) != 1 and int(parser.mode) != 2 and int(parser.mode) != 3:
        logging.error("Invalid mode! Please choose mode in 1, 2 and 3!")
        sys.exit(1)
    return parser


def build_matrix(bigwig, outpre, mode, db):
    sh("mkdir -p {}/signal/".format(outpre))
    sh("mkdir -p {}/diff/".format(outpre))
    sh("bigWigAverageOverBed {1} {0}/input.bed stdout | sort -k1g | cut -f 5 > {0}/signal/input.bed".format(outpre, bigwig))
    db_info = [i.rstrip() for i in open("{}/{}/summary.info".format(db,mode))]
    result = []
    pool = Pool(processes=15)
    for _line in db_info:
        result.append(pool.apply_async(sh,("bigWigAverageOverBed {2}/{3}/{1}_rep0.bigWig {0}/input.bed stdout | sort -k1g | cut -f5 > {0}/signal/{1}.bed"\
            .format(outpre, _line, db, mode),)))
        # sh("bigWigAverageOverBed {2}/{3}/{1}_rep0.bigWig {0}/input.bed stdout | sort -k1g | cut -f5 > {0}/signal/{1}.bed"\
        #    .format(outpre, _line, db, mode))
    pool.close()
    pool.join()

    outbed = os.listdir("{0}/signal/".format(outpre))
    out = []
    _tmp = [i.rstrip() for i in open("{0}/signal/{1}".format(outpre, "input.bed"))]
    _tmp2 = ["input"] + _tmp
    out.append(_tmp2)
    for _line in outbed:
        if _line.endswith("bed") and _line!="input.bed":
            _tmp = [i.rstrip() for i in open("{0}/signal/{1}".format(outpre,_line))]

            _tmp2 = [_line[:-4]] + _tmp
            out.append(_tmp2)
    out2 = map(list, zip(*out))
    with open("{0}/diff/whole_table.txt".format(outpre),"w") as f:
        outtab = ["\t".join(i) for i in out2]
        for _line in outtab:
            print >> f, _line


def adjust_peak(mode, expr_file, inputfile,enhancer='', outpre="./"):
    mymodel = Diff_model(expr_file, inputfile,enhancer,outpre)
    mymodel.search_db()
    if mode ==3:
        mymodel.weight_fc()
    elif mode == 1 or mode == 2:
        mymodel.weight_zscore()
    else:
        logging.error("Invalid mode! Please choose mode in 1, 2 and 3!")
        sys.exit(1)


def main():
    parser = arg_validate(prepare_parser())
    MODE = ("NULL","ATAC-seq","H3K27ac","DNase-seq")
    sh("mkdir -p {}".format(parser.outpre))
    peak_process(parser.peak, parser.bigwig, MODE[int(parser.mode)], parser.resize, parser.outpre)
    build_matrix(parser.bigwig, parser.outpre, MODE[int(parser.mode)], parser.db)
    adjust_peak(int(parser.mode),"{0}/diff/whole_table.txt".format(parser.outpre),
                "{0}/input.bed".format(parser.outpre),parser.enhancer, parser.outpre)


if __name__=="__main__":
    main()

