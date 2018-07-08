from bx.bbi.bigwig_file import BigWigFile
import os,sys

def sh(args):
    return os.system(args)



def peak_process(peaks,bigwigs, mode, resize, outpre):
    blacklist_file = "mm10.blacklist.bed.gz"
    tss_file = "mm10_tss.bed"
    if not resize:
        awk_args = '{printf "%s\\t%d\\t%d\\t%d\\n",$1,$2,$3,NR}'
        # Need revise with more genome versions.
        ##
        sh("awk '{0}' {1} | grep -v '-' | intersectBed -a stdin -b $DIFF_PRED/lib/genome_file/{4} -v\
         | intersectBed -a stdin -b $DIFF_PRED/lib/genome_file/{3} -v > {2}/input.bed"\
           .format(awk_args,peaks,outpre,blacklist_file, tss_file))
    else:
        sh("mkdir -p {}/resize/".format(outpre))
        if mode=="DNase-seq":
            _width=150
            # _width = 1000
        elif mode=="ATAC-seq" or mode=="H3K27ac":
            _width=1000
        else:
            print "Unknown error in mode! Exit!"
            sys.exit(1)
        # _get_resized(peaks, bigwigs, _width, "{}/resize/raw".format(outpre))
        _get_resized2(peaks, _width, "{}/resize/raw".format(outpre))


        sh("export LC_ALL=C; sort -k1,1 -k 2g,2g -k 3g,3g {0}/resize/raw.bed| \
        mergeBed -i stdin -c 4 -o min | sort -k 4g,4g | intersectBed -a stdin -b $DIFF_PRED/lib/genome_file/{2} -v\
         | intersectBed -a stdin -b $DIFF_PRED/lib/genome_file/{1} -v| grep -v '-'| grep -v 'GL' >\
         {0}/resize/process.bed".format(outpre, blacklist_file, tss_file))

        _get_resized2("{}/resize/process.bed".format(outpre), _width
                    , "{}/input".format(outpre))

        #  Unknown Error!!
        #  Use mid point instead temporarily
        # _get_resized("{}/resize/process.bed".format(outpre), bigwigs, _width
        #             , "{}/input".format(outpre))


        # Delete after debugging
        # sh("rm -rf {}/resize".format(outpre))

def _get_resized2(peakfile, width, outprefix):
    peaklis = [i.rstrip().split('\t') for i in open(peakfile)]
    peaklist = _check_file(peaklis)
    with open("{0}.bed".format(outprefix), "w") as f:
        for _ids, line in enumerate(peaklist):
            out_line = [line[0], max(1,str((int(line[1])+int(line[2]))/2-width)),
                        str((int(line[1])+int(line[2]))/2+width-1), str(_ids+1)]
            outline2 = "\t".join(out_line)
            print >>f, outline2


def _get_resized(peakfile, bigwigs,width, outprefix):
    # Input is peak file name which is sort by rank
    peaklis =[i.rstrip().split('\t') for i in open(peakfile)]
    peaklist = _check_file(peaklis)
    _temp_bw = open(bigwigs)
    _bw = BigWigFile(file=_temp_bw)

    # print peaklist[0]
    with open("{0}.bed".format(outprefix), "w") as f:
        for _ids, line in enumerate(peaklist):
            vals = _bw.get(line[0],int(line[1]),int(line[2]))
            vals =tuple(vals)
            if len(vals)>0:
                maxs = 0
                for _key in vals:
                    if float(_key[2])>maxs:
                        maxs = float(_key[2])
                        summit = _key[:2]
                summit_p=int((float(summit[0])+float(summit[1]))/2)
                # out_summit.append([line[0], str(summit_p)])
                print >>f, "{0}\t{1}\t{2}\t{3}"\
                    .format(line[0], str(summit_p-width),str(summit_p+width-1),str(_ids+1))


def _check_file(peaklis):
    _len = len(peaklis[0])
    peaklist = []
    for _line in peaklis:
        if len(_line)!=_len:
            print("{0} has removed".format(_line))
        else:
            peaklist.append(_line)
    return peaklist