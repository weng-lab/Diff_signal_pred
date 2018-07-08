#/usr/bin/env python

import os,sys,subprocess
from multiprocessing import Process, Pool
import lib.caller.macs as Macs
import lib.caller.music as Music
import lib.caller.fseq as Fseq
import lib.caller.homer as Homer
import lib.caller.bcp as Bcp
import lib.caller.dfilter as Dfilter
import lib.caller.hotspot as Hotspot
import lib.caller.rseg as Rseg
import lib.caller.finder as Finder
import lib.caller.mosaics as Mosaics
from lib.pythonIO import *
import lib.caller.spp as Spp

MARKS = ("dnase", "H3K4me1","H3K4me2","H3K4me3","H3K27ac","H3K9ac")
TISSUE = ("cranioface", "hindbrain", "limb","midbrain", "neural_tube","forebrain","liver","heart")
SOFTWARE = ("macs2","music","fseq","homer","hotspot","rseg","finder","bcp","mosaics","dfilter")
sh=os.system


def bamtotagalin(marks, bamfile, outfileprefix):

    if marks=='ChIP-seq_Control':
        sh('bash /home/adam/dnase/encode/Epicaller_v3/bin/bam2tagalign.part.sh {0} {1}'\
         .format(bamfile, outfileprefix))

        if outfileprefix.endswith('rep2'):
           
            sh('bash /home/adam/dnase/encode/Epicaller_v3/bin/pool.rep.part.sh {0}'\
               .format(outfileprefix))

    elif marks.startswith('H'):

        sh('bash /home/adam/dnase/encode/Epicaller_v3/bin/bam2tagalign.full.sh {0} {1}'\
          .format(bamfile, outfileprefix))
        if outfileprefix.endswith('rep2'):
             sh('bash /home/adam/dnase/encode/Epicaller_v3/bin/pool.rep.full.sh {0}'.format(outfileprefix))

    elif marks=='dnase':
        sh('bash /home/adam/dnase/encode/Epicaller_v3/bin/bam2tagalign.full.sh {0} {1}'\
         .format(bamfile, outfileprefix))
        if outfileprefix.endswith('rep2'):
            pass
            sh('bash /home/adam/dnase/encode/Epicaller_v3/bin/pool.rep.full.sh {0}'.format(outfileprefix))

    else:
        pass


def finderbam(marks, bamfile, outfileprefix):
    if outfileprefix.endswith('rep2'):
        sh('bash /home/adam/dnase/encode/Epicaller_v3/bin/finder_rep0.sh {0} {1}'\
      .format(outfileprefix, marks))


def check_folder():
    for _soft in SOFTWARE:
        _path = '/home/adam/dnase/encode/data_out/{0}'.format(_soft)
        if not os.path.isdir(_path):
            sh('mkdir {0}'.format(_path))
        for _mark in MARKS:
            _path = '/home/adam/dnase/encode/data_out/{0}/{1}'.format(_soft,_mark)
            if not os.path.isdir(_path):
                sh('mkdir {0}'.format(_path))
            for _tissue in TISSUE:
                _path = '/home/adam/dnase/encode/data_out/{0}/{1}/{2}'.format(_soft,_mark,_tissue)
                if not os.path.isdir(_path):
                    sh('mkdir {0}'.format(_path))


if __name__=='__main__':
    sh = os.system
    os.chdir('/home/adam/dnase/encode/')
    bam_table = mr('/home/adam/dnase/encode/Epicaller_v3/lib/config/bam_table.txt')
    check_folder()

    pool = Pool(processes=5)
    result = []
    for i in bam_table:
        tissue, marks, filepath, repID = i[1:]

        _ofprefix = '/home/adam/dnase/encode/data/{0}/{1}/{0}.{1}.{2}'.format(marks, tissue, repID)
        result.append(pool.apply_async(bamtotagalin, (marks, filepath, _ofprefix, )))
        finder_ofprefix = '/home/adam/dnase/encode/finder_data/{0}/{1}/{0}.{1}.{2}'.format(marks, tissue, repID)
        result.append(pool.apply_async(finderbam, (marks, filepath, finder_ofprefix, )))


    for _mark in MARKS[1:]:

        for _tissue in TISSUE:
        for _tissue in ['hindbrain']:
            sppfile='./data/{0}/{1}/{0}.{1}.rep1.filt.nodup.sample.15.tagAlign.gz.spp.ccscores'\
                .format(_mark,_tissue)

            frag=[i.rstrip().split('\t') for i in open(sppfile)]
            frag_len = frag[0][2]
            chip_prefix = '/home/adam/dnase/encode/data/{0}/{1}/{0}.{1}'.format(_mark, _tissue)
            control_prefix='/home/adam/dnase/encode/data/ChIP-seq_Control/{0}/ChIP-seq_Control.{0}'.format(_tissue)
            finder_chip_prefix='/home/adam/dnase/encode/finder_data/{0}/{1}/{0}.{1}'.format(_mark, _tissue)
            finder_control_prefix='/home/adam/dnase/encode/finder_data/ChIP-seq_Control/{0}/ChIP-seq_Control.{0}'.format(_tissue)
            result.append(pool.apply_async(Bcp.histone, (chip_prefix, control_prefix, _mark,_tissue, frag_len, )))
            result.append(pool.apply_async(Macs.histone, (chip_prefix,control_prefix,_mark,_tissue, frag_len, )))
            result.append(pool.apply_async(Hotspot.hotspot, (chip_prefix, _mark, _tissue, )))
            result.append(pool.apply_async(Dfilter.histone, (chip_prefix,control_prefix,_mark,_tissue, frag_len, )))
            result.append(pool.apply_async(Finder.histone, (finder_chip_prefix,finder_control_prefix,_mark,_tissue, )))
            result.append(pool.apply_async(Mosaics.histone, (chip_prefix,control_prefix,_mark,_tissue, )))
            result.append(pool.apply_async(Fseq.fseq, (chip_prefix,_mark,_tissue, )))
            result.append(pool.apply_async(Homer.histone, (chip_prefix,control_prefix,_mark,_tissue, frag_len, )))
            result.append(pool.apply_async(Music.histone, (chip_prefix,control_prefix,_mark,_tissue, frag_len, )))
            result.append(pool.apply_async(Rseg.histone, (chip_prefix,control_prefix,_mark,_tissue, frag_len, )))
            
    for _tissue in TISSUE[:5]:
    
        dnase_prefix = '/home/adam/dnase/encode/data/dnase/{0}/dnase.{0}'.format(_tissue)
        result.append(pool.apply_async(Macs.dnase, (dnase_prefix, _tissue, )))
        result.append(pool.apply_async(Hotspot.hotspot, (dnase_prefix, 'dnase', _tissue, )))
        result.append(pool.apply_async(Dfilter.dnase, (dnase_prefix, _tissue, )))
        result.append(pool.apply_async(Fseq.fseq, (dnase_prefix,'dnase',_tissue, )))
        result.append(pool.apply_async(Homer.dnase, (dnase_prefix,_tissue, )))
        
    pool.close()
    pool.join()

    print "DONE!"