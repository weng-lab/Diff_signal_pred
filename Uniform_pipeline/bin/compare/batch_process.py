#!/usr/bin/python
import os, sys

os.chdir("/home/adam/dnase/encode/data_out")
sh = os.system

# sh("ls ./*/*/*/*.rep0.filt.nodup.pr1*.gz > ../splits/halflist1.txt")
# sh("ls ./*/*/*/*.rep0.filt.nodup.pr2*.gz > ../splits/halflist2.txt")
# sh("ls ./*/*/*/*.rep1.filt.nodup.pr1*.gz > ../splits/quarterlist1.txt")
# sh("ls ./*/*/*/*.rep1.filt.nodup.pr2*.gz > ../splits/quarterlist2.txt")
# sh("ls ./*/*/*/*.rep2.filt.nodup.pr1*.gz > ../splits/quarterlist3.txt")
# sh("ls ./*/*/*/*.rep2.filt.nodup.pr2*.gz > ../splits/quarterlist4.txt")
# sh('ls ./*/*/*/*.rep0*.gz | grep -v "filt.nodup" > ../splits/rep0_list.txt')
# sh('ls ./*/*/*/*.rep1*.gz | grep -v "filt.nodup" > ../splits/rep1_list.txt')
# sh('ls ./*/*/*/*.rep2*.gz | grep -v "filt.nodup" > ../splits/rep2_list.txt')
#
#
# ## get counts and length file of rep1 and rep2
# sh('mkdir ../rep_basic')
# rep1_file = [i.rstrip() for i in open('../splits/rep1_list.txt')]
# rep2_file = [i.rstrip() for i in open('../splits/rep2_list.txt')]
# for _file in rep1_file:
#     _mark = _file.split('/')[2]
#     _soft = _file.split('/')[1]
#     sh('bash ../summary.sh {0} {1} rep1 {2}'.format(_mark, _soft, _file))
# for _file in rep2_file:
#     _mark = _file.split('/')[2]
#     _soft = _file.split('/')[1]
#     sh('bash ../summary.sh {0} {1} rep2 {2}'.format(_mark, _soft, _file))
#
# file = [i.rstrip() for i in open('../rep_basic/counts_rep.txt')]
# file_out = [
#     i.replace('bcp', 'BCP').replace('dfilter', "DFilter").replace('finder', 'FindER').replace('fseq', 'F-seq').replace(
#         'homer', 'HOMER').replace('hotspot', 'HotSpot').replace('macs2', 'MACS2').replace('mosaics', 'MOSAiCS').replace(
#         'music', 'MUSIC').replace('rseg', 'RSEG') for i in file]
# with open('counts_rep.txt', 'w') as f:
#     for line in file_out:
#         print >> f, line
#
# file = [i.rstrip() for i in open('../rep_basic/len_rep.txt')]
# file_out = [
#     i.replace('bcp', 'BCP').replace('dfilter', "DFilter").replace('finder', 'FindER').replace('fseq', 'F-seq').replace(
#         'homer', 'HOMER').replace('hotspot', 'HotSpot').replace('macs2', 'MACS2').replace('mosaics', 'MOSAiCS').replace(
#         'music', 'MUSIC').replace('rseg', 'RSEG') for i in file]
# with open('len_rep.txt', 'w') as f:
#     for line in file_out:
#         print >> f, line

# os.chdir("../splits/")
rep0_file = [i.rstrip() for i in open('../splits/rep0_list.txt')]
rep1_file = [i.rstrip() for i in open('../splits/rep1_list.txt')]
rep2_file = [i.rstrip() for i in open('../splits/rep2_list.txt')]
rep_files = [rep1_file, rep2_file]
half_file1 = [i.rstrip() for i in open('../splits/halflist1.txt')]
half_file2 = [i.rstrip() for i in open('../splits/halflist2.txt')]
half_file = [half_file1, half_file2]
quarter_file1 = [i.rstrip() for i in open('../splits/quarterlist1.txt')]
quarter_file2 = [i.rstrip() for i in open('../splits/quarterlist2.txt')]
quarter_file3 = [i.rstrip() for i in open('../splits/quarterlist3.txt')]
quarter_file4 = [i.rstrip() for i in open('../splits/quarterlist4.txt')]
quarter_file = [quarter_file1, quarter_file2, quarter_file3, quarter_file4]
# if len(rep0_file) == len(half_file1) == len(quarter_file3):
#     print "Same files. Keep going!"
# else:
#     print "Not equal in file counts!"
#     sys.exit(1)


awk_args = "awk '{sum += $3-$2};END {print sum}'"

# sh('mkdir ../split_result/')
total = len(rep0_file)
for i in range(total):
    _soft = rep0_file[i].split('/')[1]
    _mark = rep0_file[i].split('/')[2]
    if _soft == 'macs2':

        if "broadPeak" in rep0_file[i]:
            if not os.path.isfile('../split_result/{0}.broad.{1}.split.txt'.format(_soft,_mark)):
                sh('touch ../split_result/{0}.broad.{1}.split.txt'.format(_soft,_mark))
            sh("zcat {0} | {3} >> ../split_result/{1}.broad.{2}.split.txt"\
               .format(rep0_file[i], _soft, _mark, awk_args))
            for j in range(2):
                temp_file = half_file[j][i]
                sh("intersectBed -a {0} -b {1} | uniq | {4} >> ../split_result/{2}.broad.{3}.split.txt"\
                   .format(rep0_file[i], temp_file,_soft, _mark, awk_args ))
            for k in range(4):
                temp_file = quarter_file[k][i]
                sh("intersectBed -a {0} -b {1} | uniq | {4} >> ../split_result/{2}.broad.{3}.split.txt"\
                   .format(rep0_file[i], temp_file,_soft, _mark, awk_args ))
            for l in range(2):
                temp_file = rep_files[l][i]
                sh("intersectBed -a {0} -b {1} | uniq | {4} >> ../split_result/{2}.broad.{3}.split.txt"\
                   .format(rep0_file[i], temp_file,_soft, _mark, awk_args ))

        elif "narrowPeak" in rep0_file[i]:
            if not os.path.isfile('../split_result/{0}.narrow.{1}.split.txt'.format(_soft,_mark)):
                sh('touch ../split_result/{0}.narrow.{1}.split.txt'.format(_soft,_mark))
            sh("zcat {0} | {3} >> ../split_result/{1}.narrow.{2}.split.txt"\
               .format(rep0_file[i], _soft, _mark, awk_args))
            for j in range(2):
                temp_file = half_file[j][i]
                sh("intersectBed -a {0} -b {1} | uniq | {4} >> ../split_result/{2}.narrow.{3}.split.txt"\
                   .format(rep0_file[i], temp_file,_soft, _mark, awk_args ))
            for k in range(4):
                temp_file = quarter_file[k][i]
                sh("intersectBed -a {0} -b {1} | uniq | {4} >> ../split_result/{2}.narrow.{3}.split.txt"\
                   .format(rep0_file[i], temp_file,_soft, _mark, awk_args ))
            for l in range(2):
                temp_file = rep_files[l][i]
                sh("intersectBed -a {0} -b {1} | uniq | {4} >> ../split_result/{2}.narrow.{3}.split.txt"\
                   .format(rep0_file[i], temp_file,_soft, _mark, awk_args ))

        elif "gappedPeak" in rep0_file[i]:
            if not os.path.isfile('../split_result/{0}.gapped.{1}.split.txt'.format(_soft,_mark)):
                sh('touch ../split_result/{0}.gapped.{1}.split.txt'.format(_soft,_mark))
            sh("zcat {0} | {3} >> ../split_result/{1}.gapped.{2}.split.txt"\
               .format(rep0_file[i], _soft, _mark, awk_args))
            for j in range(2):
                temp_file = half_file[j][i]
                sh("intersectBed -a {0} -b {1} | uniq | {4} >> ../split_result/{2}.gapped.{3}.split.txt"\
                   .format(rep0_file[i], temp_file,_soft, _mark, awk_args ))
            for k in range(4):
                temp_file = quarter_file[k][i]
                sh("intersectBed -a {0} -b {1} | uniq | {4} >> ../split_result/{2}.gapped.{3}.split.txt"\
                   .format(rep0_file[i], temp_file,_soft, _mark, awk_args ))
            for l in range(2):
                temp_file = rep_files[l][i]
                sh("intersectBed -a {0} -b {1} | uniq | {4} >> ../split_result/{2}.gapped.{3}.split.txt"\
                   .format(rep0_file[i], temp_file,_soft, _mark, awk_args ))
        else:
            print "Error in macs2 input file name! Exit!"
            sys.exit(1)
    elif _soft == 'music':
        if "broad" in rep0_file[i]:
            if not os.path.isfile('../split_result/{0}.broad.{1}.split.txt'.format(_soft,_mark)):
                sh('touch ../split_result/{0}.broad.{1}.split.txt'.format(_soft,_mark))
            sh("zcat {0} | {3} >> ../split_result/{1}.broad.{2}.split.txt"\
               .format(rep0_file[i], _soft, _mark, awk_args))
            for j in range(2):
                temp_file = half_file[j][i]
                sh("intersectBed -a {0} -b {1} | uniq | {4} >> ../split_result/{2}.broad.{3}.split.txt"\
                   .format(rep0_file[i], temp_file,_soft, _mark, awk_args ))
            for k in range(4):
                temp_file = quarter_file[k][i]
                sh("intersectBed -a {0} -b {1} | uniq | {4} >> ../split_result/{2}.broad.{3}.split.txt"\
                   .format(rep0_file[i], temp_file,_soft, _mark, awk_args ))
            for l in range(2):
                temp_file = rep_files[j][i]
                sh("intersectBed -a {0} -b {1} | uniq | {4} >> ../split_result/{2}.broad.{3}.split.txt"\
                   .format(rep0_file[i], temp_file,_soft, _mark, awk_args ))

        elif "punctuate" in rep0_file[i]:
            if not os.path.isfile('../split_result/{0}.punctuate.{1}.split.txt'.format(_soft,_mark)):
                sh('touch ../split_result/{0}.punctuate.{1}.split.txt'.format(_soft,_mark))
            sh("zcat {0} | {3} >> ../split_result/{1}.punctuate.{2}.split.txt"\
               .format(rep0_file[i], _soft, _mark, awk_args))
            for j in range(2):
                temp_file = half_file[j][i]
                sh("intersectBed -a {0} -b {1} | uniq | {4} >> ../split_result/{2}.punctuate.{3}.split.txt"\
                   .format(rep0_file[i], temp_file,_soft, _mark, awk_args ))
            for k in range(4):
                temp_file = quarter_file[k][i]
                sh("intersectBed -a {0} -b {1} | uniq | {4} >> ../split_result/{2}.punctuate.{3}.split.txt"\
                   .format(rep0_file[i], temp_file,_soft, _mark, awk_args ))
            for l in range(2):
                temp_file = rep_files[j][i]
                sh("intersectBed -a {0} -b {1} | uniq | {4} >> ../split_result/{2}.punctuate.{3}.split.txt"\
                   .format(rep0_file[i], temp_file,_soft, _mark, awk_args ))
        else:
            print "Error in Music!"
            sys.exit(1)


    else:
        if not os.path.isfile('../split_result/{0}.{1}.split.txt'.format(_soft,_mark)):
            sh('touch ../split_result/{0}.{1}.split.txt'.format(_soft,_mark))

        sh("zcat {0} | {3} >> ../split_result/{1}.{2}.split.txt"\
           .format(rep0_file[i], _soft, _mark, awk_args))
        for j in range(2):
            temp_file = half_file[j][i]
            sh("intersectBed -a {0} -b {1} | uniq | {4} >> ../split_result/{2}.{3}.split.txt"\
               .format(rep0_file[i], temp_file,_soft, _mark, awk_args ))
        for k in range(4):
            temp_file = quarter_file[k][i]
            sh("intersectBed -a {0} -b {1} | uniq | {4} >> ../split_result/{2}.{3}.split.txt"\
               .format(rep0_file[i], temp_file,_soft, _mark, awk_args ))
        for l in range(2):
            temp_file = rep_files[l][i]
            sh("intersectBed -a {0} -b {1} | uniq | {4} >> ../split_result/{2}.{3}.split.txt"\
               .format(rep0_file[i], temp_file,_soft, _mark, awk_args ))

peaksum_file = os.listdir('../split_result/')
os.chdir('../split_result/')
s_file=open('../consistency_summary.txt','w')
for _file in peaksum_file:
    # out_half_file = _file.rstrip('txt')+'half.txt'
    # out_quarter_file = _file.rstrip('txt')+'quarter.txt'
    _soft = _file.split('.')[0]
    _mark = _file.split('.')[1]

    temps_file = [i.rstrip() for i in open(_file)]
    temp_file = []
    for line in temps_file:
        if line=="":
            temp_file.append(int(0))
        else:
            temp_file.append(int(line))


    # h = open(out_half_file,'w')
    # q = open(out_quarter_file,'w')
    for line_idx in range(0,len(temp_file),9):

        if int(temp_file[line_idx])>0:
            half_temp = temp_file[line_idx+1:line_idx+3]
            flag = True
            for item in half_temp:
                if item >0:
                    pass
                else:
                    flag =False
            if flag:
                half_out = str(float(temp_file[line_idx+1]+float(temp_file[line_idx+2]))/(2*float(temp_file[line_idx])))
                print >>s_file , "{0}\t{1}\t{2}\thalf_split".format(half_out,_mark,_soft)

            quarter_temp = temp_file[line_idx+3:line_idx+7]
            flag = True
            for item in quarter_temp:
                if item >0:
                    pass
                else:
                    flag =False
            if flag:
                quarter_out = str(float(temp_file[line_idx+3]+float(temp_file[line_idx+4])+float(temp_file[line_idx+5])\
                                     +float(temp_file[line_idx+6]))/(4*float(temp_file[line_idx])))
                print >>s_file , "{0}\t{1}\t{2}\tquarter_split".format(quarter_out,_mark,_soft)

            rep_temp = temp_file[line_idx+7:]
            flag = True
            for item in rep_temp:
                if item >0:
                    pass
                else:
                    flag =False
            if flag:
                rep_out = str(float(temp_file[line_idx+7]+float(temp_file[line_idx+8]))/(2*float(temp_file[line_idx])))
                print >>s_file , "{0}\t{1}\t{2}\tmerge_up".format(rep_out,_mark,_soft)

s_file.close()
