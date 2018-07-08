import os,sys

infile=sys.argv[1]
outfile=sys.argv[2]
status=sys.argv[3]
with open(outfile, "w") as f:
    myf=[i.rstrip().split("\t") for i in open(infile)]
    mydict={}
    for line in myf:
        keys="\t".join(line[:3])
        value=line[3]
        if not mydict.has_key(keys):
            mydict[keys]=value
        else:
            pass
    for _key in mydict.keys():
        print >>f, "{0}\t{1}".format(mydict[_key],str(status))
