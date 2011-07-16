import os, sys
maxgap=3
minctglen=100
f1n = sys.argv[1]
f2n = f1n+".contigs"
f1 = open(f1n,'r')
f2 = open(f2n,'w')
d1 = f1.read()
d1 = d1.split(">")
for scaf in d1:
    if len(scaf) < 2:
        continue
    ctgcnt = 1
    hdr,seq = scaf.split("\n",1)
    #only split 2 Ns or more
    seqs = seq.split(maxgap*'N')
    ctgs = []
    for s in seqs:
        if len(s.replace("\n","")) > minctglen:
            ctgs.append(s.replace("\n",""))
    
    for ctg in ctgs:
        f2.write(">%s_%d\n"%(hdr,ctgcnt))
        f2.write(ctg+"\n")
        ctgcnt +=1

f1.close()
f2.close()
