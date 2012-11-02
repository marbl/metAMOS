import os, string, sys

if __name__ == "__main__":
    r = open(sys.argv[1],'r')
    seqs = {}
    prevseq = ""
    header = ""
    for line in r:
        if ">" in line:
            if prevseq != "":
                seqs[header] = prevseq
            prevseq = ""
            header = line
        else:
            prevseq += line
            
    if len(prevseq) > 0:
       seqs[header] = prevseq

    cnt = 1
    curlen = 0
    maxlen = int(sys.argv[2])
    prefix = sys.argv[1]
    if (len(sys.argv) > 3):
       prefix = sys.argv[3]
       cnt = int(sys.argv[4])

    r = open(prefix+"_part%d.fa"%(cnt),'w')            
    for id in seqs.keys():
        seq = seqs[id] 
        if curlen <= maxlen:
            r.write(id + seq)
            curlen+= len(seq)
        else:
            cnt +=1
            curlen = 0
            r.close()
            r = open(prefix+"_part%d.fa"%(cnt),'w')
            r.write(id+seq)
            curlen+= len(seq)
    r.close()
