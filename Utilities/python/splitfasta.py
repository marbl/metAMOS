import os, string, sys

if __name__ == "__main__":
    r = open(sys.argv[1],'r')
    seqs = []
    prevseq = ""
    for line in r:
        if ">" in line:
            if prevseq != "":
                seqs.append(prevseq)
            prevseq = line
        else:
            prevseq += line
            

    cnt = 1
    curlen = 0

    r = open(sys.argv[1]+"_part%d.fa"%(cnt),'w')            
    for seq in seqs:

        if curlen+len(seq) <= 6000000:
            r.write(seq)
            curlen+= len(seq)-2
        else:
            cnt +=1
            curlen = 0
            r = open(sys.argv[1]+"_part%d.fa"%(cnt),'w')
            r.write(seq)
            curlen+= len(seq)-2
    
            
            
        
    
    
