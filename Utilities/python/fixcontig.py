import os, sys, string

if __name__ == "__main__":
    r = open("./out/SRS023915.denovo_duplicates_marked.trimmed.contig",'r')
    of = open("./out/s12.contig",'w')
    incontig = 0
    consensus = ""
    for line in r:
        if "##" in line:
            line = line.replace("bases","bases,")
            line = line.replace("contig_","")
            incontig = 1
            consensus = ""
            of.write(line)
        elif "#" in line and incontig:
            consensus = consensus.replace("\n","")
            width = 70
            i = 0
            newcon = ""
            while i + width < len(consensus):
                newcon += consensus[i:i+width]+"\n"
                i+=width
            newcon += consensus[i:]+"\n"
            of.write(newcon)
            incontig = 0
            of.write(line)
            consensus = ""
        elif "#" in line:
            of.write(line)
        else:
            #should be consensus
            consensus += line
    r.close()
    of.close()
