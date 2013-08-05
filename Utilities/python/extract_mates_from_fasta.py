import string, sys

#if __name__ == "__main__":
def extract_mates_from_fasta(infile):
    f1 = open(infile,'r')
    f2 = open("%s.mates"%(infile),'w')
    first = 1
    second = 0
    firstmate = ""
    linecnt = 0
    for line in f1.xreadlines():
        #if linecnt % 2 == 0:#">" not in line:
        if ">" in line:

            line = line.replace(">","")
            line = line.replace("\n","")
            data = line.split(" ")
            mate= data[0]
            mate = mate.strip()
        else:
            linecnt +=1
            continue
        if first:
            firstmate = mate
            first = 0        
            second = 1
        elif second:
            f2.write(firstmate+"\t"+mate+"\n")
            f2.flush()
            first = 1
            second = 0
        else:
            linecnt+=1
            continue
        linecnt +=1
    f1.close()
    f2.close()
