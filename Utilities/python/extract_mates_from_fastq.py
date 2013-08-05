import string, sys

#if __name__ == "__main__":
def extract_mates_from_fastq(infile):
    f1 = open(infile,'r')
    f2 = open("%s.mates"%(infile),'w')
    while 1:
        s1 = f1.readline()
        s2 = f1.readline()
        s3 = f1.readline()
        s4 = f1.readline()

        l1 = f1.readline()
        l2 = f1.readline()
        l3 = f1.readline()
        l4 = f1.readline()
        if l4 == "":
            break
        m1 = s1.split(" ")[0]
        m2 = l1.split(" ")[0]
        m1 = m1.replace("\n","")
        m2 = m2.replace("\n","")

        f2.write(m1+"\t"+m2+"\n")
        f2.flush()
    f1.close()
    f2.close()
