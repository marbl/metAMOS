import sys, string
#3 formats, 33, 56, 59?
#determine hdr string
#
if __name__ == "__main__":
    if len(sys.argv) < 3:
        print "usage: out_prefix input_file"
        sys.exit(1)
    seq = sys.argv[1]
    #"abaylyi_sep"
    dd_file = sys.argv[2]
    try:
        dd = open (dd_file,'r')
    except IOError:
        print dd, " does not exist!"
        sys.exit(1)
    of = open("./out/%s_reads.fasta"%(seq),'w')
    oq = open("./out/%s_reads.qual"%(seq),'w')    
    olist = []
    olist2 = []
    ddi = dd.readlines()
    pline1 = 0
    pline2 = 0
    cnt = 1
    seqid = ""

    for line in ddi:

        if "@" in line:
                pline1 =1
                of.write(line.replace("@",">"))

        elif pline1 == 1:
                of.write(line)
                pline1 = 0
        elif "+" in line:
                pline2 = 1
                oq.write(line.replace("+",">"))
                cnt +=1
        elif pline2 == 1:
                pline2 = 0
                qualstr = ""
                for char in line[:-1]:
                    qualstr += "%d "%(ord(char)-62)
                oq.write(qualstr+"\n")
        else:
                pline1 = 0
                pline2 = 0

    of.close()
    oq.close()

