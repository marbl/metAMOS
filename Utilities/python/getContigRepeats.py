import sys, string, os


def concatContig(ctgfile):
    if len(sys.argv) < 3:
        print "usage: contig_file out_file" 
    contig_file = open(ctgfile,'r')
    out_file = open(ctgfile+".merged",'w')
    out_data = ""
    for line in contig_file.xreadlines():
        if ">" not in line:
             out_data += line.replace("\n","")
    width = 60
    pp = 0
    out_file.write(">seq\n")
    while pp+60 < len(out_data):
        out_file.write(out_data[pp:pp+60]+"\n")
        pp +=60
    out_file.write(out_data[pp:]+"\n")
    out_file.close()
    contig_file.close()



if __name__ == "__main__":
    contig_repeats = ""
    contig_file = ""
    if len(sys.argv) < 3:
        print "usage: getContigRepeats.py <contig_file> <out_file>"
        sys.exit(1)
#    contig_repeats = open("myreps.out",'w')
    try:
        contig_repeats = open(sys.argv[2],'w')
    except IOError, errmsg:
        print "Error creating output file %s "%(sys.argv[2]), errmsg
        sys.exit(1)

    try:
        contig_file = open(sys.argv[1],'r')
    except IOError, errmsg:
        print "Error openinig input file %s "%(sys.argv[1]), errmsg
        sys.exit(1)
    contig_file.close()
    contig_file = open(sys.argv[1],'r')
    concatContig(sys.argv[1])

    if 1:
        os.system("/fs/szdevel/metAMOS/Utilities/cpp/repeatoire --minreplen=200 --z=17 --sequence=%s.merged --xmfa=%s.xmfa"%(sys.argv[1],sys.argv[1]))
#    repeat_file = open(sys.argv[2],'r')
    repeat_file = open(sys.argv[1]+".xmfa",'r')
    ctg_dict = {}
    seq_map = {}
    contig_data = contig_file.read()
    num_contigs = contig_data.count(">")
    contig_data = contig_data.split(">")[1:]
    prev_pos = 0
    eid = ""
    iid = 1
    for contig in contig_data:
        hdr,seq = contig.split("\n",1)
        id = hdr.split(" ")[0]
        hdr = hdr.replace(">","").replace("\n","")
        start = prev_pos
        clen = len(seq.replace("\n",""))
        end = prev_pos+clen
        ctg_dict[iid] = [start, end, seq]
        i = 0
        while i < clen:
            seq_map[prev_pos+i] = hdr#iid
            i+=1
        prev_pos = end+1
        iid +=1

    repeat_data = repeat_file.readlines()
    repfam = 1
    reppos = []
    clc = 1
    for line in repeat_data:
        if "=" in line:
          repfam +=1
          ctg_list = []
          for copy in reppos:
             try:
                #print seq_map[int(copy[0])]
                if seq_map[int(copy[0])] == seq_map[int(copy[1])]:
                    ctg_list.append(seq_map[int(copy[0])])
                    #ctg_list.append(seq_map[copy[1]])
             except KeyError:
                 continue
          #print ctg_list

          if len(ctg_list) > 1 and ctg_list.count(ctg_list[0]) != len(ctg_list):
              for item in ctg_list:
                   contig_repeats.write("%d:"%repfam+str(item)+"\n")
          clc +=1
          reppos = []
        if ">" not in line:
            continue
        gg, info = line.split(":",1)
        spos,info = info.split("-",1)
        epos,info = info.split(" ",1)
        orient, info = info.split(" ",1)
#        print spos, epos, orient
        reppos.append([spos,epos])
        
