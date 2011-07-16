import string
from operator import itemgetter

def map2contig(self):
    bowtie_mapping = 1
    
    #infile = open("sample2.seq",'r')
    #locfile = open("sample2.loc2",'r')
    #contigfile = open("sample2.contigs",'r')
    #matefile = open("sample2.mates2",'r')
    readDir = ""
    asmDir = ""
    threads = 0
    #System.err.println("Launching bowtie aligner: perl /fs/wrenhomes/sergek/Utils/get_singles.pl -reads " + readDir + " -assembly " + asmDir + " --threads 2");    
system "$bowtie-build $asmdir/contig.fa $resultDir/IDX";
system("$bowtie -v 1 -M 2 $resultDir/IDX $resultDir/$prefix.trim.fastq >& $resultDir/$prefix.bout");
    infile = open("seq100.bowtie",'r')
    contigfile = open("seq100.fna",'r')
    seqfile = open("seq100.seq",'w')
    matefile = open("allmates.txt",'r')
    tigr_file = open("seq100.tigr",'w')
    new_matefile = open("mappedmates_500.txt",'w')
    seqdict = {}
    hdr = ""
    cnt = 0
    matedict = {}
    contigdict = {}
    contigdict2 = {}
    readdict = {}
    matedict = {}
    ctgmates = 0
    matectgdict = {}
    mateotdict = {}
    for line in matefile.xreadlines():
        line = line.replace("\n","")
        mate1, mate2 = line.split("\t")
        matedict[mate2] = mate1


            
        #matedict[mate2] = mate1
        #print mate1, mate2

    if bowtie_mapping:
        for line1 in infile.xreadlines():
            line1 = line1.replace("\n","")
            ldata = line1.split("\t")
            if len(ldata) < 6:
                continue
            read, strand, contig, spos,read_seq, read_qual  = ldata[:6]
            epos = int(spos)+len(read_seq)
            try:
                contigdict[contig].append([int(spos), int(spos)+epos, strand, read])
            except KeyError:
                contigdict[contig] = [[int(spos),int(spos)+epos,strand,read]]
            seqdict[read] = read_seq
            seqfile.write(">%s\n%s\n"%(read,read_seq))
            
    else:

        for line in locfile.xreadlines():
            line = line.replace("\n","")
            read, contig, spos, epos, strand = line.split("\t")
            try:
                contigdict[contig].append([int(spos), epos, strand, read])
            except KeyError:
                contigdict[contig] = [[int(spos),epos,strand,read]]
            
        for line1 in infile.xreadlines():
            line1 = line1.replace("\n","")
            if line1 == "":
                cnt = 0
                continue
            if cnt == 0:
                hdr = line1.replace(">","")
                cnt = 1
            else:
                seqdict[hdr] = line1
                cnt  = 0

    contig_data = contigfile.read()
    contig_data = contig_data.split(">")
    errfile = open("contigs_wo_location_info.txt",'w')
    new_ctgfile = open("seq100.contig",'w')
    ctgcnt = 1
    ctgseq = 0
    ctgsizes = []
    n50_size = 0
    n50_mid = 955,000
    for item in contig_data:
        if item == '':
            continue

        item = item.split("\n",1)
        ref = item[0]
        ref = ref.replace("\n","")
        cseq = item[1].replace("\n","")
        ctgseq+=len(cseq)
        ctgsizes.append(len(cseq))
        i = 0
        cpos = 0
        width = 70
        cseq_fmt = ""
        while i+width < len(cseq):
            cseq_fmt += cseq[i:i+width]+"\n"
            i+= width
        cseq_fmt += cseq[i:]+"\n"
        ctgslen = len(item[1])
        #contigdict2[ref] = item[1]
        try:
            tigr_file.write("##%d %d %d bases, 00000000 checksum.\n"%(ctgcnt,len(contigdict[ref]), len(item[1])))
        except KeyError:
            print "oops, not in mapping file\n"
            errfile.write("%s\n"%ref)
            continue
        new_ctgfile.write(">%d\n%s"%(ctgcnt,cseq_fmt))#item[1]))
        ctgcnt +=1
        tigr_file.write(cseq_fmt)#item[1])
        contigdict[ref].sort()
        #print contigdict[ref]
        for read in contigdict[ref]:
            
            try:
                if read[0] <= 500 and ctgslen - (int(read[1])) <= 500:
                    matectgdict[read[-1]] = ref
                    mateotdict[read[-1]] = read[2]
            except KeyError:
                pass
            if read[2] == "-":
#                matedict[read[-1]]
                tigr_file.write("#%s(%d) [RC] %d bases, 00000000 checksum. {%d 1} <%d %s>\n"%(read[-1],read[0]-1, len(seqdict[read[-1]]), len(seqdict[read[-1]]), read[0], read[1]))
            else:
                tigr_file.write("#%s(%d) [] %d bases, 00000000 checksum. {1 %d} <%d %s>\n"%(read[-1],read[0]-1, len(seqdict[read[-1]]), len(seqdict[read[-1]]), read[0], read[1]))
            tigr_file.write(seqdict[read[-1]]+"\n")

   
    nonctgmates = 0
    sffcnt = 0
    sfrcnt = 0
    srfcnt = 0
    srrcnt = 0
    nffcnt = 0
    nfrcnt = 0
    nrfcnt = 0
    nrrcnt = 0
    new_matefile.write("library\t110110\t240.0\t300.00\n")
    linked_contigs = {}
    for mate in matedict.keys():
        try:
            matectgdict[mate]
            matectgdict[matedict[mate]]
            mateotdict[mate]

        except KeyError:
            continue
        new_matefile.write("%s\t%s\t110110\n"%(mate,matedict[mate]))
        if matectgdict[mate] == matectgdict[matedict[mate]]:
            ctgmates +=1
            if mateotdict[mate] == "+" and mateotdict[matedict[mate]] == "+":
                sffcnt +=1
            elif mateotdict[mate] == "+" and mateotdict[matedict[mate]] == "-":
                sfrcnt +=1
            elif mateotdict[mate] == "-" and mateotdict[matedict[mate]] == "-":
                srrcnt +=1
            elif mateotdict[mate] == "-" and mateotdict[matedict[mate]] == "+":
                srfcnt +=1
        else:
            nonctgmates +=1
            ctgid = [matectgdict[mate],matectgdict[matedict[mate]]]
            ctgid.sort()
            try:
                linked_contigs[ctgid[0]+"*"+ctgid[1]] += 1
            except KeyError:
                linked_contigs[ctgid[0]+"*"+ctgid[1]] = 1

            if mateotdict[mate] == "+" and mateotdict[matedict[mate]] == "+":
                nffcnt +=1
            elif mateotdict[mate] == "+" and mateotdict[matedict[mate]] == "-":
                nfrcnt +=1
            elif mateotdict[mate] == "-" and mateotdict[matedict[mate]] == "-":
                nrrcnt +=1
            elif mateotdict[mate] == "-" and mateotdict[matedict[mate]] == "+":
                nrfcnt +=1

    ffcnt = 0
    frcnt = 0
    rfcnt = 0
    rrcnt = 0


    
    for mate in matedict.keys():
        try:
            matectgdict[mate]
            matectgdict[matedict[mate]]
        except KeyError:
            continue
        if mateotdict[mate] == "+" and mateotdict[matedict[mate]] == "+":
            ffcnt +=1
        elif mateotdict[mate] == "+" and mateotdict[matedict[mate]] == "-":
            frcnt +=1
        elif mateotdict[mate] == "-" and mateotdict[matedict[mate]] == "-":
            rrcnt +=1
        elif mateotdict[mate] == "-" and mateotdict[matedict[mate]] == "+":
            rfcnt +=1

    #ctgcnt = 1
    #ctgseq = 0
    #ctgsizes = []
    #n50_size = 0
    #n50_mid = 955,000
    print "total contigs: ", ctgcnt
    print "max contig size: ", max(ctgsizes)
    print "min contig size: ", min(ctgsizes)
    print "average contig size: ",float(sum(ctgsizes))/float(ctgcnt)
    i = 0
    cp = 0
    ctgsizes.sort()
    ctgsizes = ctgsizes[::-1]
    while i < len(ctgsizes):
        cp += ctgsizes[i]
        if cp >= 955000:
            n50_size = ctgsizes[i]
            break
        i+=1
    print "N50 size: ", n50_size
    print "total pairs of contigs with links:", len(linked_contigs.keys())
    items = linked_contigs.items()
    items.sort(key=itemgetter(1),reverse=True)
    print items[0:10]
    print "---->    ---->", ffcnt
    print "---->    <----", frcnt
    print "<----    ---->", rfcnt
    print "<----    <----", rrcnt
    print "mates located in same contig: ", ctgmates
    print "---->    ---->", sffcnt
    print "---->    <----", sfrcnt
    print "<----    ---->", srfcnt
    print "<----    <----", srrcnt
    print "mates located in diff contig: ", nonctgmates
    print "---->    ---->", nffcnt
    print "---->    <----", nfrcnt
    print "<----    ---->", nrfcnt
    print "<----    <----", nrrcnt

    
    tigr_file.close()
        
    
        
        

    
    
