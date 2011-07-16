import sys, os, string

#os.system(" promer --maxmatch -c 4 -l 4 --coords ./CRISPR/CYANOCRISPRs.fasta")
#os.system(" promer --maxmatch -c 4 -l 4 --coords ./CRISPR/roseiCRISPRhits.fasta")

def GC(seq):
    gc =0
    gcnt = 0
    sequ = string.upper(seq)
    totc = 0
    for char in sequ:
        if char == "G" or char == "C":
            gcnt +=1
        if char != "\n" and char != ">":
            totc +=1
    return int(100*(float(gcnt)/float(totc)))
            
if __name__ == "__main__":
    if len(sys.argv) < 3:
        print "usage: annotate_contig.py <454AllContigs.fna> <prefix>"
        sys.exit(0)    
    asm_file = sys.argv[1]
    asm1 = open(asm_file,'r')
    prefix = sys.argv[2]
    ccnt = 0
    seq = ""
    asmseq = ""
    gene_dict = {}
    for line in asm1:
        asmseq += line
        if ">" not in line:
            ccnt += len(line.replace("\n",""))
            seq += line.replace("\n","")

    asm1.close()
    cov_file = open(asm_file[:-3]+"cov",'r')
    coverage = {}

    for line in cov_file:
        line = line.replace("\n","")
        data = line.split("\t")
        if data[0] == 'C':
            break
        else:
            coverage[data[1]] = int(float(data[3]))
            #reads[data[1]] = int(data[2])
            
    #run glimmer!

#    os.system("glimmer3 -i ignore.txt -g110 -t30 -o5 %s phage.icm proba"%(asm_file+".fix"))
#    os.system("glimmer3 -l  --separate_genes -g110 -t20 -o30 %s phage.icm proba"%(asm_file))
#    os.system("./g3-iterated.csh %s proba"%(asm_file+".fix"))
    os.system("./gmhmmp -o %s.orfs -m MetaGeneMark_v1.mod -a %s"%(prefix,asm_file))
    coords = open("%s.orfs"%(prefix),'r')
    coords.readline()
#    outf = open("proba.orfs",'w')
    prevhdr = 0
    curcontig = ""
    curseq = ""
    reads = {}    
    for line in coords:
        if ">gene" in line:
            if curseq != "":
                try:
                    gene_dict[curcontig].append(curseq)
                except KeyError:
                    gene_dict[curcontig] = []
                    gene_dict[curcontig].append(curseq)
                curseq = ""

            lined = line.replace("\n","")
            data = line.split(">contig",1)[1]
            curcontig = "contig"+data.split(" ")[0]
            #print curcontig
            #print int(data.split(" ")[-1].split("=",1)[-1])
            #reads[curcontig] = int(data.split(" ")[-1].split("=",1)[-1])

            prevhdr = 1
        elif len(line) > 1 and prevhdr == 1:
            curseq += line
        elif len(line) <= 1 and prevhdr == 1:
            prevhdr = 0
            
        else:
            continue
    outf = open("%s.faa"%(prefix),'w')
    
    print len(gene_dict.keys())
    orfs = {}
    for key in gene_dict.keys():
        genecnt = 1        
        for gene in gene_dict[key]:
            try:
                #print "contig"+key
                orfs["%s"%(key)] +=1
            except KeyError:
                orfs["%s"%(key)] =1                
            outf.write(">%s_gene%d\n%s"%(key,genecnt,gene))
            genecnt +=1
#        print gene_dict[key][0]
    outf.close()
    #os.system(" phmmer --cpu 10 -E 0.1 -o %s.phm.out --tblout %s.phm.tbl --notextw %s.faa ./DBS/allprots.faa"%(prefix,prefix,prefix))
    hit_dict = {}
    phmout = open("%s.phm.tbl"%(prefix),'r')
    phmmer_hits = {}
    ctghits = {}
    annot = {}
    for line in phmout:
        line = line.replace("\n","")
        
        if "contig" in line:
            line,phage_annot = line.split("[",1)
            phage_annot = phage_annot.replace("]","")
            data = line.split(" ")
            data2 = []
            for item in data:
                if item == "" or item == "-" or item == "\n":
                    continue
                else:
                    data2.append(item)
            try:
                data2[16]        
                for git in data2[16:]:
                    phage_annot += " "+git + " "
                
            except IndexError:
                pass

            data2 = data2[:15]
            #print phage_annot
            #print data2
            #print data2[1].split("_",1)[0]
            try:
                ctghits[data2[1]]
                continue
            except KeyError:
                ctghits[data2[1]] = 1
                pass
            phage_annot = phage_annot.replace(",","")
            try:            
                annot[data2[1].split("_",1)[0]] += phage_annot
            except KeyError:
                annot[data2[1].split("_",1)[0]] = phage_annot
            try:
                phmmer_hits[data2[1].split("_",1)[0]] +=1
            except KeyError:
                phmmer_hits[data2[1].split("_",1)[0]] = 1                
            try:
                hit_dict[data2[1]]
            except KeyError:
                hit_dict[data2[1]] = [float(data2[2]),int(float(data2[3])),phage_annot]
    print len(hit_dict.keys())
    #for key in hit_dict.keys():
    #    print hit_dict[key]

    trna = {}
    os.system("tRNAscan-SE -o trna.out -B %s > /dev/null"%(asm_file))
    trnaf = open("trna.out",'r')
    for line in trnaf:
        if "contig" in line:
            data = line.split("\t")
            ckey = data[0].strip()
            #print len(ckey)
            trna[ckey] = int(data[1])
            #sys.exit(0)
    os.system("rm trna.out")
    overlaps = {}
    os.system("python getContigRepeats.py %s %s.overlaps"%(asm_file,prefix))
    overlapf = open("%s.overlaps"%(prefix),'r')
    for line in overlapf:
        line = line.replace("\n","")
        data = line.split(",")
        ctgids = []
        for item in data:
            if item != "":
                
                ctgids.append(item)
        ctgids.sort()
        #width = 5
        #zeros = 5-len(str(ctgids[0]))
#        ckey = "contig" + '0'*(5-len(str(ctgids[0])))+str(ctgids[0])
        ckey = ctgids[0]
        #print ckey
        #sys.exit(1)
        try:
            overlaps[ckey]
            for id in ctgids[1:]:
                if id not in overlaps[ckey]:
                    overlaps[ckey] += id+";"
        except KeyError:
            overlaps[ckey] = ""            
            for id in ctgids[1:]:
                overlaps[ckey] += id+";"
                
    
    data = asmseq.split(">")[1:]
    cyano_hits = {}
    rosei_hits = {}
    cyano = {}
    rosei = {}
    fhdr = ""
    repeats = {}
    contigs = {}
    for item in data:
        hdr,seq = item.split("\n",1)
        hdr = hdr.replace("\n","").replace(">","")
        fhdr = hdr.split(" ")[0]
        print fhdr
        #print int(hdr.split(" ")[-1].split("=",1)[-1])
        reads[fhdr] = int(hdr.split(" ")[-1].split("=",1)[-1])
        contigs[fhdr] = seq
        ff = open("t1.fna",'w')
        ff.write(">proba\n%s"%(seq))
        ff.close()
        os.system("./repeatoire --minreplen=20 --z=11 --extend=0 --allow-redundant=0 --sequence=t1.fna --output=reps.out >& test.out")
        repfile = open("reps.out",'r')
        repd = repfile.read()
        repeats[fhdr] = repd.count("Alignment #")
        
        os.system(" promer --maxmatch -c 4 -l 4 --coords t1.fna ./CRISPR/CYANOCRISPRs.fasta > t1.out 2> t2.err")
        os.system("show-coords -I 70 -o -k -L 10 -T -c -q out.delta > out.coords")
        ffi = open("out.coords",'r')

        for line in ffi:
            line = line.replace("\n","")
            
            if "proba" in line:
                data = line.split("\t")
                pid = int(float(data[6]))
                cov = int(float(data[10]))
                crisp = data[-1]
                if pid < 80 or cov < 80:
                    continue
                try:
                    cyano[fhdr] +=1
                except KeyError:
                    cyano[fhdr] = 1                
                try:
                    cyano_hits[fhdr].append([pid,cov,crisp])
                except KeyError:
                    cyano_hits[fhdr] = []
                    cyano_hits[fhdr].append([pid,cov,crisp])                    

        os.system(" promer --maxmatch -c 4 -l 4 --coords t1.fna ./CRISPR/roseiCRISPRhits.fasta > t1.out 2> t2.err")
        os.system("show-coords -I 70 -o -k -L 12 -T -c -r out.delta > out.coords")
        ffi = open("out.coords",'r')        
        for line in ffi:
            line = line.replace("\n","")
            
            if "proba" in line:
                data = line.split("\t")
                pid = int(float(data[6]))
                cov = int(float(data[10]))
                crisp = data[-1]
                if pid < 80 or cov < 80:
                    continue
                try:
                    rosei[fhdr] +=1
                except KeyError:
                    rosei[fhdr] = 1
                try:
                    rosei_hits[fhdr].append([pid,cov,crisp])
                except KeyError:
                    rosei_hits[fhdr] = []
                    rosei_hits[fhdr].append([pid,cov,crisp])                    


        #sys.exit(1)        
    #os.system(" promer --maxmatch -c 4 -l 4 --coords ./CRISPR/CYANOCRISPRs.fasta")
    #os.system(" promer --maxmatch -c 4 -l 4 --coords ./CRISPR/roseiCRISPRhits.fasta")
    report_file = "%s.report"%(prefix)
    report = open(report_file,'w')
    report.write("asm,contig,reads,len,cov,gc,orfs,phmmerhits,cyanoCR,roseiCR,overlaps,repeats,tRNA,annot\n")
    ckeys = contigs.keys()
    ckeys.sort()
    for contig in ckeys:
        ctg = contigs[contig]
        ctg_id = contig.replace("\n","").replace(">","").split(" ")[0]
        #print ctg_id
        try:
            reads[ctg_id]
        except KeyError:
            reads[ctg_id] = 0

        try:
            coverage[ctg_id]
        except KeyError:
            coverage[ctg_id] = 0

        try:
            orfs[ctg_id]
        except KeyError:
            orfs[ctg_id] = 0

        try:
            phmmer_hits[ctg_id]
        except KeyError:
            phmmer_hits[ctg_id] = 0

        try:
            rosei[ctg_id]
        except KeyError:
            rosei[ctg_id] = 0

        try:
            cyano[ctg_id]
        except KeyError:
            cyano[ctg_id] = 0

        try:
            annot[ctg_id]
        except KeyError:
            annot[ctg_id] = "NA"                                                


        try:
            overlaps[ctg_id]
        except KeyError:
            overlaps[ctg_id] = "NA"

        try:
            repeats[ctg_id]
        except KeyError:
            repeats[ctg_id] = 0

        try:
            trna[ctg_id]
        except KeyError:
            trna[ctg_id] = 0                                                            
            
          
        try:
            report.write("%s,%s,%d,%d,%d,%d,%d,%d,%d,%d,%s,%d,%d,%s\n"%(prefix,contig,reads[ctg_id],len(ctg),coverage[ctg_id],GC(ctg),orfs[ctg_id],phmmer_hits[ctg_id],cyano[ctg_id],rosei[ctg_id],overlaps[ctg_id],repeats[ctg_id],trna[ctg_id],annot[ctg_id]))
        except KeyError:
            print ctg_id
    report.close()
    #for each contig:
    #get ORFs
    #    blastx/phmmer
    #promer CRISPRs
    #nucmer CRISPRs
    #record coverage (lookup in Newbler output)
    #record GC content
    #run repeatoire, list other contigs with shared repeats?
    #make assemblies/contigs available for Michelle
    #make phage DB available
    
