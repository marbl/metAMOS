#!python

import os, sys, math, string, time, BaseHTTPServer, getopt, re, subprocess, webbrowser, shelve
from operator import itemgetter
from datetime import datetime

from utils import *
from preprocess import Preprocess
from assemble import Assemble
sys.path.append(INITIAL_UTILS)
from ruffus import *
import pysam

MAX_ENTRIES = 10000

_readlibs = []
_skipsteps = []
_settings = Settings()
_asm = None
_mapper = "bowtie"
_ctgbpcov = False
_lowmem= False
def init(reads, skipsteps, asm,mapper,savebtidx,ctgbpcov,lowmem):
   global _readlibs
   global _asm
   global _skipsteps
   global _savebtidx
   global _ctgbpcov
   global _mapper
   global _lowmem
   _mapper = mapper
   _readlibs = reads
   _skipsteps = skipsteps
   _asm = asm
   _savebtidx = savebtidx
   _ctgbpcov = ctgbpcov
   _lowmem = lowmem
def meanstdv(x):
    n, mean, std = len(x), 0, 0
    for a in x:
       mean = mean + a
    mean = mean / float(n)
    for a in x:
        std = std + (a - mean)**2
    std = math.sqrt(std / float(n-1))
    return mean, std

def map2contig():
    bowtie_mapping = 1
    
    readDir = ""
    asmDir = ""
    #threads = 0
    #run_process(_settings, "cp %s/Assemble/out/%s.asm.contig2 %s/Assemble/out/%s.asm.contig"%(_settings.rundir,_settings.PREFIX,_settings.rundir,_settings.PREFIX))
    tigr_file = open("%s/Assemble/out/%s.asm.tigr"%(_settings.rundir,_settings.PREFIX),'w')
    contigfile = open("%s/Assemble/out/%s.asm.contig"%(_settings.rundir,_settings.PREFIX),'r')

    seqdict = {}
    hdr = ""
    cnt = 0

    contigdict = {}
    dictcnt = 0
    if _lowmem:
        contigdict = shelve.open("%s/Assemble/out/%s.contigdict"%(_settings.rundir,_settings.PREFIX),'c', writeback=True)
    contigdict2 = {}
    readdict = {}
    matedict = {}
    ctgmates = 0
    matectgdict = {}
    mateotdict = {}
    read_lookup = {}
    readcnt = 1
    mapped_reads = {}
    readcontig_dict  = {}
    strand_dict = {}
    fiveprimeend_dict = {}
    for lib in _readlibs:
         

        matefile = open("%s/Preprocess/out/lib%d.seq.mates"%(_settings.rundir,lib.id),'r')
        matedict[lib.id] = {}
        for line in matefile.xreadlines():
            line = line.replace("\n","")
            mate1, mate2 = line.split("\t")
            mate1 = mate1.replace("@","").replace(">","")
            mate2 = mate2.replace("@","").replace(">","")
            matedict[lib.id][mate2] = mate1
            #matedict[lib.id][mate1] = mate2
            read_lookup[readcnt] = mate1
            read_lookup[readcnt+1] = mate2
            readcnt += 2

    if bowtie_mapping == 1:
        if _mapper == "bowtie":
           if "bowtie" not in _skipsteps and not _savebtidx:# and not os.path.exists("%s/Assemble/out/IDX.1.ebwt"%(_settings.rundir)):
               run_process(_settings, "%s/bowtie-build -o 2 %s/Assemble/out/%s.asm.contig %s/Assemble/out/IDX"%(_settings.BOWTIE, _settings.rundir,_settings.PREFIX,_settings.rundir),"MapReads")
        elif _mapper == "bowtie2":
           if not os.path.exists("%s/Assemble/out/IDX.1.bt2"%(_settings.rundir)):
               run_process(_settings, "%s/bowtie2-build -o 2 %s/Assemble/out/%s.asm.contig %s/Assemble/out/IDX"%(_settings.BOWTIE, _settings.rundir,_settings.PREFIX,_settings.rundir),"MapReads")
        for lib in _readlibs:
            seqfile = open("%s/Preprocess/out/lib%d.seq.btfilt"%(_settings.rundir,lib.id),'w')


            #trim to 25bp
            trim = 0
            if trim:
                f1 = open("%s/Preprocess/out/lib%d.seq"%(_settings.rundir,lib.id))
                f2 = open("%s/Preprocess/out/lib%d.seq.trim"%(_settings.rundir,lib.id),'w')
                linecnt = 1
                for line in f1.xreadlines():
                    if linecnt % 2 == 0:
                        f2.write(line[0:25]+"\n")
                    else:
                        f2.write(line)
                    linecnt +=1
                f1.close()
                f2.close()
            if _mapper == "bowtie":
                if "bowtie" not in _skipsteps and (lib.format == "fasta" or lib.format == "sff"):
                    if trim:
                        run_process(_settings, "%s/bowtie -p %d -f -v 1 -M 2 --un %s/Assemble/out/lib%d.unaligned.fasta %s/Assemble/out/IDX %s/Preprocess/out/lib%d.seq.trim > %s/Assemble/out/lib%d.bout"%(_settings.BOWTIE,_settings.threads,_settings.rundir,lib.id_settings.rundir,_settings.rundir,lib.id,_settings.rundir,lib.id),"MapReads")
                    else:
                        run_process(_settings, "%s/bowtie -p %d -f -l 25 -e 140 --best --strata -m 10 -k 1 --un %s/Assemble/out/lib%d.unaligned.fasta %s/Assemble/out/IDX %s/Preprocess/out/lib%d.seq > %s/Assemble/out/lib%d.bout"%(_settings.BOWTIE,_settings.threads,_settings.rundir,lib.id,_settings.rundir,_settings.rundir,lib.id,_settings.rundir,lib.id),"MapReads")
                elif "bowtie" not in _skipsteps and lib.format != "fasta":
                    if trim:
                        run_process(_settings, "%s/bowtie  -p %d -v 1 -M 2 --un %s/Assemble/out/lib%d.unaligned.seq %s/Assemble/out/IDX %s/Preprocess/out/lib%d.seq.trim > %s/Assemble/out/lib%d.bout"%(_settings.BOWTIE,_settings.threads,_settings.rundir,lib.id,_settings.rundir,_settings.rundir,lib.id,_settings.rundir,lib.id),"MapReads")
                    else:
                        run_process(_settings, "%s/bowtie  -p %d -l 25 -e 140 --best --strata -m 10 -k 1 --un %s/Assemble/out/lib%d.unaligned.seq %s/Assemble/out/IDX %s/Preprocess/out/lib%d.seq > %s/Assemble/out/lib%d.bout"%(_settings.BOWTIE,_settings.threads,_settings.rundir,lib.id,_settings.rundir,_settings.rundir,lib.id,_settings.rundir,lib.id),"MapReads")
                    run_process(_settings, "java -cp %s convertFastqToFasta %s/Assemble/out/lib%d.unaligned.seq %s/Assemble/out/lib%d.unaligned.fasta %s/Assemble/out/lib%d.unaligned.fasta.qual"%(_settings.METAMOS_JAVA, _settings.rundir, lib.id, _settings.rundir, lib.id, _settings.rundir, lib.id), "MapReads")

                infile = open("%s/Assemble/out/lib%d.bout"%(_settings.rundir,lib.id),'r')
                readctgfile =  open("%s/Assemble/out/%s.lib%dcontig.reads"%(_settings.rundir,_settings.PREFIX, lib.id),'w')

                for line1 in infile.xreadlines():
                    line1 = line1.replace("\n","")
                    ldata = line1.split("\t")
                    if "Warning" in line1 or "warning" in line1:
                        continue 
                    if len(ldata) < 6:
                        continue
                    read = ldata[0]
                    strand = ldata[1]
                    contig = ldata[2]
                    spos = ldata[3]
                    try:
                        int(spos)
                    except ValueError:
                        #bowtie output for this line malformed, skip
                        continue 
                    read_seq = ldata[4]
                    read_qual = ldata[5]
                    read = read.split(" ")[0]
                    try:
                        epos = int(spos)+len(read_seq)
                    except ValueError:
                        continue
                    mapped_reads[read] = 1
                    strand_dict[read] = strand
                    readcontig_dict[read] = contig
                    readctgfile.write("%s\t%s\n"%(read, contig))
                    if strand == "+":
                        fiveprimeend_dict[read] = int(spos)
                    else:
                        fiveprimeend_dict[read] = int(epos)
                    try:
                        #if lowmem:
                        contigdict[contig].append([int(spos), int(epos), strand, read,len(read_seq),lib.id])
                        #contigdict_disk[contig].append([int(spos), int(epos), strand, read,len(read_seq),lib.id])
                    except KeyError:
                        contigdict[contig] = [[int(spos),int(epos),strand,read,len(read_seq),lib.id]]
                        #contigdict_disk[contig] = [[int(spos),int(epos),strand,read,len(read_seq),lib.id]]
                    if _lowmem and dictcnt % MAX_ENTRIES == 0:
                       #print "Processed %d entries %s"%(dictcnt, str(datetime.now()))
                       contigdict.sync()
                    dictcnt+=1
                    #print contig
                    seqdict[read] = read_seq
                    seqfile.write(">%s\n%s\n"%(read,read_seq))
                    seqfile.flush()
                readctgfile.close()
            elif _mapper == "bowtie2":
                if "bowtie" not in _skipsteps and (lib.format == "fasta" or lib.format == "sff"):
                    run_process(_settings, "%s/bowtie2 -p %d -f -D 15 -R 2 -N 0 -L 20 -i S,1,1.10 --un %s/Assemble/out/lib%d.unaligned.fasta %s/Assemble/out/IDX %s/Preprocess/out/lib%d.seq -S %s/Assemble/out/lib%d.sam"%(_settings.BOWTIE2,_settings.threads,_settings.rundir,lib.id,_settings.rundir,_settings.rundir,lib.id,_settings.rundir,lib.id),"MapReads")
                elif "bowtie" not in _skipsteps and lib.format != "fasta":
                    run_process(_settings, "%s/bowtie2 -p %d -D 15 -R 2 -N 0 -L 20 -i S,1,1.10 --un %s/Assemble/out/lib%d.unaligned.seq %s/Assemble/out/IDX %s/Preprocess/out/lib%d.seq -S %s/Assemble/out/lib%d.sam"%(_settings.BOWTIE2,_settings.threads,_settings.rundir,lib.id,_settings.rundir,_settings.rundir,lib.id,_settings.rundir,lib.id),"MapReads") 
                    run_process(_settings, "java -cp %s convertFastqToFasta %s/Assemble/out/lib%d.unaligned.seq %s/Assemble/out/lib%d.unaligned.fasta %s/Assemble/out/lib%d.unaligned.fasta.qual"%(_settings.METAMOS_JAVA, _settings.rundir, lib.id, _settings.rundir, lib.id, _settings.rundir, lib.id), "MapReads") 

                if not os.path.exists("%s/Assemble/out/lib%d.sam"%(_settings.rundir,lib.id)):
                    pass
                else:

                   samfile = pysam.Samfile("%s/Assemble/out/lib%d.sam"%(_settings.rundir,lib.id),'r')
                   readctgfile =  open("%s/Assemble/out/%s.lib%dcontig.reads"%(_settings.rundir,_settings.PREFIX, lib.id),'w') 

                   allreads = samfile.fetch()
                   #print len(allreads)
                   cnt = 1
                   for read in allreads:
                       if read == None:
                           continue

                       #print cnt
                       #print read
                       cnt +=1
                       read_seq = read.seq
                       read_qual = read.qual
                       readname = read.qname
                       #print readname
                       epos = read.pos+read.rlen
                       spos = read.pos
                       contig = ""
                       if read.tid >= 0:
                           contig = samfile.getrname(read.tid)
                           readctgfile.write("%s\t%s\n"%(readname, contig))
                       else:
                           continue
                       strand = "+"
                       if read.is_reverse:
                           strand = "-"
                       mapped_reads[readname] = 1
                       strand_dict[readname] = strand
                       readcontig_dict[readname] = contig
                       if strand == "+":
                           fiveprimeend_dict[readname] = int(spos)
                       else:
                           fiveprimeend_dict[readname] = int(epos)
                       try:
                           contigdict[contig].append([int(spos), int(epos), strand, readname,len(read_seq),lib.id])
                           #contigdict_disk[contig].append([int(spos), int(epos), strand, readname,len(read_seq),lib.id])
                       except KeyError:
                           contigdict[contig] = [[int(spos),int(epos),strand,readname,len(read_seq),lib.id]]
                           #contigdict_disk[contig] = [[int(spos),int(epos),strand,readname,len(read_seq),lib.id]]
                       if _lowmem and dictcnt % MAX_ENTRIES == 0:
                          #print "Processed %d entries"%(dictcnt, str(datetime.now()))
                          contigdict.sync()
                       dictcnt+=1
                       seqdict[readname] = read_seq
                       seqfile.write(">%s\n%s\n"%(readname,read_seq))
                       seqfile.flush()
                       print readname, "end"
                   print "close sam"
                   samfile.close()
                   readctgfile.close()
                   #pileup_columns = samfile.pileup( refg,spos,epos)
                   #for pc in pileup_columns :
                   #    perbp_cov.write('%s,%s\n' % (pc.pos, pc.n))
                   #alignedreads = samfile.fetch(refg,spos,epos)



    else:
        if 0:
 
            #open soap ReadOnContig
            #some contigs are missing!
            infile = open("%s/Assemble/out/%s.asm.readOnContig"%(_settings.rundir,_settings.PREFIX),'r')
            #readID, ContigID, startpos, strand
            hdr = infile.readline()
            linecnt = 1
            for line in infile.xreadlines():
                if linecnt % 100000 == 0:
                    #print linecnt,
                    sys.stdout.flush()
                data = line.replace("\n","").split("\t")
                #print data
                if len(data) < 4:
                    continue
                contig = data[1]
                spos = int(data[2])
                if spos < 0:
                    spos = 0
                epos = spos+readlen
                strand = data[3]
                read = int(data[0])

                try:
                    contigdict[contig].append([int(spos), int(spos)+epos, strand, read_lookup[read]])
                    #contigdict_disk[contig].append([int(spos), int(spos)+epos, strand, read_lookup[read]])
                except KeyError:
                    contigdict[contig] = [[int(spos),int(spos)+epos,strand,read_lookup[read]]]
                    #contigdict_disk[contig] = [[int(spos),int(spos)+epos,strand,read_lookup[read]]]
                read_seq = "TEST"
            
                seqdict[read_lookup[read]] = read_seq
                linecnt +=1
        
    contig_data = contigfile.read()
    contig_data = contig_data.split(">")
    errfile = open("%s/Assemble/out/contigs_wo_location_info.txt"%(_settings.rundir),'w')
    new_ctgfile = open("%s/Assemble/out/%s.seq100.contig"%(_settings.rundir,_settings.PREFIX),'w')
    ctgcnt = 1
    ctgseq = 0
    ctgsizes = []
    n50_size = 0
    n50_mid = 955,000
    ctg_cvg_file = open("%s/Assemble/out/%s.contig.cvg"%(_settings.rundir,_settings.PREFIX),'w')
    ctg_cnt_file = open("%s/Assemble/out/%s.contig.cnt"%(_settings.rundir,_settings.PREFIX),'w')
    libcov_dict = {}
    if _ctgbpcov:
        for lib in _readlibs:
            libcov_dict["lib%d"%(lib.id)] = {}
    for item in contig_data:
        if item == '':
            continue

        item = item.split("\n",1)
        ref = item[0].split(" ")[0]
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
        ctgslen = len(cseq)
        #libcov_dict[ref] = {}
        if _ctgbpcov:
            for lib in _readlibs:
                ii = 0
                while ii < ctgslen:
                    try:
                        libcov_dict["lib%d"%(lib.id)][ref][ii] = 0
                    except KeyError:
                        libcov_dict["lib%d"%(lib.id)][ref] = {}
                        libcov_dict["lib%d"%(lib.id)][ref][ii] = 0
                    ii+=1
        try:
            tigr_file.write("##%s %d %d bases, 00000000 checksum.\n"%(ref.replace(">",""),len(contigdict[ref]), len(item[1])))
            tigr_file.flush()
        except KeyError:
            #print "oops, not in mapping file\n"
            errfile.write("%s\n"%ref)
            continue
        new_ctgfile.write(">%d\n%s"%(ctgcnt,cseq_fmt))#item[1]))
        ctgcnt +=1
        tigr_file.write(cseq_fmt)#item[1])
        contigdict[ref].sort()
        #if lowmem
        #contigdict_disk[ref].sort()
        #print contigdict[ref]
        try:
            ctg_cvg_file.write("%s\t%.2f\n"%(ref,(float(len(contigdict[ref])*len(seqdict[contigdict[ref][0][3]]))/float(ctgslen))))
            ctg_cnt_file.write("%s\t%d\n"%(ref,len(contigdict[ref])))
        except KeyError:
            #no reads map to this contig, skip?
            continue
        
        for read in contigdict[ref]:
            #libcovfile = open("%s/Assemble/out/%s.%s.contig.cov"%(_settings.rundir,_settings.PREFIX,read[3][0:4]),'w')
            try:
                #if read[0] <= 500 and ctgslen - (int(read[1])) <= 500:
                matectgdict[read[-2]] = ref
                mateotdict[read[-2]] = read[2]
            except KeyError:
                pass
            ii = 0
            #contigdict2[ref] = item[1]
            if _ctgbpcov:
                while ii < read[-2]:
                    try:     
                        libcov_dict["lib%d"%(read[-1])][ref][read[0]+ii]+=1
                    except KeyError:
                        libcov_dict["lib%d"%(read[-1])][ref][read[0]+ii] = 1
                    ii+=1
            if read[2] == "-":
                tigr_file.write("#%s(%d) [RC] %d bases, 00000000 checksum. {%d 1} <%d %d>\n"%(read[3],read[0]-1, read[-2], read[-2], read[0], read[1]))
            else:
                tigr_file.write("#%s(%d) [] %d bases, 00000000 checksum. {1 %d} <%d %d>\n"%(read[3],read[0]-1, read[-2], read[-2], read[0], read[1]))
            tigr_file.write(seqdict[read[3]]+"\n")

   

    for lib in _readlibs:
       
        libid = "lib%d"%(lib.id)
        if _ctgbpcov:#for libid in libcov_dict.keys():
            libcovfile = open("%s/Assemble/out/lib%d.contig.cov"%(_settings.rundir,lib.id),'w')
    
            for ctgid in libcov_dict[libid].keys():
                libcovfile.write(">%s\n"%(ctgid))
                for pos in libcov_dict[libid][ctgid].keys():
                    libcovfile.write("%d,%d\n"%(pos,libcov_dict[libid][ctgid][pos]))
            libcovfile.close()
        mateheader = open("%s/Assemble/out/%s.lib%d.hdr"%(_settings.rundir,_settings.PREFIX,lib.id),'w')
        new_matefile = open("%s/Assemble/out/%s.lib%d.mappedmates"%(_settings.rundir,_settings.PREFIX,lib.id),'w')
        badmatefile = open("%s/Assemble/out/%s.lib%d.badmates"%(_settings.rundir,_settings.PREFIX,lib.id),'w')
        ctgmatefile = open("%s/Assemble/out/%s.lib%d.mates_in_diff_contigs"%(_settings.rundir,_settings.PREFIX,lib.id),'w')

        #    for lib in _readlibs:
        linked_contigs = {}
        insertlens = []
        oldstdev = 0
        oldstdev = (lib.mmin+lib.mmax)/6
        oldmean = (lib.mmin+lib.mmax)/2
        oldmax = oldmean+oldstdev
        oldmin = oldmean-oldstdev
        if oldmin < 0:
           oldmin = 0
        for mate in matedict[lib.id].keys():
            matepair = matedict[lib.id][mate]
            try:
                mapped_reads[mate]
                mapped_reads[matepair]
            except KeyError:
                continue

            if readcontig_dict[mate] != readcontig_dict[matepair]:
                new_matefile.write("%s\t%s\t%d\n"%(mate,matepair,lib.id))
                ctgmatefile.write("%s\t%s\t%d\n"%(mate,matepair,lib.id))
            elif strand_dict[mate] == "+" and strand_dict[matepair] == "-":
                new_matefile.write("%s\t%s\t%d\n"%(mate,matepair,lib.id))
                ilen = fiveprimeend_dict[matepair]-fiveprimeend_dict[mate]
                if ilen < oldmax and ilen > oldmin:
                    insertlens.append(ilen)
            elif strand_dict[mate] == "-" and strand_dict[matepair] == "+":
                new_matefile.write("%s\t%s\t%d\n"%(matepair,mate,lib.id))
                ilen = fiveprimeend_dict[matepair]-fiveprimeend_dict[mate]
                if ilen < oldmax and ilen > oldmin:
                    insertlens.append(ilen)

            elif strand_dict[mate] == "+" and strand_dict[matepair] == "+":
                #output to file
                badmatefile.write("%s\t%s\t%d\n"%(mate,matepair,lib.id))
            elif strand_dict[mate] == "-" and strand_dict[matepair] == "-":
                #output to file
                badmatefile.write("%s\t%s\t%d\n"%(mate,matepair,lib.id))
            new_matefile.flush()
            badmatefile.flush()
            ctgmatefile.flush()
            continue

        if len(insertlens) > 0:
           lmin = min(insertlens)
           lmax = max(insertlens)
           lavg = sum(insertlens)/len(insertlens)
           lmean,lstdev = meanstdv(insertlens)
           lib.mmin = lmin
           lib.mmax = lmax
           lib.mean = lmean
           lib.stdev = lstdev
           #if lavg * 1.2 < lmax or lavg * 0.8 > lmin:
           #    lmin = 
           
           #if verbose, output
           #print "Old insert length min: ", lib.mmin
           #print "New insert length min: ", lmin
           #print "Old insert length max: ", lib.mmax
           #print "New insert length max: ", lmax
        else:
           lmin = lib.mmin
           lmax = lib.mmax
        mateheader.write("library\t%d\t%d\t%d\n"%(lib.id,lmin,lmax))
        new_matefile.close()
        badmatefile.close()
        mateheader.close()
        run_process(_settings, "cat %s/Assemble/out/%s.lib%d.mappedmates >> %s/Assemble/out/%s.lib%d.hdr "%(_settings.rundir,_settings.PREFIX, lib.id,_settings.rundir,_settings.PREFIX,lib.id),"MapReads")
        run_process(_settings, "cp %s/Assemble/out/%s.lib%d.hdr %s/Assemble/out/%s.lib%d.mappedmates "%(_settings.rundir,_settings.PREFIX, lib.id,_settings.rundir,_settings.PREFIX,lib.id),"MapReads")
    ctg_cvg_file.close()
    ctg_cnt_file.close()
    tigr_file.close()
    if _lowmem:
        contigdict.close()
    else:
        contigdict = {}

@files("%s/Assemble/out/%s.asm.contig"%(_settings.rundir,_settings.PREFIX),"%s/Assemble/out/mapreads.success"%(_settings.rundir))
#@posttask(create_symlink,touch_file("completed.flag"))
@posttask(touch_file("%s/Logs/mapreads.ok"%(_settings.rundir)))
@follows(Assemble)
def MapReads(input,output):

   if "MapReads" in _skipsteps or "mapreads" in _skipsteps:
      run_process(_settings, "touch %s/Logs/mapreads.skip"%(_settings.rundir), "MapReads")
      return 0
   if _mapper == "bowtie" or _mapper == "bowtie2":
       map2contig()
   else:
       print "Read mapper not supported, time to exit"
       sys.exit(1)
   run_process(_settings, "touch %s/Assemble/out/mapreads.success"%(_settings.rundir), "MapReads")
   #stop here, for now
   #sys.exit(0)
   #check if sucessfully completed   
