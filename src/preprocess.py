#!python

import os, sys, string, time, BaseHTTPServer, getopt, re, subprocess, webbrowser
from operator import itemgetter

from utils import *

sys.path.append(INITIAL_UTILS)
from ruffus import *

from extract_mates_from_fasta import *
from extract_mates_from_fastq import *

_filter = None
_readlibs = []
_asmcontigs = []
_skipsteps = []
_asm = None
_run_fastqc = False
_settings = Settings() 

def init(reads, asmcontigs, skipsteps, asm, run_fastqc,filter):
   global _readlibs
   global _asmcontigs
   global _skipsteps
   global _asm
   global _run_fastqc
   global _filter
   _readlibs = reads
   _asmcontigs = asmcontigs
   _skipsteps = skipsteps
   _asm = asm
   _run_fastqc = run_fastqc
   _filter = filter

def LCS(S1, S2):
    M = [[0]*(1+len(S2)) for i in xrange(1+len(S1))]
    longest, x_longest = 0, 0
    for x in xrange(1,1+len(S1)):
        for y in xrange(1,1+len(S2)):
            if S1[x-1] == S2[y-1]:
                M[x][y] = M[x-1][y-1] + 1
                if M[x][y]>longest:
                    longest = M[x][y]
                    x_longest  = x
            else:
                M[x][y] = 0
    return S1[x_longest-longest: x_longest]

def convertFastaToFastq(fastaFile, qualFile, outputFile, stepName):
   run_process(_settings, "java -cp %s convertFastaAndQualToFastq %s %s > %s"%(_settings.METAMOS_JAVA, fastaFile, qualFile, outputFile), stepName)

def convertInputFastaToFastq(libID, mated):
   run_process(_settings, "ln %s/Preprocess/out/lib%d.seq %s/Preprocess/out/lib%d.fasta"%(_settings.rundir, libID, _settings.rundir, libID), "Preprocess")
   convertFastaToFastq("%s/Preprocess/out/lib%d.seq"%(_settings.rundir, libID), "%s/Preprocess/out/lib%d.seq.qual"%(_settings.rundir, libID), "%s/Preprocess/out/lib%d.fastq"%(_settings.rundir, libID), "Preprocess")
   if mated:
      convertFastaToFastq("%s/Preprocess/out/lib%d.1.fasta"%(_settings.rundir, libID), "%s/Preprocess/out/lib%d.1.fasta.qual"%(_settings.rundir, libID), "%s/Preprocess/out/lib%d.1.fastq"%(_settings.rundir, libID), "Preprocess")
      convertFastaToFastq("%s/Preprocess/out/lib%d.2.fasta"%(_settings.rundir, libID), "%s/Preprocess/out/lib%d.2.fasta.qual"%(_settings.rundir, libID), "%s/Preprocess/out/lib%d.2.fastq"%(_settings.rundir, libID), "Preprocess")

def convertFastqToFasta(fastqFile, fastaFile, qualFile, stepName):
   run_process(_settings, "java -cp %s convertFastqToFasta %s %s %s"%(_settings.METAMOS_JAVA, fastqFile, fastaFile, qualFile), stepName)

def convertInputFastqToFasta(libID, mated):
   run_process(_settings, "ln %s/Preprocess/out/lib%d.seq %s/Preprocess/out/lib%d.fastq"%(_settings.rundir, libID, _settings.rundir, libID), "Preprocess")
   convertFastqToFasta("%s/Preprocess/out/lib%d.seq"%(_settings.rundir, libID), "%s/Preprocess/out/lib%d.fasta"%(_settings.rundir, libID), "%s/Preprocess/out/lib%d.fasta.qual"%(_settings.rundir, libID), "Preprocess")
   if mated:
      run_process(_settings, "java -cp %s convertFastqToFasta %s/Preprocess/out/lib%d.1.fastq %s/Preprocess/out/lib%d.1.fasta %s/Preprocess/out/lib%d.1.fasta.qual"%(_settings.METAMOS_JAVA, _settings.rundir, libID, _settings.rundir, libID, _settings.rundir, libID), "Preprocess")
      run_process(_settings, "java -cp %s convertFastqToFasta %s/Preprocess/out/lib%d.2.fastq %s/Preprocess/out/lib%d.2.fasta %s/Preprocess/out/lib%d.2.fasta.qual"%(_settings.METAMOS_JAVA, _settings.rundir, libID, _settings.rundir, libID, _settings.rundir, libID), "Preprocess")

def parseInterleaved(rf,wf,fastq=True):
    if 1:
        if 1:
            if 1:
                   #this means we have this entire lib in one file
                   #parse out paired record (8 lines), rename header to be filename + "/1" or "/2", and remove reads with N
                   rf = open(read.path,'r')
                   wf = open(read.path.replace("/in/","/out/"),'w')
                   start = 1
                   rcnt = 0
                   recordcnt = 0
                   record = []
                   shdr = ""
                   for line in rf.xreadlines():
                       if start:
                           s1hdr = line
                           record.append(line)
                           start = 0
                           rcnt =1

                       else:
                           if rcnt == 7:
                               #end of record
                               record.append(line)
                               rcnt +=1
                               if len(record) != 8:
                                   #something went wrong
                                   continue
                               rseq = string.upper(record[0]+record[5])                               
                               if "N" in rseq:
                                   #skip both, dont' want Ns
                                   continue
                               #update hdrs to be filename /1 or /2
                               recordcnt +=1
                               hdr = read.sid+"r"+str(recordcnt)+"/"
                               #hdr2 = lib.sid[0:3]+str(int(lib.sid[3:])+1)+str(recordcnt)+"/"
                               if fastq == True:
                                   wf.writelines("@"+hdr+"1\n")
                                   wf.writelines(record[1])
                                   wf.writelines("+"+hdr+"1\n")
                                   wf.writelines(record[3])
                                   wf.writelines("@"+hdr+"2\n")
                                   wf.writelines(record[5])
                                   wf.writelines("+"+hdr+"2\n")
                                   wf.writelines(record[7])
                               else:
                                   wf.writelines(">"+hdr+"1\n")
                                   wf.writelines(record[1])
                                   wf.writelines(">"+hdr+"1\n")
                                   wf.writelines(record[3])
                                   wf.writelines(">"+hdr+"2\n")
                                   wf.writelines(record[5])
                                   wf.writelines(">"+hdr+"2\n")
                                   wf.writelines(record[7])
                           elif rcnt % 4 == 0:
                               s2hdr = line
                               rlcs = LCS(s1hdr,s2hdr)
                               #these should almost identical
                               if float(len(rlcs))/float(len(s1hdr)) < 0.9:
                                   #missing record somewhere, start over with this one
                                   s1hdr = line
                                   record = [line]
                                   start = 0
                                   rcnt = 1
                               else:
                                   record.append(line)
                                   rcnt +=1
                           elif rcnt % 2 == 0:
                               #quality hdr
                               record .append(line)
                               rcnt +=1
                           else:
                               record .append(line)
                               rcnt +=1
                   #update to new path
                   read.path = read.path.replace("/in/","/out/")            

#@transform(readpaths,["%s/Preprocess/out/all.seq"%(_settings.rundir),"%s/Preprocess/out/all.seq.mates"%(_settings.rundir)])
@posttask(touch_file("%s/Logs/preprocess.ok"%(_settings.rundir)))
@files(_settings.readpaths,"%s/Preprocess/out/preprocess.success"%(_settings.rundir))
#filtreadpaths)
def Preprocess(input,output):
   global _run_fastqc

   for contig in _asmcontigs:
      run_process(_settings, "ln %s/Preprocess/in/%s %s/Preprocess/out/%s.asm.contig"%(_settings.rundir, contig, _settings.rundir, os.path.splitext(contig)[0]), "Preprocess")

   # update file names if necessary to avoid conflicts and create qual files
   for lib in _readlibs:
      for read in lib.reads:
         if "lib%d"%(lib.id) in os.path.basename(read.path):
            if lib.mated and not lib.interleaved:
                readpair = lib.getPair(read.id)
                if readpair == -1:
                    #not interleaved and mated, yet do not have 2nd file..
                    continue
                (nprefix, nmatch, nsuffix) = readpair.path.rpartition("lib%d"%(lib.id))
                if len(nmatch) == 0:
                    print "Error: Could not find expected input file %s\n"%(readpair.path)
                    raise(JobSignalledBreak)
                npath = nprefix + "inputLib%d"%(lib.id) + nsuffix
                run_process(_settings, "cp %s %s"%(readpair.path, npath), "Preprocess")
                readpair.path = npath
                readpair.fname = os.path.basename(readpair.path)

            (nprefix, nmatch, nsuffix) = read.path.rpartition("lib%d"%(lib.id))
            if (len(nmatch) == 0):
                print "Error: Could not find expected input file %s\n"%(readpair.path)
                raise(JobSignalledBreak)
            npath = nprefix + "inputLib%d"%(lib.id) + nsuffix
            run_process(_settings, "cp %s %s"%(read.path, npath), "Preprocess")
            read.path = npath
            read.fname = os.path.basename(read.path)

         if lib.format == "fasta" and not os.path.isfile("%s/Preprocess/in/%s.qual"%(_settings.rundir, read.fname)):
            run_process(_settings, "java -cp %s:. outputDefaultQuality %s/Preprocess/in/%s > %s/Preprocess/in/%s.qual"%(_settings.METAMOS_JAVA, _settings.rundir, read.fname, _settings.rundir, read.fname), "Preprocess")
            read.qformat = "sanger"
            if lib.mated and not lib.interleaved:
                readpair = lib.getPair(read.id)
                if readpair == -1:
                    #not interleaved and mated, yet do not have 2nd file..
                    continue
                run_process(_settings, "java -cp %s:. outputDefaultQuality %s/Preprocess/in/%s > %s/Preprocess/in/%s.qual"%(_settings.METAMOS_JAVA, _settings.rundir, readpair.fname, _settings.rundir, readpair.fname), "Preprocess")
                readpair.qformat = "sanger"

   #move input files into Preprocess ./in dir
   #output will either be split fastq files in out, or AMOS bank
   if "Preprocess" in _skipsteps or "preprocess" in _skipsteps:
       for lib in _readlibs:
           for read in lib.reads:
               run_process(_settings, "ln %s/Preprocess/in/%s %s/Preprocess/out/"%(_settings.rundir,read.fname,_settings.rundir),"Preprocess")
       return 0
   run_process(_settings, "rm %s/Preprocess/out/all.seq.mates"%(_settings.rundir), "Preprocess")

   global _filter
   if _filter == "eautils" and not os.path.exists(_settings.EAUTILS + "fastq-mcf"):
      _filter = "metamos"

   if _filter == "metamos":
       #print "filtering.."
     
       #for reads+libs
       cnt = 1
       for lib in _readlibs:
           for read in lib.reads:
               if not read.filtered and read.format == "fastq" and read.mated and read.interleaved:
                   #this means we have this entire lib in one file
                   #parse out paired record (8 lines), rename header to be filename + "/1" or "/2", and remove reads with N
                   rf = open(read.path,'r')
                   npath = read.path.replace("/in/","/out/")
                   #readpath,base = os.path.split(npath)
                   #newpath = readpath+"lib%d"%(lib.id)
                   wf = open(npath,'w')
                   #wf = open(read.path.replace("/in/","/out/"),'w')
                   #wf = open(readpath+"lib%d"%(lib.id),'w')
                   start = 1
                   rcnt = 0
                   recordcnt = 0
                   record = []
                   shdr = ""
                   for line in rf.xreadlines():
                       if start:
                           s1hdr = line
                           record.append(line)
                           start = 0
                           rcnt =1

                       else:
                           if rcnt == 7:
                               #end of record
                               record.append(line)
                               rcnt +=1
                               if len(record) != 8:
                                   #something went wrong
                                   continue
                               rseq = string.upper(record[0]+record[5])                               
                               if "N" in rseq:
                                   #skip both, dont' want Ns
                                   continue
                               #update hdrs to be filename /1 or /2
                               recordcnt +=1
                               hdr = read.sid+"r"+str(recordcnt)+"/"
                               #hdr2 = lib.sid[0:3]+str((int(lib.sid[3:])+1))+str(recordcnt)+"/"
                               wf.writelines("@"+hdr+"1\n")
                               wf.writelines(record[1])
                               wf.writelines("+"+hdr+"1\n")
                               wf.writelines(record[3])
                               wf.writelines("@"+hdr+"2\n")
                               wf.writelines(record[5])
                               wf.writelines("+"+hdr+"2\n")
                               wf.writelines(record[7])
                           elif rcnt % 4 == 0:
                               s2hdr = line
                               rlcs = LCS(s1hdr,s2hdr)
                               #these should almost identical
                               if len(rlcs)+2 != len(s1hdr) and float(len(rlcs))/float(len(s1hdr)) < 0.9:
                                   #missing record somewhere, start over with this one
                                   s1hdr = line
                                   record = [line]
                                   start = 0
                                   rcnt = 1
                               else:
                                   record.append(line)
                                   rcnt +=1
                           elif rcnt % 2 == 0:
                               #quality hdr
                               record .append(line)
                               rcnt +=1
                           else:
                               record .append(line)
                               rcnt +=1
                   #update to new path
                   read.path = read.path.replace("/in/","/out/")            
                   #read.fname = "lib%d"%(lib.id)
                   read.filtered = True
                   wf.close()
               elif not read.filtered and read.format == "fastq" and read.mated and not read.interleaved:
                   readpair = lib.getPair(read.id)
                   if readpair == -1:
                       #not interleaved and mated, yet do not have 2nd file..
                       continue
                   rf1 = open(read.path,'r')
                   wf1 = open(read.path.replace("/in/","/out/"),'w')
                   rf2 = open(readpair.path,'r')
                   wf2 = open(readpair.path.replace("/in/","/out/"),'w')
                   recordcnt = 0
                   f1cnt = 0
                   f2cnt = 0
                   while 1:                   
                       rs1 = rf1.readline()
                       rs2 = rf1.readline()
                       rs3 = rf1.readline()
                       rs4 = rf1.readline()

                       rp1 = rf2.readline()
                       rp2 = rf2.readline()
                       rp3 = rf2.readline()
                       rp4 = rf2.readline()

                       if rs1 == "" or rs2 == "" or rs3 == "" or rs4 == "":
                           #EOF or something went wrong, break
                           break 
                       rseq = string.upper(rs2)                               
                       if "N" in rseq:
                           continue
                       if rp1 == "" or rp2 == "" or rp3 == "" or rp4 == "":
                           #EOF or something went wrong, break
                           break 
                       f1cnt +=4
                       f2cnt +=4
                       rseq = string.upper(rp2)                               
                       if "N" in rseq:
                           continue

                       #rlcs = LCS(rs1,rp1)
                       if 0:#float(len(rlcs))/float(len(rs1)) < 0.9:
                           #not aligned! something needs to be removed?
                           #go on to next
                           continue
                       else:
                           #record.append(line)
                           #rcnt +=1
                           recordcnt +=1
                           hdr = read.sid+"r"+str(recordcnt)+"/"
                           wf1.writelines("@"+hdr+"1\n")
                           wf1.writelines(rs2)
                           wf1.writelines("+"+hdr+"1\n")
                           wf1.writelines(rs4)
                           hdr2 = readpair.sid[0:3]+str(int(lib.sid[3:])+1)+str(recordcnt)+"/"
                           wf2.writelines("@"+hdr+"2\n")
                           wf2.writelines(rp2)
                           wf2.writelines("+"+hdr+"2\n")
                           wf2.writelines(rp4)
                           wf1.flush()
                           wf2.flush()
                   if f1cnt != f2cnt:
                       print "Error: error in library, read files not of equal length!"
                       raise(JobSignalledBreak)
                   readpair.filtered = True
                   read.filtered = True
                   read.path = read.path.replace("/in/","/out/")
                   readpair.path = readpair.path.replace("/in/","/out/")
                   wf1.close()
                   wf2.close()
               elif not read.filtered and read.format == "fastq" and not read.mated:
                   #this is easy, just throw out reads with Ns
                   rf = open(read.path,'r')
                   wf = open(read.path.replace("/in/","/out/"),'w')
                   while 1:
                       rs1 = rf.readline()
                       rs2 = rf.readline()
                       rs3 = rf.readline()
                       rs4 = rf.readline()
                       if rs1 == "" or rs2 == "" or rs3 == "" or rs4 == "":
                           #EOF or something went wrong, break
                           break 
                       rseq = string.upper(rs2)                               
                       if "N" in rseq:
                           continue
                       wf.writelines(rs1)
                       wf.writelines(rs2)
                       wf.writelines(rs3)
                       wf.writelines(rs4)
                   read.path = read.path.replace("/in/","/out/")
                   read.filtered = True
                   wf.close()
               elif not read.filtered and read.format == "fasta" and read.mated and read.interleaved:
                   #this means we have this entire lib in one file
                   #parse out paired record (4 lines), rename header to be filename + "/1" or "/2", and remove reads with N
                   rf = open(read.path,'r')
                   rq = open(read.path+".qual", 'r')
                   npath = read.path.replace("/in/","/out/")
                   #print npath
                   #readpath,base = os.path.split(npath)
                   #newpath = readpath+"lib%d"%(lib.id)
                   wf = open(npath,'w')
                   wq = open(npath+".qual", 'w')
                   #wf = open(read.path.replace("/in/","/out/"),'w')
                   start = 1
                   rcnt = 0
                   recordcnt = 0
                   record = []
                   shdr = ""
                   reads = rf.read().split(">")[1:]
                   quals = rq.read().split(">")[1:]
                   if len(reads) % 2 != 0:
                       print "Read file corrupted, please fix and restart!"
                       sys.exit(1)

                   prevok = False
                   first = True
                   second = False
                   prevseq = ""
                   prevqual = ""
                   readcnt = 1
                   currIndex = 0
                   for currIndex, rd in enumerate(reads):
                       if first:
                           hdr,seq = rd.split("\n",1)
                           if "N" in string.upper(seq) or len(seq) < 2:
                               prevok = False
                           else: 
                               prevok = True
                               prevseq = seq
                               try:
                                   prevqual = quals[currIndex].split("\n",1)[1]
                               except IndexError:
                                   prevqual = ["40\n"]

                           second = True
                           first = False
                       elif second:
                           hdr,seq = rd.split("\n",1)
                           if "N" in string.upper(seq) or len(seq) < 2:
                               pass
                           elif prevok:
                               hdr = read.sid+"r"+str(readcnt)+"/"
                               wf.writelines(">"+hdr+"1\n")
                               wf.writelines(prevseq)
                               wf.writelines(">"+hdr+"2\n")
                               wf.writelines(seq)

                               wq.writelines(">"+hdr+"1\n")
                               wq.writelines(prevqual)
                               wq.writelines(">"+hdr+"2\n")
                               try:
                                   wq.writelines(quals[currIndex].split("\n",1)[1])
                               except IndexError:
                                   wq.writelines(["40\n"])

                               readcnt +=1
                           second = False
                           first = True

                   #update to new path
                   read.path = read.path.replace("/in/","/out/")            
                   read.filtered = True
                   wf.close()
                   wq.close()
               elif not read.filtered and read.format == "fasta" and read.mated and not read.interleaved:
                   readpair = lib.getPair(read.id)
                   if readpair == -1:
                       #not interleaved and mated, yet do not have 2nd file..
                       continue
                   rf1 = open(read.path,'r')
                   qf1 = open(read.path+".qual",'r')
                   wf1 = open(read.path.replace("/in/","/out/"),'w')
                   wq1 = open(read.path.replace("/in/","/out/")+".qual", 'w')
                   rf2 = open(readpair.path,'r')
                   qf2 = open(readpair.path+".qual", 'r')
                   wf2 = open(readpair.path.replace("/in/","/out/"),'w')
                   wq2 = open(readpair.path.replace("/in/","/out/")+".qual", 'w')

                   recordcnt = 0
                   while 1:                   
                       rs1 = rf1.readline()
                       rs2 = ""
                       for line in rf1:
                          if ">" in line:
                             break
                          rs2 += line.rstrip()
                       qs1 = qf1.readline()
                       qs2 = ""
                       for line in qf1:
                          if ">" in line:
                             break;
                          qs2.append(line.rstrip()) 

                       if rs1 == "" or rs2 == "":
                           #EOF or something went wrong, break
                           break 
                       rseq = string.upper(rs2)                               
                       if "N" in rseq:
                           continue
                       rp1 = rf2.readline()
                       rp2 = ""
                       for line in rf2:
                          if ">" in line:
                             break
                          rp2.append(line.rstrip()) 
                       qp1 = qf2.readline()
                       qp2 = ""
                       for line in qf2:
                          if ">" in line:
                             break;
                          qs2.append(line.rstrip())

                       if rp1 == "" or rp2 == "":
                           #EOF or something went wrong, break
                           break 
                       rseq = string.upper(rp2)                               
                       if "N" in rseq:
                           continue

                       rlcs = LCS(rs1,rp1)
                       if 0:#float(len(rlcs))/float(len(rs1)) < 0.9:
                           #not aligned! something needs to be removed?
                           #go on to next
                           continue
                       else:
                           #record.append(line)
                           #rcnt +=1
                           recordcnt +=1
                           hdr = read.sid+"r"+str(recordcnt)+"/"
                           wf1.writelines(">"+hdr+"1\n")
                           wf1.writelines(rs2)
                           wq1.writelines(">"+hdr+"1\n")
                           wq2.writelines(qs2)
                           wf2.writelines(">"+hdr+"2\n")
                           wf2.writelines(rp2)
                           wq2.writelines(">"+hdr+"2\n")
                           wq2.writelines(qp2)

                   readpair.filtered = True
                   read.filtered = True
                   read.path = read.path.replace("/in/","/out/")
                   readpair.path = readpair.path.replace("/in/","/out/")
                   wf1.close()
                   wf2.close()
                   wq1.close()
                   wq2.close()
               elif not read.filtered and read.format == "fasta" and not read.mated:
                   #easiest case, check for Ns
                   rf = open(read.path,'r')
                   rq = open(read.path+".qual", 'r')
                   wf = open(read.path.replace("/in/","/out/"),'w')
                   wq = open(read.path.replace("/in/","/out/")+".qual", 'w')

                   reads = rf.read().split(">")[1:]
                   quals = rq.read().split(">")[1:]
                   if len(reads) % 2 != 0:
                       print "Read file corrupted, please fix and restart!"
                       sys.exit(1)

                   readcnt = 1
                   for currIndex, rd in enumerate(reads):
                      hdr,seq = rd.split("\n",1)
                      if "N" in seq:
                         continue
                      hdr = read.sid+"r"+str(readcnt)
                      wf.writelines(">"+hdr+"\n")
                      wf.writelines(seq)
                      wq.writelines(">"+hdr+"\n")
                      wq.writelines(quals[currIndex].split("\n",1)[1])
                      readcnt += 1
                   read.path = read.path.replace("/in/","/out/")
                   read.filtered = True
                   wf.close()
                   wq.close()
           cnt +=1
   elif _filter == "eautils":
      offset = ""
      commandOptions = ""
      for lib in _readlibs:
         firstName = ""
         secondName = ""
         for read in lib.reads:
            if offset == "":
               offset = read.qformat
            elif not offset == read.qformat:
               print "Error: inconsistent PHRED offsets in libraries. Previous library had %s and current library %d has %s\n"%(offset, lib.id, read.qformat)
               raise(JobSignalledBreak)

         if lib.format == "fasta":
            # convert to fastq uninterleaved
            if lib.mated:
               if lib.interleaved:
                  firstName = "%s.1.fastq"%(os.path.splitext(lib.f1.fname)[0])
                  secondName = "%s.2.fastq"%(os.path.splitext(lib.f1.fname)[0])
                  convertFastaToFastq("%s/Preprocess/in/%s"%(_settings.rundir, lib.f1.fname), "%s/Preprocess/in/%s.qual"%(_settings.rundir, lib.f1.fname), "%s/Preprocess/in/%s.fastq"%(_settings.rundir, os.path.splitext(lib.f1.fname)[0]), "Preprocess")
                  run_process(_settings, "perl %s/perl/split_fastq.pl %s/Preprocess/in/%s.fastq %s/Preprocess/in/%s %s/Preprocess/in/%s"%(_settings.METAMOS_UTILS,_settings.rundir,os.path.splitext(lib.f1.fname)[0],_settings.rundir,firstName,_settings.rundir, secondName) ,"Preprocess")
               else:
                  firstName = "%s.fastq"%(os.path.splitext(lib.f1.fname)[0])
                  secondName = "%s.fastq"%(os.path.splitext(lib.f2.fname)[0])
                  convertFastaToFastq("%s/Preprocess/in/%s"%(_settings.rundir, lib.f1.fname), "%s/Preprocess/in/%s.qual"%(_settings.rundir, lib.f1.fname), "%s/Preprocess/in/%s"%(_settings.rundir, firstName), "Preprocess")
                  convertFastaToFastq("%s/Preprocess/in/%s"%(_settings.rundir, lib.f2.fname), "%s/Preprocess/in/%s.qual"%(_settings.rundir, lib.f2.fname), "%s/Preprocess/in/%s"%(_settings.rundir, secondName), "Preprocess")
               commandOptions += " -o %s/Preprocess/out/%s -o %s/Preprocess/out/%s %s/Preprocess/in/%s %s/Preprocess/in/%s"%(_settings.rundir, firstName, _settings.rundir, secondName, _settings.rundir, firstName, _settings.rundir, secondName)
            else:
               firstName = "%s.fastq"%(os.path.splitext(lib.f1.fname)[0])
               secondName = ""
               convertFastaToFastq("%s/Preprocess/in/%s"%(_settings.rundir, lib.f1.fname), "%s/Preprocess/in/%s.qual"%(_settings.rundir, lib.f1.fname), "%s/Preprocess/in/%s"%(_settings.rundir, firstName), "Preprocess")
               commandOptions += " -o %s/Preprocess/out/%s %s/Preprocess/in/%s"%(_settings.rundir, firstName, _settings.rundir, firstName)
         elif lib.format == "fastq":
            if lib.mated:
               if lib.interleaved:
                  firstName = "%s.1.fastq"%(os.path.splitext(lib.f1.fname)[0])
                  secondName = "%s.2.fastq"%(os.path.splitext(lib.f1.fname)[0])
                  run_process(_settings, "perl %s/perl/split_fastq.pl %s/Preprocess/in/%s %s/Preprocess/in/%s %s/Preprocess/in/%s"%(_settings.METAMOS_UTILS,_settings.rundir,lib.f1.fname,_settings.rundir,firstName,_settings.rundir,secondName),"Preprocess")
               else:
                  firstName = lib.f1.fname
                  secondName = lib.f2.fname
               commandOptions += " -o %s/Preprocess/out/%s -o %s/Preprocess/out/%s %s/Preprocess/in/%s %s/Preprocess/in/%s"%(_settings.rundir, firstName, _settings.rundir, secondName, _settings.rundir, firstName, _settings.rundir, secondName)
            else:
               firstName = lib.f1.fname
               secondName = ""
               commandOptions += " -o %s/Preprocess/out/firstName %s/Preprocess/out/%s"%(_settings.rundir, firstName, _settings.rundir, firstName)
         else:
            print "*** metAMOS: EA-UTILS filtering does not currently support non fasta/fastq files, will not filter library %d"%(lib.id) 
            for read in lib.reads:
               run_process(_settings, "ln %s/Preprocess/in/%s %s/Preprocess/out/"%(_settings.rundir,read.fname,_settings.rundir),"Preprocess")
            continue

         lib.format = "fastq"
         lib.f1.fname = firstName
         lib.f1.path = "%s%s%s"%(os.path.dirname(lib.f1.path.replace("/in/","/out/")), os.sep, firstName)
         if lib.mated and lib.interleaved:
            lib.interleaved = False
            lib.f2 = Read("fastq","%s%s%s"%(os.path.dirname(lib.f1.path), os.sep,secondName), lib.mated,lib.interleaved)
            lib.reads.append(lib.f2)
         elif lib.mated:
               lib.f2.fname = secondName
               lib.f2.path = "%s%s%s"%(os.path.dirname(lib.f2.path.replace("/in/","/out/")), os.sep, secondName)
         # update the library info in our config file as well, otherwise restarts won't realize reads are fastq now
         updateLibInfo(_settings.rundir + os.sep + "pipeline.ini", lib)

      # done now run EA-UTILS
      run_process(_settings, "%s/fastq-mcf -q 10 -t 0.01 -P %d %s/adapter.fna %s"%(_settings.EAUTILS, 33 if offset.lower() == "sanger" else 64, _settings.DB_DIR, commandOptions), "Preprocess")

      # finally rename reads
      for lib in _readlibs:
         if lib.format != "fastq":
            continue
         for read in lib.reads:
            counter = 1
            rq = open(read.path, 'r')
            wq = open(read.path.replace("fastq", "renamed.fastq"),'w')
            lines = rq.xreadlines()
            for line in lines:
               seq = lines.next()
               qltheader = lines.next()
               qlt = lines.next()
               header = ""
               if lib.mated:
                  header = "lib%dr%d/%d"%(lib.id, counter, 1 if "1.fastq" in read.path else 2)
               else:
                  header = "lib%dr%d"%(lib.id, counter)
               
               wq.write("@%s\n"%(header))
               wq.write("%s"%(seq))
               wq.write("+%s\n"%(header))
               wq.write("%s"%(qlt))
               counter += 1
            rq.close()
            wq.close()
            read.fname = read.fname.replace("fastq", "renamed.fastq")
            read.path = read.path.replace("fastq", "renamed.fastq")

   elif _filter == "pbcr":
      # correct the sequences using self-correction
      # for now, we do not support correction using other libraries, we should in the future though

      if not os.path.exists("%s/pacBioToCA"%(_settings.CA)) or not os.path.exists("%s/blasr"%(_settings.BLASR)):
         print "*** metAMOS: Error: cannot correct without Celera Assembler and BLASR"
         raise(JobSignalledBreak)

      for lib in _readlibs:
         if lib.format != "fastq":
            print "*** metAMOS: Warning: can only correct fastq libraries"
            continue
         else:
            if lib.mated:
               print "*** metAMOS: Warning can only correct unpaired libraries"
               continue

            # first estimate genome size
            genomeSize = 5000000
            if os.path.exists("%s%skmergenie"%(_settings.KMERGENIE, os.sep)) and os.path.exists("%s%sR"%(_settings.R, os.sep)):
               run_process(_settings, "%s/kmergenie --one-pass -l 5 -k 20 -t %d -s 5 -o kmergenie %s"%(_settings.KMERGENIE, _settings.threads, lib.f1.path), "Preprocess")
               result = open("%s/Preprocess/out/%s_report.html"%(_settings.rundir, "kmergenie"), 'r')
               for line in result.xreadlines():
                  if "Predicted assembly size:" in line:
                     genomeSize = int(line.replace("<p><h4>Predicted assembly size:", "").replace("</h4></p>", "").replace("bp", "").strip())

            if _settings.VERBOSE:
               print "*** metAMOS: Correcting library %d with genome size estimate %s"%(lib.id, genomeSize)

            # run correction
            availableMem = getAvailableMemory(_settings)
            if availableMem == 0:
               availableMem = 8192
            else:
               availableMem *= 1024
            
            oldPath = os.environ["PATH"]
            os.environ["PATH"] = _settings.AMOS + os.pathsep + _settings.BLASR + os.pathsep + oldPath

            run_process(_settings, "%s/pacBioToCA -l lib%d -s %s/config/pacbio.blasr.spec -t %d -partitions 100 fastqFile=%s genomeSize=%s longReads=1 ovlThreads=%d merylThreads=%d cnsConcurrency=%d merylMemory=%s ovlStoreMemory=%s"%(_settings.CA, lib.id, _settings.METAMOS_UTILS, _settings.threads, lib.f1.path, genomeSize, _settings.threads, _settings.threads, _settings.threads, availableMem, availableMem), "Preprocess")
            # subset longest 25X?
            # update library format and file names
            lib.f1.path = "%s/Preprocess/out/lib%d.fastq"%(_settings.rundir, lib.id)
            lib.f1.fname = os.path.basename(lib.f1.path)
            lib.format = "fastq"
            lib.offset = "sanger"
            
   elif _filter == "khmer":
      print "*** metAMOS: Not supported yet, stay tuned"
      for lib in _readlibs:
          for read in lib.reads:
              run_process(_settings, "ln %s/Preprocess/in/%s %s/Preprocess/out/"%(_settings.rundir,read.fname,_settings.rundir),"Preprocess")
              if (lib.format == "fasta"):
                 run_process(_settings, "ln %s/Preprocess/in/%s.qual %s/Preprocess/out/"%(_settings.rundir,read.fname,_settings.rundir),"Preprocess")
   else:
       for lib in _readlibs:
           for read in lib.reads:
               run_process(_settings, "ln %s/Preprocess/in/%s %s/Preprocess/out/"%(_settings.rundir,read.fname,_settings.rundir),"Preprocess")
               if (lib.format == "fasta"): 
                  run_process(_settings, "ln %s/Preprocess/in/%s.qual %s/Preprocess/out/"%(_settings.rundir,read.fname,_settings.rundir),"Preprocess")

   #PUNT HERE
   for lib in _readlibs:
      if 1:
           #this means interleaved, single file
           if lib.format == "sff":
               run_process(_settings, "unlink %s/Preprocess/out/lib%d.sff"%(_settings.rundir, lib.id), "Preprocess")
               run_process(_settings, "ln %s/Preprocess/in/%s %s/Preprocess/out/lib%d.sff"%(_settings.rundir, lib.f1.fname, _settings.rundir, lib.id), "Preprocess")

               if "newbler" in asm:
                  if _run_fastqc:
                     print "Warning: FastQC cannot run on SFF files, skipping."
                     _run_fastqc = false
                  run_process(_settings, "ln %s/Preprocess/out/lib%d.sff %s/Preprocess/out/lib%d.seq"%(_settings.rundir, lib.id,_settings.rundir, lib.id), "Preprocess")
               else:
                  if not os.path.exists(_settings.CA + os.sep + "sffToCA"):
                     print "Error: CA not found in %s. It is needed to convert SFF files to fasta.\n"%(_settings.CA)
                     raise(JobSignalledBreak)

                  # generate the fasta files from the sff file
                  run_process(_settings, "rm -rf %s/Preprocess/out/%s.tmpStore"%(_settings.rundir, _settings.PREFIX), "Preprocess")
                  run_process(_settings, "rm -rf %s/Preprocess/out/%s.gkpStore"%(_settings.rundir, _settings.PREFIX), "Preprocess")
                  run_process(_settings, "unlink %s/Preprocess/out/lib%d.frg"%(_settings.rundir, lib.id), "Preprocess")
                  sffToCACmd = "%s/sffToCA -clear discard-n "%(_settings.CA)
                  if lib.linkerType != "flx":
                     sffToCACmd += "-clear 454 "
                  sffToCACmd += "-trim hard -libraryname lib%d -output %s/Preprocess/out/lib%d"%(lib.id, _settings.rundir, lib.id)
                  if (read.mated == True):
                      run_process(_settings, "%s -linker %s -insertsize %d %d %s/Preprocess/in/%s"%(sffToCACmd, lib.linkerType, lib.mean, lib.stdev, _settings.rundir, lib.f1.fname),"Preprocess")
                  else:
                      run_process(_settings, "%s %s/Preprocess/in/%s"%(sffToCACmd, _settings.rundir, lib.f1.fname),"Preprocess")
                  run_process(_settings, "%s/gatekeeper -T -F -o %s/Preprocess/out/%s.gkpStore %s/Preprocess/out/lib%d.frg"%(_settings.CA, _settings.rundir, _settings.PREFIX, _settings.rundir, lib.id),"Preprocess")
                  run_process(_settings, "%s/gatekeeper -dumpnewbler %s/Preprocess/out/lib%d %s/Preprocess/out/%s.gkpStore"%(_settings.CA, _settings.rundir, lib.id, _settings.rundir, _settings.PREFIX),"Preprocess")
                  run_process(_settings, "%s/gatekeeper -dumpfastq   %s/Preprocess/out/lib%d %s/Preprocess/out/%s.gkpStore"%(_settings.CA, _settings.rundir, lib.id, _settings.rundir, _settings.PREFIX), "Preprocess")
                  run_process(_settings, "%s/gatekeeper -dumplibraries -tabular %s/Preprocess/out/%s.gkpStore |awk '{if (match($3, \"U\") == 0 && match($1, \"UID\") == 0) print \"library\t\"$1\"\t\"$4-$5*3\"\t\"$4+$5*3}' >> %s/Preprocess/out/all.seq.mates"%(_settings.CA, _settings.rundir, _settings.PREFIX, _settings.rundir),"Preprocess")
                  run_process(_settings, "%s/gatekeeper -dumpfragments -tabular %s/Preprocess/out/%s.gkpStore|awk '{if ($3 != 0 && match($1, \"UID\")==0 && $1 < $3) print $1\"\t\"$3\"\t\"$5}' >> %s/Preprocess/out/all.seq.mates"%(_settings.CA, _settings.rundir, _settings.PREFIX, _settings.rundir),"Preprocess")
                  run_process(_settings, "%s/gatekeeper -dumpfragments -tabular %s/Preprocess/out/%s.gkpStore|awk '{if ($3 != 0 && match($1, \"UID\")==0 && $1 < $3) print $1\"\t\"$3}' > %s/Preprocess/out/lib%d.seq.mates"%(_settings.CA, _settings.rundir, _settings.PREFIX, _settings.rundir, lib.id), "Preprocess")
                  run_process(_settings, "unlink %s/Preprocess/out/lib%d.fasta"%(_settings.rundir,lib.id),"Preprocess")
                  run_process(_settings, "unlink %s/Preprocess/out/lib%d.fasta.qual"%(_settings.rundir,lib.id),"Preprocess")
                  run_process(_settings, "unlink %s/Preprocess/out/lib%d.fastq"%(_settings.rundir,lib.id),"Preprocess")
                  run_process(_settings, "unlink %s/Preprocess/out/lib%d.seq"%(_settings.rundir,lib.id),"Preprocess")
                  run_process(_settings, "ln %s/Preprocess/out/lib%d.fna %s/Preprocess/out/lib%d.fasta"%(_settings.rundir,lib.id,_settings.rundir,lib.id),"Preprocess")
                  run_process(_settings, "ln %s/Preprocess/out/lib%d.fna.qual %s/Preprocess/out/lib%d.fasta.qual"%(_settings.rundir,lib.id,_settings.rundir,lib.id),"Preprocess")
                  run_process(_settings, "ln %s/Preprocess/out/lib%d.fasta %s/Preprocess/out/lib%d.seq"%(_settings.rundir, lib.id, _settings.rundir,lib.id),"Preprocess")
                  if lib.mated:
                     run_process(_settings, "perl %s/perl/shuffleSequences_fastq.pl  %s/Preprocess/out/lib%d.1.fastq %s/Preprocess/out/lib%d.2.fastq %s/Preprocess/out/lib%d.seq"%(_settings.METAMOS_UTILS,_settings.rundir,lib.id, _settings.rundir,lib.id,_settings.rundir, lib.id), "Preprocess")
                  else:
                     run_process(_settings, "ln %s/Preprocess/out/lib%d.unmated.fastq %s/Preprocess/out/lib%d.fastq"%(_settings.rundir, lib.id, _settings.rundir, lib.id), "Preproces")

                  run_process(_settings, "rm -rf %s/Preproces/out/%s.gkpStore"%(_settings.rundir, _settings.PREFIX),"Preprocess")
                  run_process(_settings, "cat %s/Preprocess/out/lib%d.seq.mates >> %s/Preprocess/out/all.seq.mates"%(_settings.rundir, lib.id, _settings.rundir), "Preprocess")

                  if _asm.lower() != "ca" and _asm.lower() != "newbler":
                     #flip the type to fastq
                     lib.format = "fastq"
                     lib.interleaved = False
                     if lib.mated:
                        lib.f1 = Read(lib.format,"%s/Preprocess/out/lib%d.1.fastq"%(_settings.rundir, lib.id),lib.mated,lib.interleaved) 
                        lib.f2 = Read(lib.format,"%s/Preprocess/out/lib%d.2.fastq"%(_settings.rundir, lib.id),lib.mated,lib.interleaved) 
                     else:
                        run_process(_settings, "unlink %s/Preprocess/out/lib%d.seq"%(_settings.rundir,lib.id),"Preprocess")
                        run_process(_settings, "ln %s/Preprocess/out/lib%d.unmated.fastq %s/Preprocess/out/lib%d.seq"%(_settings.rundir, lib.id, _settings.rundir, lib.id), "Preproces")
                        lib.f1 = Read(lib.format,"%s/Preprocess/out/lib%d.seq"%(_settings.rundir, lib.id),lib.mated,lib.interleaved)  


           elif lib.format == "fasta" and not lib.mated:
               run_process(_settings, "ln %s/Preprocess/out/%s %s/Preprocess/out/lib%d.seq"%(_settings.rundir,lib.f1.fname,_settings.rundir,lib.id),"Preprocess")
               run_process(_settings, "ln %s/Preprocess/out/%s.qual %s/Preprocess/out/lib%d.seq.qual"%(_settings.rundir,lib.f1.fname,_settings.rundir,lib.id),"Preprocess")
               run_process(_settings, "touch %s/Preprocess/out/lib%d.seq.mates"%(_settings.rundir,lib.id),"Preprocess")
               convertInputFastaToFastq(lib.id, lib.mated)

           elif lib.format == "fastq" and not lib.mated:
               run_process(_settings, "ln %s/Preprocess/out/%s %s/Preprocess/out/lib%d.seq"%(_settings.rundir, lib.f1.fname, _settings.rundir, lib.id), "Preprocess")
               run_process(_settings, "touch %s/Preprocess/out/lib%d.seq.mates"%(_settings.rundir, lib.id), "Preprocess")
               convertInputFastqToFasta(lib.id, lib.mated)

           elif lib.format == "fasta" and lib.mated and not lib.interleaved:
               #FIXME, make me faster!filter
               run_process(_settings, "perl %s/perl/shuffleSequences_fasta.pl  %s/Preprocess/out/%s %s/Preprocess/out/%s %s/Preprocess/out/lib%d.seq"%(_settings.METAMOS_UTILS,_settings.rundir,lib.f1.fname, _settings.rundir,lib.f2.fname,_settings.rundir,lib.id),"Preprocess")
               run_process(_settings, "perl %s/perl/shuffleSequences_fasta.pl  %s/Preprocess/out/%s.qual %s/Preprocess/out/%s.qual %s/Preprocess/out/lib%d.seq.qual"%(_settings.METAMOS_UTILS,_settings.rundir,lib.f1.fname, _settings.rundir,lib.f2.fname,_settings.rundir,lib.id),"Preprocess")
               run_process(_settings, "ln %s/Preprocess/out/%s %s/Preprocess/out/lib%d.1.fasta"%(_settings.rundir, lib.f1.fname, _settings.rundir, lib.id), "Preprocess")
               run_process(_settings, "ln %s/Preprocess/out/%s %s/Preprocess/out/lib%d.2.fasta"%(_settings.rundir, lib.f2.fname, _settings.rundir, lib.id), "Preprocess")
               run_process(_settings, "ln %s/Preprocess/out/%s.qual %s/Preprocess/out/lib%d.1.fasta.qual"%(_settings.rundir, lib.f1.fname, _settings.rundir, lib.id), "Preprocess")
               run_process(_settings, "ln %s/Preprocess/out/%s.qual %s/Preprocess/out/lib%d.2.fasta.qual"%(_settings.rundir, lib.f2.fname, _settings.rundir, lib.id), "Preprocess")
               run_process(_settings, "which python","Preprocess")
               run_process(_settings, "echo $PYTHONPATH","Preprocess")
               #run_process(_settings, "python %s/python/extract_mates_from_fasta.py %s/Preprocess/out/lib%d.seq"%(_settings.METAMOS_UTILS,_settings.rundir,lib.id),"Preprocess")
               extract_mates_from_fasta("%s/Preprocess/out/lib%d.seq"%(_settings.rundir,lib.id))
               convertInputFastaToFastq(lib.id, lib.mated)

           elif lib.format == "fastq" and lib.mated and not lib.interleaved:
               #extract mates from fastq
               run_process(_settings, "perl %s/perl/shuffleSequences_fastq.pl  %s/Preprocess/out/%s %s/Preprocess/out/%s %s/Preprocess/out/lib%d.seq"%(_settings.METAMOS_UTILS,_settings.rundir,lib.f1.fname, _settings.rundir,lib.f2.fname,_settings.rundir,lib.id),"Preprocess")
               run_process(_settings, "ln %s/Preprocess/out/%s %s/Preprocess/out/lib%d.1.fastq"%(_settings.rundir, lib.f1.fname, _settings.rundir, lib.id), "Preprocess")
               run_process(_settings, "ln %s/Preprocess/out/%s %s/Preprocess/out/lib%d.2.fastq"%(_settings.rundir, lib.f2.fname, _settings.rundir, lib.id), "Preprocess")
               run_process(_settings, "which python","Preprocess")
               run_process(_settings, "echo $PYTHONPATH","Preprocess")
               #run_process(_settings, "python %s/python/extract_mates_from_fastq.py %s/Preprocess/out/lib%d.seq"%(_settings.METAMOS_UTILS,_settings.rundir,lib.id),"Preprocess")
               extract_mates_from_fastq("%s/Preprocess/out/lib%d.seq"%(_settings.rundir,lib.id))
               convertInputFastqToFasta(lib.id, lib.mated)

           elif lib.mated and lib.interleaved:
               run_process(_settings, "cp %s/Preprocess/out/%s %s/Preprocess/out/lib%d.seq"%(_settings.rundir,lib.f1.fname,_settings.rundir,lib.id),"Preprocess")
               if lib.format == "fastq":
                   run_process(_settings, "which python","Preprocess")
                   run_process(_settings, "echo $PYTHONPATH","Preprocess")
                   #run_process(_settings, "python %s/python/extract_mates_from_fastq.py %s/Preprocess/out/lib%d.seq"%(_settings.METAMOS_UTILS,_settings.rundir,lib.id),"Preprocess")
                   extract_mates_from_fastq("%s/Preprocess/out/lib%d.seq"%(_settings.rundir,lib.id))
                   # unshuffle the sequences
                   run_process(_settings, "perl %s/perl/split_fastq.pl %s/Preprocess/out/lib%d.seq %s/Preprocess/out/lib%d.1.fastq %s/Preprocess/out/lib%d.2.fastq"%(_settings.METAMOS_UTILS, _settings.rundir, lib.id, _settings.rundir, lib.id, _settings.rundir, lib.id), "Preprocess")
                   convertInputFastqToFasta(lib.id, lib.mated)
               else:
                   run_process(_settings, "cp %s/Preprocess/out/%s.qual %s/Preprocess/out/lib%d.seq.qual"%(_settings.rundir,lib.f1.fname,_settings.rundir,lib.id),"Preprocess")
 
                   run_process(_settings, "which python","Preprocess")
                   run_process(_settings, "echo $PYTHONPATH","Preprocess")
                   #run_process(_settings, "python %s/python/extract_mates_from_fasta.py %s/Preprocess/out/lib%d.seq"%(_settings.METAMOS_UTILS,_settings.rundir,lib.id),"Preprocess")
                   extract_mates_from_fasta("%s/Preprocess/out/lib%d.seq"%(_settings.rundir,lib.id))
                   # unshuffle the sequences
                   run_process(_settings, "perl %s/perl/split_fasta.pl %s/Preprocess/out/lib%d.seq %s/Preprocess/out/lib%d.1.fasta %s/Preprocess/out/lib%d.2.fasta"%(_settings.METAMOS_UTILS, _settings.rundir, lib.id, _settings.rundir, lib.id, _settings.rundir, lib.id), "Preprocess")
                   run_process(_settings, "perl %s/perl/split_fasta.pl %s/Preprocess/out/lib%d.seq.qual %s/Preprocess/out/lib%d.1.fasta.qual %s/Preprocess/out/lib%d.2.fasta.qual"%(_settings.METAMOS_UTILS, _settings.rundir, lib.id, _settings.rundir, lib.id, _settings.rundir, lib.id), "Preprocess")
                   convertInputFastaToFastq(lib.id, lib.mated)

           #update_soap_config()
           #elif _asm == "ca":
           #    #useful for 454, need to get SFF to FRG?
           #    #/fs/wrenhomes/sergek/wgs-assembler/Linux-amd64/bin/sffToCA
           #    pass
           #elif _asm == "amos":
           #    #call toAmos_new              
           #    pass
#   print "%s\n"%filter;
   if _run_fastqc:
       if not os.path.exists(_settings.FASTQC + os.sep + "fastqc"):
           print "Error: FastQC not found in %s. Please check your path and try again.\n"%(_settings.FASTQC)
           raise(JobSignalledBreak)
       fastqFiles = []
       reportNames = []
#       print "format: %s"%lib.format
       for lib in _readlibs:
           if lib.format == "fastq":
               fastqFiles.append("%s/Preprocess/out/%s"%(_settings.rundir, lib.f1.fname))
               if lib.mated == True and not lib.interleaved:
                   fastqFiles.append("%s/Preprocess/out/%s"%(_settings.rundir, lib.f2.fname))
           elif lib.format == "sff":
               # fastq files will have been generated but not assigned to f1 and f2
               if lib.mated == True:
                   fastqFiles.append("%s/Preprocess/out/lib%d.1.fastq"%(_settings.rundir, lib.id))
                   fastqFiles.append("%s/Preprocess/out/lib%d.2.fastq"%(_settings.rundir, lib.id))
               else:
                   fastqFiles.append("%s/Preprocess/out/lib%d.unmated.fastq"%(_settings.rundir, lib.id))
           elif lib.format == "fasta":
              # fastq files will have been generated
              if lib.mated and not lib.interleaved:
                 fastqFiles.append("%s/Preprocess/out/lib%d.1.fastq"%(_settings.rundir, lib.id))
                 fastqFiles.append("%s/Preprocess/out/lib%d.2.fastq"%(_settings.rundir, lib.id))
              else:
                 fastqFiles.append("%s/Preprocess/out/lib%d.fastq"%(_settings.rundir, lib.id))

           if lib.format == "fastq" or lib.format == "sff" or lib.format == "fasta":
               reportNames.append("lib%d.1"%(lib.id))
               if lib.mated == True and not lib.interleaved:
                   reportNames.append("lib%d.2"%(lib.id))
       if len(fastqFiles) == 0:
           print "Warning: no fastq files found in libs; FastQC will not run.\n"
       else:
           fastQCCmd = "%s%s%s -t %s %s\n"%(_settings.FASTQC, os.sep, "fastqc", _settings.threads, ' '.join(fastqFiles))
           #print "%s\n"%fastQCCmd
           run_process(_settings, fastQCCmd, "Preprocess")
           for i, fastqFile in enumerate(fastqFiles):
               reportFile = "%s_fastqc"%fastqFile.replace(".fastq", "")
               if os.path.exists(reportFile):
                   run_process(_settings, "mv %s %s/Preprocess/out/%s.fastqc"%(reportFile, _settings.rundir, reportNames[i]), "Preprocess")
               else:
                   print "Warning: FastQC did not generate %s.\n"%reportFile
                   #raise(JobSignalledBreak)

   # before quitting make sure all the reads were not filtered out
   haveReads = False
   for lib in _readlibs:
      fileSize = os.stat("%s/Preprocess/out/lib%d.seq"%(_settings.rundir, lib.id)).st_size
      if fileSize > 0:
         haveReads = True
      else:
         print "Warning: library %d has no sequences\n"%(lib.id)
   if haveReads == False:
      print "**ERROR**"
      print "All input sequences were empty\n"
      print "**ERROR**"
      print ""
      print ""
      raise (JobSignalledBreak)

   run_process(_settings, "touch %s/Preprocess/out/preprocess.success"%(_settings.rundir),"Preprocess")
