#!python

import os, sys, string, time, BaseHTTPServer, getopt, re, subprocess, webbrowser
from operator import itemgetter

from utils import *

sys.path.append(INITIAL_UTILS)
from ruffus import *

_readlibs = []
_skipsteps = []
_asm = None
_run_fastqc = False
_settings = Settings() 

def init(reads, skipsteps, asm, run_fastqc):
   global _readlibs
   global _skipsteps
   global _asm
   global _run_fastqc

   _readlibs = reads
   _skipsteps = skipsteps
   _asm = asm
   _run_fastqc = run_fastqc

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
@files(_settings.readpaths,"%s/Preprocess/out/preprocess.success"%(_settings.rundir))
#filtreadpaths)
def Preprocess(input,output):
   global _run_fastqc

   #move input files into Preprocess ./in dir
   #output will either be split fastq files in out, or AMOS bank
   if "Preprocess" in _skipsteps or "preprocess" in _skipsteps:
       for lib in _readlibs:
           for read in lib.reads:
               run_process(_settings, "ln -s -t %s/Preprocess/out/ %s/Preprocess/in/%s"%(_settings.rundir,rundir,read.fname),"Preprocess")
       return 0
   run_process(_settings, "rm %s/Preprocess/out/all.seq.mates"%(_settings.rundir), "Preprocess")
   if filter == True:
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
                   #read.fname = "lib%d"%(lib.id)
                   read.filtered = True
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
               elif not read.filtered and read.format == "fasta" and read.mated and read.interleaved:
                   #this means we have this entire lib in one file
                   #parse out paired record (4 lines), rename header to be filename + "/1" or "/2", and remove reads with N
                   rf = open(read.path,'r')
                   npath = read.path.replace("/in/","/out/")
                   #print npath
                   #readpath,base = os.path.split(npath)
                   #newpath = readpath+"lib%d"%(lib.id)
                   wf = open(npath,'w')
                   #wf = open(read.path.replace("/in/","/out/"),'w')
                   start = 1
                   rcnt = 0
                   recordcnt = 0
                   record = []
                   shdr = ""
                   reads = rf.read().split(">")[1:]
                   if len(reads) % 2 != 0:
                       print "Read file corrupted, please fix and restart!"
                       sys.exit(1)

                   prevok = False
                   first = True
                   second = False
                   prevseq = ""
                   readcnt = 1
                   for rd in reads:
                       if first:
                           hdr,seq = rd.split("\n",1)
                           if "N" in string.upper(seq) or len(seq) < 2:
                               prevok = False
                           else: 
                               prevok = True
                               prevseq = seq

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
                               readcnt +=1
                           second = False
                           first = True

                   #update to new path
                   read.path = read.path.replace("/in/","/out/")            

                   #read.path = newpath#read.path.replace("/in/","/out/")            
                   #read.fname = "lib%d"%(lib.id)
               elif not read.filtered and read.format == "fasta" and read.mated and not read.interleaved:
                   readpair = lib.getPair(read.id)
                   if readpair == -1:
                       #not interleaved and mated, yet do not have 2nd file..
                       continue
                   rf1 = open(read.path,'r')
                   wf1 = open(read.path.replace("/in/","/out/"),'w')
                   rf2 = open(readpair.path,'r')
                   wf2 = open(readpair.path.replace("/in/","/out/"),'w')
                   recordcnt = 0
                   while 1:                   
                       rs1 = rf1.readline()
                       rs2 = rf1.readline()

                       if rs1 == "" or rs2 == "":
                           #EOF or something went wrong, break
                           break 
                       rseq = string.upper(rs2)                               
                       if "N" in rseq:
                           continue
                       rp1 = rf2.readline()
                       rp2 = rf2.readline()

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
                           wf2.writelines(">"+hdr+"2\n")
                           wf2.writelines(rp2)

                   readpair.filtered = True
                   read.filtered = True
                   read.path = read.path.replace("/in/","/out/")
                   readpair.path = readpair.path.replace("/in/","/out/")
               elif not read.filtered and read.format == "fasta" and not read.mated:
                   #easiest case, check for Ns
                   rf = open(read.path,'r')
                   wf = open(read.path.replace("/in/","/out/"),'w')
                   while 1:
                       rs1 = rf.readline()
                       rs2 = rf.readline()
                       if rs1 == "" or rs2 == "":
                           #EOF or something went wrong, break
                           break
                       rseq = string.upper(rs2)                               
                       if "N" in rseq:
                           continue
                       wf.writelines(rs1)
                       wf.writelines(rs2)
                   read.path = read.path.replace("/in/","/out/")
           cnt +=1
   else:
       for lib in _readlibs:
           for read in lib.reads:
               run_process(_settings, "ln -s -t %s/Preprocess/out/ %s/Preprocess/in/%s"%(_settings.rundir,_settings.rundir,read.fname),"Preprocess")
   #PUNT HERE
   for lib in _readlibs:
      if 1:
           #this means interleaved, single file
           if lib.format == "sff":
               run_process(_settings, "unlink %s/Preprocess/out/lib%d.sff"%(_settings.rundir, lib.id), "Preprocess")
               run_process(_settings, "ln -s %s %s/Preprocess/out/lib%d.sff"%(read.path, _settings.rundir, lib.id), "Preprocess")

               if _asm == "newbler":
                  if _run_fastqc:
                     print "Warning: FastQC cannot run on SFF files, skipping.";
                     _run_fastqc = false
                  run_process(_settings, "touch %s/Preprocess/out/lib%d.seq"%(_settings.rundir, lib.id), "Preprocess");
               else:
                  if not os.path.exists(_settings.CA + os.sep + "sffToCA"):
                     print "Error: CA not found in %s. It is needed to convert SFF files to fasta.\n"%(_settings.CA)
                     raise(JobSignalledBreak)

                  # generate the fasta files from the sff file
                  run_process(_settings, "rm -rf %s/Preprocess/out/%s.tmpStore"%(_settings.rundir, _settings.PREFIX), "Preprocess")
                  run_process(_settings, "rm -rf %s/Preprocess/out/%s.gkpStore"%(_settings.rundir, _settings.PREFIX), "Preprocess")
                  run_process(_settings, "unlink %s/Preprocess/out/%s.frg"%(_settings.rundir, _settings.PREFIX), "Preprocess")
                  sffToCACmd = "%s/sffToCA -clear discard-n "%(_settings.CA)
                  if lib.linkerType != "flx":
                     sffToCACmd += "-clear 454 "
                  sffToCACmd += "-trim hard -libraryname lib%d -output %s/Preprocess/out/lib%d"%(lib.id, _settings.rundir, lib.id)
                  if (read.mated == True):
                      run_process(_settings, "%s -linker %s -insertsize %d %d %s"%(sffToCACmd, lib.linkerType, lib.mean, lib.stdev, read.path),"Preprocess")
                  else:
                      run_process(_settings, "%s %s"%(sffToCACmd, read.path),"Preprocess")
                  run_process(_settings, "%s/gatekeeper -T -F -o %s/Preprocess/out/%s.gkpStore %s/Preprocess/out/lib%d.frg"%(_settings.CA, _settings.rundir, _settings.PREFIX, _settings.rundir, lib.id),"Preprocess")
                  run_process(_settings, "%s/gatekeeper -dumpnewbler %s/Preprocess/out/lib%d %s/Preprocess/out/%s.gkpStore"%(_settings.CA, _settings.rundir, lib.id, _settings.rundir, _settings.PREFIX),"Preprocess")
                  run_process(_settings, "%s/gatekeeper -dumpfastq   %s/Preprocess/out/lib%d %s/Preprocess/out/%s.gkpStore"%(_settings.CA, _settings.rundir, lib.id, _settings.rundir, _settings.PREFIX), "Preprocess")
                  run_process(_settings, "%s/gatekeeper -dumplibraries -tabular %s/Preprocess/out/%s.gkpStore |awk '{if (match($3, \"U\") == 0 && match($1, \"UID\") == 0) print \"library\t\"$1\"\t\"$4-$5*3\"\t\"$4+$5*3}' >> %s/Preprocess/out/all.seq.mates"%(_settings.CA, _settings.rundir, _settings.PREFIX, _settings.rundir),"Preprocess")
                  run_process(_settings, "%s/gatekeeper -dumpfragments -tabular %s/Preprocess/out/%s.gkpStore|awk '{if ($3 != 0 && match($1, \"UID\")==0 && $1 < $3) print $1\"\t\"$3\"\t\"$5}' >> %s/Preprocess/out/all.seq.mates"%(_settings.CA, _settings.rundir, _settings.PREFIX, _settings.rundir),"Preprocess")
                  run_process(_settings, "%s/gatekeeper -dumpfragments -tabular %s/Preprocess/out/%s.gkpStore|awk '{if ($3 != 0 && match($1, \"UID\")==0 && $1 < $3) print $1\"\t\"$3}' > %s/Preprocess/out/lib%d.seq.mates"%(_settings.CA, _settings.rundir, _settings.PREFIX, _settings.rundir, lib.id), "Preprocess")
                  run_process(_settings, "unlink %s/Preprocess/out/lib%d.seq"%(_settings.rundir,lib.id),"Preprocess")
                  run_process(_settings, "ln -s %s/Preprocess/out/lib%d.fna %s/Preprocess/out/lib%d.seq"%(_settings.rundir,lib.id,_settings.rundir,lib.id),"Preprocess")
                  run_process(_settings, "ln -s %s/Preprocess/out/lib%d.fna.qual %s/Preprocess/out/lib%d.seq.qual"%(_settings.rundir,lib.id,_settings.rundir,lib.id),"Preprocess")
                  run_process(_settings, "rm -rf %s/Preproces/out/%s.gkpStore"%(_settings.rundir, _settings.PREFIX),"Preprocess")
                  run_process(_settings, "cat %s/Preprocess/out/lib%d.seq.mates >> %s/Preprocess/out/all.seq.mates"%(_settings.rundir, lib.id, _settings.rundir), "Preprocess")

                  if _asm != "CA" and _asm != "newbler":
                     #flip the type to fastq
                     lib.format = "fastq"
                     lib.interleaved = False
                     if lib.mated:
                        lib.f1 = Read(lib.format,"%s/Preprocess/out/lib%d.1.fastq"%(_settings.rundir, lib.id),lib.mated,lib.interleaved) 
                        lib.f2 = Read(lib.format,"%s/Preprocess/out/lib%d.2.fastq"%(_settings.rundir, lib.id),lib.mated,lib.interleaved) 
                     else:
                        lib.f1 = Read(lib.format,"%s/Preprocess/out/lib%d.unmated.fastq"%(_settings.rundir, lib.id),lib.mated,lib.interleaved)  
           elif lib.format == "fasta" and not lib.mated:
               run_process(_settings, "ln -s %s/Preprocess/in/%s %s/Preprocess/out/lib%d.seq"%(_settings.rundir,lib.f1.fname,rundir,lib.id),"Preprocess")
               run_process(_settings, "ln -s %s/Preprocess/in/%s.qual %s/Preprocess/out/lib%d.seq.qual"%(_settings.rundir,lib.f1.fname,_settings.rundir,lib.id),"Preprocess")
               run_process(_settings, "touch %s/Preprocess/out/lib%d.seq.mates"%(_settings.rundir,lib.id),"Preprocess")
           elif lib.format == "fastq" and not lib.mated:
               run_process(_settings, "ln -s %s/Preprocess/in/%s %s/Preprocess/out/lib%d.seq"%(_settings.rundir, lib.f1.fname, _settings.rundir, lib.id), "Preprocess")
               run_process(_settings, "touch %s/Preprocess/out/lib%d.seq.mates"%(_settings.rundir, lib.id), "Preprocess")
           elif lib.format == "fasta" and lib.mated and not lib.interleaved:
               #FIXME, make me faster!filter
               run_process(_settings, "perl %s/perl/shuffleSequences_fasta.pl  %s/Preprocess/out/%s %s/Preprocess/out/%s %s/Preprocess/out/lib%d.seq"%(_settings.METAMOS_UTILS,_settings.rundir,lib.f1.fname, _settings.rundir,lib.f2.fname,rundir,lib.id),"Preprocess")
               run_process(_settings, "python %s/python/extract_mates_from_fasta.py %s/Preprocess/out/lib%d.seq"%(_settings.METAMOS_UTILS,_settings.rundir,lib.id),"Preprocess")
               run_process(_settings, "unlink %s/Preprocess/out/lib%d.seq.mates"%(_settings.rundir, lib.id),"Preprocess")
               run_process(_settings, "ln -t %s/Preprocess/out/ -s %s/Preprocess/in/lib%d.seq.mates"%(_settings.rundir,rundir,lib.id),"Preprocess")
           elif lib.format == "fastq" and lib.mated and not lib.interleaved:
               #extract mates from fastq
               run_process(_settings, "perl %s/perl/shuffleSequences_fastq.pl  %s/Preprocess/out/%s %s/Preprocess/out/%s %s/Preprocess/out/lib%d.seq"%(_settings.METAMOS_UTILS,_settings.rundir,lib.f1.fname, _settings.rundir,lib.f2.fname,_settings.rundir,lib.id),"Preprocess")
               run_process(_settings, "python %s/python/extract_mates_from_fastq.py %s/Preprocess/out/lib%d.seq"%(_settings.METAMOS_UTILS,_settings.rundir,lib.id),"Preprocess")
           elif lib.mated and lib.interleaved:
               run_process(_settings, "cp %s/Preprocess/out/%s %s/Preprocess/out/lib%d.seq"%(_settings.rundir,lib.f1.fname,_settings.rundir,lib.id),"Preprocess")
               if format == "fastq":
                   run_process(_settings, "python %s/python/extract_mates_from_fastq.py %s/Preprocess/out/lib%d.seq"%(_settings.METAMOS_UTILS,_settings.rundir,lib.id),"Preprocess")
               else:
                   run_process(_settings, "python %s/python/extract_mates_from_fasta.py %s/Preprocess/out/lib%d.seq"%(_settings.METAMOS_UTILS,_settings.rundir,lib.id),"Preprocess")
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
           if lib.format == "fastq" or lib.format == "sff":
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

   run_process(_settings, "touch %s/Preprocess/out/preprocess.success"%(_settings.rundir),"Preprocess")
