#!python

import os, sys, string, time, BaseHTTPServer, getopt, re, subprocess, webbrowser
from operator import itemgetter

from utils import *
from preprocess import Preprocess
sys.path.append(INITIAL_UTILS)
from ruffus import *

_readlibs = []
_skipsteps = []
_settings = Settings()
_asm = None
_usecontigs = False

def init(reads, skipsteps, asm, usecontigs):
   global _readlibs
   global _asm
   global _skipsteps
   global _usecontigs
   _readlibs = reads
   _skipsteps = skipsteps
   _asm = asm
   _usecontigs = usecontigs
   if _usecontigs:
       _asm = None

def extractNewblerReads():
   run_process(_settings, "unlink %s/Preprocess/out/all.seq.mates"%(_settings.rundir), "Assemble")

   # prepare trim and pair information
   run_process(_settings, "cat %s/Assemble/out/assembly/454TrimStatus.txt |grep -v Trimpoints | grep -v left |grep -v right |awk '{print $0}' | awk '{print $1\" \"$2}' > %s/Assemble/out/454TrimNoPairs.txt"%(_settings.rundir, _settings.rundir), "Assemble")
   run_process(_settings, "cat %s/Assemble/out/assembly/454TrimStatus.txt |grep left |sed s/_left//g |awk '{print $0}' | awk '{print $1\" \"$2}' > %s/Assemble/out/454TrimLeftPairs.txt"%(_settings.rundir, _settings.rundir), "Assemble")
   run_process(_settings, "cat %s/Assemble/out/assembly/454TrimStatus.txt |grep right |sed s/_right//g | awk '{print $0}' | awk '{print $1\" \"$2}' > %s/Assemble/out/454TrimRightPairs.txt"%(_settings.rundir, _settings.rundir), "Assemble")

   for lib in _readlibs:
       run_process(_settings, "unlink %s/Preprocess/out/lib%d.seq"%(_settings.rundir, lib.id), "Assemble")
       run_process(_settings, "%s/sfffile -i %s/Assemble/out/454TrimNoPairs.txt -t %s/Assemble/out/454TrimNoPairs.txt -o %s/Preprocess/out/lib%d.noPairs.sff %s/Preprocess/out/lib%d.sff"%(_settings.NEWBLER, _settings.rundir, _settings.rundir, _settings.rundir, lib.id, _settings.rundir, lib.id), "Assemble")
       run_process(_settings, "%s/sffinfo -s %s/Preprocess/out/lib%d.noPairs.sff > %s/Preprocess/out/lib%d.seq"%(_settings.NEWBLER, _settings.rundir, lib.id, _settings.rundir, lib.id), "Assemble")
       run_process(_settings, "%s/sfffile -i %s/Assemble/out/454TrimLeftPairs.txt -t %s/Assemble/out/454TrimLeftPairs.txt -o %s/Preprocess/out/lib%d.noPairs.sff %s/Preprocess/out/lib%d.sff"%(_settings.NEWBLER, _settings.rundir, _settings.rundir, _settings.rundir, lib.id, _settings.rundir, lib.id), "Assemble")
       run_process(_settings, "%s/sffinfo -s %s/Preprocess/out/lib%d.noPairs.sff |awk '{if (match($1, \">\") == 1) { print $1\"_left\"; } else { print $0; }}' >> %s/Preprocess/out/lib%d.seq"%(_settings.NEWBLER, _settings.rundir, lib.id, _settings.rundir, lib.id), "Assemble")
       run_process(_settings, "%s/sfffile -i %s/Assemble/out/454TrimRightPairs.txt -t %s/Assemble/out/454TrimRightPairs.txt -o %s/Preprocess/out/lib%d.noPairs.sff %s/Preprocess/out/lib%d.sff"%(_settings.NEWBLER, _settings.rundir, _settings.rundir, _settings.rundir, lib.id, _settings.rundir, lib.id), "Assemble")
       run_process(_settings, "%s/sffinfo -s %s/Preprocess/out/lib%d.noPairs.sff |awk '{if (match($1, \">\") == 1) { print $1\"_right\"; } else { print $0; }}' >> %s/Preprocess/out/lib%d.seq"%(_settings.NEWBLER, _settings.rundir, lib.id, _settings.rundir, lib.id), "Assemble")

       run_process(_settings, "cat %s/Assemble/out/454TrimLeftPairs.txt |awk '{print $1\"_left\t\"$1\"_right\"}' > %s/Preprocess/out/lib%d.seq.mates"%(_settings.rundir, _settings.rundir, lib.id), "Assemble")
       run_process(_settings, "echo \"library\t%s\t%d\t%d\" >> %s/Preprocess/out/all.seq.mates"%(lib.sid, lib.mmin, lib.mmax, _settings.rundir), "Assemble")
       run_process(_settings, "cat %s/Assemble/out/454TrimLeftPairs.txt |awk '{print $1\"_left\t\"$1\"_right\t%s\"}' >> %s/Preprocess/out/all.seq.mates"%(_settings.rundir, lib.sid, _settings.rundir), "Assemble")
       run_process(_settings, "rm %s/Preprocess/out/lib%d.noPairs.sff"%(_settings.rundir, lib.id), "Assemble")
 
def runVelvet(velvetPath, name):
   if not os.path.exists(velvetPath + os.sep + "velvetg"):
      print "Error: %s not found in %s. Please check your path and try again.\n"%(name, velvetPath)
      raise(JobSignalledBreak)

   CATEGORIES = 0.0;
   p = subprocess.Popen("%s/velveth | grep CATEGORIES"%(velvetPath), shell=True, stdin=None, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
   (checkStdout, checkStderr) = p.communicate()
   if checkStderr != "":
      print "Warning: Cannot determine Velvet number of supported libraries"
   else:
      mymatch = re.split('\s+=\s+', checkStdout.strip())
      if (len(mymatch) == 2 and mymatch[1] != None):
         CATEGORIES = float(mymatch[1])

   velvethCommandLine = "%s/velveth %s/Assemble/out/ %d "%(velvetPath,_settings.rundir,_settings.kmer) 
   currLibID = 1
   currLibString = "";
   for lib in _readlibs:
      if (currLibID > CATEGORIES):
         print "Warning: Velvet only supports %d libraries, will not input any more libraries\n"%(CATEGORIES)
         break 

      if lib.format == "fasta":
         velvethCommandLine += "-fasta"
      elif lib.format == "fastq":
         velvethCommandLine += "-fastq"

      if lib.mated:
         velvethCommandLine += " -shortPaired%s "%(currLibString)
      else:
         velvethCommandLine += " -short%s "%(currLibString)
      velvethCommandLine += "%s/Preprocess/out/lib%d.seq "%(_settings.rundir, lib.id)
      currLibID += 1
      if (currLibID > 1):
         currLibString = "%d"%(currLibID) 

   # now run velveth
   run_process(_settings, "%s"%(velvethCommandLine), "Assemble")

   # now build velvetg command line
   velvetgCommandLine = "%s/velvetg %s/Assemble/out/"%(velvetPath, _settings.rundir) 

   currLibID = 1;
   currLibString = ""
   for lib in _readlibs:
      if (currLibID > CATEGORIES):
         print "Warning: Velvet only supports %d libraries, will not input any more libraries\n"%(CATEGORIES)
         break

      if lib.mated:
         velvetgCommandLine += " -ins_length%s %d -ins_length%s_sd %d"%(currLibString, lib.mean, currLibString, lib.stdev)
      currLibID += 1
      if (currLibID > 1):
         currLibString = "%d"%(currLibID)
   velvetgCommandLine += " %s"%(getProgramParams(_settings.METAMOS_UTILS, "%s.spec"%(name), "", "-"))
   velvetgCommandLine += " -read_trkg yes -scaffolding no -amos_file yes";
   run_process(_settings, "%s"%(velvetgCommandLine), "Assemble")

   # make symlinks
   run_process(_settings, "rm %s/Assemble/out/%s.afg"%(_settings.rundir, _settings.PREFIX), "Assemble")
   run_process(_settings, "ln -s %s/Assemble/out/velvet_asm.afg %s/Assemble/out/%s.afg"%(_settings.rundir, _settings.rundir, _settings.PREFIX),"Assemble")
   run_process(_settings, "rm %s/Assemble/out/%s.asm.contig"%(_settings.rundir, _settings.PREFIX),"Assemble")
   run_process(_settings, "ln -s %s/Assemble/out/contigs.fa %s/Assemble/out/%s.asm.contig"%(_settings.rundir, _settings.rundir, _settings.PREFIX), "Assemble")

        
@files(_settings.asmfiles,["%s/Assemble/out/%s.asm.contig"%(_settings.rundir,_settings.PREFIX)])
#@posttask(create_symlink,touch_file("completed.flag"))
@follows(Preprocess)
def Assemble(input,output):
   #pick assembler
   if "Assemble" in _skipsteps or "assemble" in _skipsteps:
      return 0
   if _asm == "soap":
      #open & update config
      soapf = open("%s/config.txt"%(_settings.rundir),'r')
      soapd = soapf.read()
      soapf.close()
      cnt = 1
      libno = 1
      #print libs
      for lib in _readlibs:
          if (lib.format == "fastq" or lib.format == "fasta")  and lib.mated and not lib.interleaved:
              soapd = soapd.replace("LIB%dQ1REPLACE"%(lib.id),"%s/Preprocess/out/%s"%(_settings.rundir,lib.f1.fname))
              soapd = soapd.replace("LIB%dQ2REPLACE"%(lib.id),"%s/Preprocess/out/%s"%(_settings.rundir,lib.f2.fname))

          elif lib.format == "fastq"  and lib.mated and lib.interleaved:
              #this is NOT supported by SOAP, make sure files are split into two..
              #need to update lib.f2 path
              run_process(_settings, "perl %s/perl/split_fastq.pl %s/Preprocess/out/%s %s/Assemble/in/%s %s/Assemble/in/%s.f2"%(_settings.METAMOS_UTILS,_settings.rundir,lib.f1.fname,_settings.rundir,lib.f1.fname,_settings.rundir,lib.f1.fname),"Assemble")
              soapd = soapd.replace("LIB%dQ1REPLACE"%(lib.id),"%s/Assemble/in/%s"%(_settings.rundir,lib.f1.fname))
              soapd = soapd.replace("LIB%dQ2REPLACE"%(lib.id),"%s/Assemble/in/%s"%(_settings.rundir,lib.f1.fname+".f2"))

          elif lib.format == "fasta"  and lib.mated and lib.interleaved:
              soapd = soapd.replace("LIB%dQ1REPLACE"%(lib.id),"%s/Preprocess/out/%s"%(_settings.rundir,lib.f1.fname))
          else:
              soapd = soapd.replace("LIB%dQ1REPLACE"%(lib.id),"%s/Preprocess/out/%s"%(_settings.rundir,lib.f1.fname))

      #cnt +=1
      soapw = open("%s/soapconfig.txt"%(_settings.rundir),'w')
      soapw.write(soapd)
      soapw.close()

      if not os.path.exists(_settings.SOAP + os.sep + "soap63"):
         print "Error: SOAPdenovo not found in %s. Please check your path and try again.\n"%(_settings.SOAP)
         raise(JobSignalledBreak)

      print "Running SOAPdenovo on input reads..."
      soapOptions = getProgramParams(_settings.METAMOS_UTILS, "soap.spec", "", "-") 
      #start stopwatch
      if _settings.kmer > 63:
          
          run_process(_settings, "%s/soap127 all -p %d -R -d -K %d %s -s %s/soapconfig.txt -o %s/Assemble/out/%s.asm"%(_settings.SOAP, _settings.threads, _settings.kmer, soapOptions, _settings.rundir,_settings.rundir,_settings.PREFIX),"Assemble")#SOAPdenovo config.txt
      else:
          
          run_process(_settings, "%s/SOAPdenovo-63mer all -p %d -R -d -K %d %s -s %s/soapconfig.txt -o %s/Assemble/out/%s.asm"%(_settings.SOAP, _settings.threads, _settings.kmer, soapOptions, _settings.rundir,_settings.rundir,_settings.PREFIX),"Assemble")#SOAPdenovo config.txt

      #if OK, convert output to AMOS
   elif _asm == "metaidba":
      bowtie_mapping = 1
      for lib in _readlibs:
          if lib.format != "fasta"  or (lib.mated and not lib.interleaved):
              print "ERROR: meta-IDBA requires reads to be in (interleaved) fasta format, cannot run"
              sys.exit(1)
          #apparently connect = scaffold? need to convert fastq to interleaved fasta to run, one lib per run??
          #print "%s/metaidba --read %s/Preprocess/out/%s --output  %s/Assemble/out/%s.asm --mink 21 --maxk %d --cover 1 --connect"%(_settings.METAIDBA,_settings.rundir,lib.f1.fname,_settings.rundir,_settings.PREFIX,_settings.kmer)
          run_process(_settings, "%s/metaidba --read %s/Preprocess/out/%s --output  %s/Assemble/out/%s.asm --mink 21 --maxk %d --cover 1 --connect"%(_settings.METAIDBA,_settings.rundir,lib.f1.fname,_settings.rundir,_settings.PREFIX,_settings.kmer),"Assemble")
          run_process(_settings, "mv %s/Assemble/out/%s.asm-contig.fa %s/Assemble/out/%s.asm.contig"%(_settings.rundir,_settings.PREFIX,_settings.rundir,_settings.PREFIX),"Assemble")

   elif _asm == "newbler":
      if not os.path.exists(_settings.NEWBLER + os.sep + "newAssembly"):
         print "Error: Newbler not found in %s. Please check your path and try again.\n"%(_settings.NEWBLER)
         raise(JobSignalledBreak)

      run_process(_settings, "%s/newAssembly -force %s/Assemble/out"%(_settings.NEWBLER, _settings.rundir),"Assemble")

      NEWBLER_VERSION = 0.0;
      p = subprocess.Popen("%s/newAssembly --version | head -n 1"%(_settings.NEWBLER), shell=True, stdin=None, stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
      (checkStdout, checkStderr) = p.communicate()
      if checkStderr != "":
         print "Warning: Cannot determine Newbler version"
      else:
         mymatch = re.findall('\d+\.\d+', checkStdout.strip())
         if (len(mymatch) == 1 and mymatch[0] != None):
            NEWBLER_VERSION = float(mymatch[0])

      for lib in _readlibs:
          if lib.format == "fasta":
              run_process(_settings, "%s/addRun %s/Assemble/out %s/Preprocess/out/lib%d.seq"%(_settings.NEWBLER, _settings.rundir, _settings.rundir,lib.id),"Assemble")
          elif lib.format == "sff":
              run_process(_settings, "%s/addRun %s %s/Assemble/out %s/Preprocess/out/lib%d.sff"%(_settings.NEWBLER, ("-p" if lib.mated else ""), _settings.rundir, _settings.rundir, lib.id), "Assemble")
          elif lib.format == "fastq" and lib.interleaved:
              if (NEWBLER_VERSION < 2.6):
                 print "Error: FASTQ + Newbler only supported in Newbler version 2.6+. You are using version %s."%(_settings.NEWBLER_VERSION)
                 raise(JobSignalledBreak)
              run_process(_settings, "%s/addRun %s/Assemble/out %s/Preprocess/out/lib%d.seq"%(_settings.NEWBLER, _settings.rundir, _settings.rundir, lib.id),"Assemble")
          elif not lib.interleaved:
              print "Error: Only interleaved fastq files are supported for Newbler"
              raise(JobSignalledBreak)

      newblerCmd = "%s%srunProject "%(_settings.NEWBLER, os.sep)
      # read spec file to input to newbler parameters
      newblerCmd += getProgramParams(_settings.METAMOS_UTILS, "newbler.spec", "", "-")
      run_process(_settings, "%s -cpu %d %s/Assemble/out"%(newblerCmd,_settings.threads,_settings.rundir),"Assemble")

      # unlike other assemblers, we can only get the preprocess info for newbler after assembly (since it has to split sff files by mates)
      extractNewblerReads()

      # convert to AMOS
      run_process(_settings, "cat %s/Assemble/out/assembly/454Contigs.ace |awk '{if (match($2, \"\\\\.\")) {STR= $1\" \"substr($2, 1, index($2, \".\")-1); for (i = 3; i <=NF; i++) STR= STR\" \"$i; print STR} else { print $0} }' > %s/Assemble/out/%s.ace"%(_settings.rundir, _settings.rundir,_settings.PREFIX), "Assemble") 
      run_process(_settings, "%s/toAmos -o %s/Assemble/out/%s.mates.afg -m %s/Preprocess/out/all.seq.mates -ace %s/Assemble/out/%s.ace"%(AMOS,_settings.rundir, _settings.PREFIX, _settings.rundir, _settings.rundir, _settings.PREFIX),"Assemble")
      # get info on EID/IIDs for contigs
      run_process(_settings, "cat %s/Assemble/out/%s.mates.afg | grep -A 3 \"{CTG\" |awk '{if (match($1, \"iid\") != 0) {IID = $1} else if (match($1, \"eid\") != 0) {print $1\" \"IID; } }'|sed s/eid://g |sed s/iid://g > %s/Assemble/out/454eidToIID"%(_settings.rundir, _settings.PREFIX, _settings.rundir),"Assemble")
      run_process(_settings, "java -cp %s convert454GraphToCTL %s/Assemble/out/454eidToIID %s/Assemble/out/assembly/454ContigGraph.txt > %s/Assemble/out/%s.graph.cte"%(_settings.METAMOS_JAVA, _settings.rundir, _settings.rundir, _settings.rundir, _settings.PREFIX),"Assemble")
      run_process(_settings, "cat %s/Assemble/out/%s.mates.afg %s/Assemble/out/%s.graph.cte > %s/Assemble/out/%s.afg"%(_settings.rundir, _settings.PREFIX, _settings.rundir, _settings.PREFIX, _settings.rundir, _settings.PREFIX),"Assemble")
    
      # make symlink for subsequent steps
      run_process(_settings, "rm %s/Assemble/out/%s.asm.contig"%(_settings.rundir, _settings.PREFIX),"Assemble")
      run_process(_settings, "ln -s %s/Assemble/out/assembly/454AllContigs.fna %s/Assemble/out/%s.asm.contig"%(_settings.rundir, _settings.rundir, _settings.PREFIX),"Assemble")
      if mated == True:
         run_process(_settings, "ln -s %s/Assemble/out/assembly/454Scaffolds.fna %s/Assemble/out/%s.asm.scafSeq"%(_settings.rundir, _settings.rundir, _settings.PREFIX),"Assemble")
      else:
         run_process(_settings, "ln -s %s/Assemble/out/assembly/454AllContigs.fna %s/Assemble/out/%s.asm.scafSeq"%(_settings.rundir, _settings.rundir, _settings.PREFIX),"Assemble")

   elif _asm == "amos":
      run_process(_settings, "%s/Minimus %s/Preprocess/out/bank"%(AMOS,_settings.rundir),"Assemble")
   elif _asm == "CA" or _asm == "ca":
      #runCA script
      frglist = ""
      matedString = ""
      for lib in _readlibs:
          for read in lib.reads:
              if read.format == "fastq":
                  if lib.mated:
                      matedString = "-insertsize %d %d -%s"%(lib.mean, lib.stdev, "innie" if lib.innie else "outtie") 
                  run_process(_settings, "%s/fastqToCA %s -libraryname %s -t illumina -fastq %s/Preprocess/in/%s > %/Preprocess/out/lib%d.frg"%(CA, matedString, lib.read.path, _settings.rundir, _settings.PREFIX, _settings.rundir, lib.id),"Assemble")
              elif read.format == "fasta":
                  if lib.mated:
                      matedString = "-mean %d -stddev %d -m %s/Preprocess/out/lib%d.seq.mates"%(lib.mean, lib.stdev, lib.id)
                  run_process(_settings, "%s/convert-fasta-to-v2.pl -l %s %s -s %s/Preprocess/in/%s -q %s/Preprocess/in/%s.qual > %s/Preprocess/out/lib%d.frg"%(_settings.CA, lib.sid, matedString, _settings.rundir, read.fname, _settings.rundir, read.fname, _settings.rundir, read.fname,_settings.rundir),"Assemble")
              frglist += "%s/Preprocess/out/lib%d.frg"%(_settings.rundir, lib.id)
      run_process(_settings, "%s/runCA -p %s -d %s/Assemble/out/ -s %s/config/asm.spec %s"%(_settings.CA,_settings.PREFIX,_settings.rundir,_settings.METAMOS_UTILS,frglist),"Assemble")
      #convert CA to AMOS
      run_process(_settings, "%s/gatekeeper -dumpfrg -allreads %s.gkpStore > %s.frg"%(_settings.CA, _settings.PREFIX, _settings.PREFIX),"Assemble")
      run_process(_settings, "%s/terminator -g %s.gkpStore -t %s.tigStore/ 2 -o %s"%(_settings.CA, _settings.PREFIX, _settings.PREFIX, _settings.PREFIX),"Assemble")
      run_process(_settings, "%s/asmOutputFasta -p %s < %s.asm"%(_settings.CA, _settings.PREFIX, _settings.PREFIX), "Assemble")
      run_process(_settings, "ln -s %s.utg.fasta %s.asm.contig"%(_settings.PREFIX, _settings.PREFIX), "Assemble")
   elif _asm == "velvet":
      runVelvet(_settings.VELVET, "velvet")
   elif _asm == "velvet-sc":
      runVelvet(_settings.VELVET_SC, "velvet-sc")
   elif _asm == "none":
      pass
   else:  
      print "Error: %s is an unknown assembler. No valid assembler specified."%(_asm)
      raise(JobSignalledBreak)

   #if _usecontigs:
   #     map2contig()
   #stop here, for now
   #sys.exit(0)
   #check if sucessfully completed   
