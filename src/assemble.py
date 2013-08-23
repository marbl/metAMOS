#!python

import os, sys, string, time, BaseHTTPServer, getopt, re, subprocess, webbrowser
from operator import itemgetter

from utils import *
from preprocess import Preprocess
sys.path.append(INITIAL_UTILS)
from ruffus import *

import generic

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
   run_process(_settings, "touch %s/Preprocess/out/all.seq.mates"%(_settings.rundir), "Assemble")

   # prepare trim and pair information
   run_process(_settings, "cat %s/Assemble/out/assembly/454TrimStatus.txt |grep -v Trimpoints | grep -v left |grep -v right |awk '{print $0}' | awk '{print $1\" \"$2}' > %s/Assemble/out/454TrimNoPairs.txt"%(_settings.rundir, _settings.rundir), "Assemble")
   run_process(_settings, "cat %s/Assemble/out/assembly/454TrimStatus.txt |grep left |sed s/_left//g |awk '{print $0}' | awk '{print $1\" \"$2}' > %s/Assemble/out/454TrimLeftPairs.txt"%(_settings.rundir, _settings.rundir), "Assemble")
   run_process(_settings, "cat %s/Assemble/out/assembly/454TrimStatus.txt |grep right |sed s/_right//g | awk '{print $0}' | awk '{print $1\" \"$2}' > %s/Assemble/out/454TrimRightPairs.txt"%(_settings.rundir, _settings.rundir), "Assemble")

   for lib in _readlibs:
       if lib.format == "sff":
          run_process(_settings, "unlink %s/Preprocess/out/lib%d.seq"%(_settings.rundir, lib.id), "Assemble")
          run_process(_settings, "%s/sfffile -i %s/Assemble/out/454TrimNoPairs.txt -t %s/Assemble/out/454TrimNoPairs.txt -o %s/Preprocess/out/lib%d.noPairs.sff %s/Preprocess/out/lib%d.sff"%(_settings.NEWBLER, _settings.rundir, _settings.rundir, _settings.rundir, lib.id, _settings.rundir, lib.id), "Assemble")
          run_process(_settings, "%s/sffinfo -s %s/Preprocess/out/lib%d.noPairs.sff > %s/Preprocess/out/lib%d.seq"%(_settings.NEWBLER, _settings.rundir, lib.id, _settings.rundir, lib.id), "Assemble")
          run_process(_settings, "%s/sfffile -i %s/Assemble/out/454TrimLeftPairs.txt -t %s/Assemble/out/454TrimLeftPairs.txt -o %s/Preprocess/out/lib%d.noPairs.sff %s/Preprocess/out/lib%d.sff"%(_settings.NEWBLER, _settings.rundir, _settings.rundir, _settings.rundir, lib.id, _settings.rundir, lib.id), "Assemble")
          run_process(_settings, "%s/sffinfo -s %s/Preprocess/out/lib%d.noPairs.sff |awk '{if (match($1, \">\") == 1) { print $1\"_left\"; } else { print $0; }}' >> %s/Preprocess/out/lib%d.seq"%(_settings.NEWBLER, _settings.rundir, lib.id, _settings.rundir, lib.id), "Assemble")
          run_process(_settings, "%s/sfffile -i %s/Assemble/out/454TrimRightPairs.txt -t %s/Assemble/out/454TrimRightPairs.txt -o %s/Preprocess/out/lib%d.noPairs.sff %s/Preprocess/out/lib%d.sff"%(_settings.NEWBLER, _settings.rundir, _settings.rundir, _settings.rundir, lib.id, _settings.rundir, lib.id), "Assemble")
          run_process(_settings, "%s/sffinfo -s %s/Preprocess/out/lib%d.noPairs.sff |awk '{if (match($1, \">\") == 1) { print $1\"_right\"; } else { print $0; }}' >> %s/Preprocess/out/lib%d.seq"%(_settings.NEWBLER, _settings.rundir, lib.id, _settings.rundir, lib.id), "Assemble")

          run_process(_settings, "echo \"library\t%s\t%d\t%d\" >> %s/Preprocess/out/all.seq.mates"%(lib.sid, lib.mmin, lib.mmax, _settings.rundir), "Assemble")
          run_process(_settings, "cat %s/Assemble/out/454TrimLeftPairs.txt |awk '{print $1\"_left\t\"$1\"_right\"}' > %s/Preprocess/out/lib%d.seq.mates"%(_settings.rundir, _settings.rundir, lib.id), "Assemble")
          run_process(_settings, "cat %s/Assemble/out/454TrimLeftPairs.txt |awk '{print $1\"_left\t\"$1\"_right\t%s\"}' >> %s/Preprocess/out/all.seq.mates"%(_settings.rundir, lib.sid, _settings.rundir), "Assemble")
          run_process(_settings, "rm %s/Preprocess/out/lib%d.noPairs.sff"%(_settings.rundir, lib.id), "Assemble")
       elif lib.mated:
          run_process(_settings, "echo \"library\t%s\t%d\t%d\" >> %s/Preprocess/out/all.seq.mates"%(lib.sid, lib.mmin, lib.mmax, _settings.rundir), "Assemble")
          run_process(_settings, "cat %s/Preprocess/out/lib%d.seq.mates |awk '{print $0\"\t%s\"}' >> %s/Preprocess/out/all.seq.mates"%(_settings.rundir, lib.id, lib.sid, _settings.rundir), "Assemble")

def getVelvetGCommand(velvetPath):
   CATEGORIES = 0.0;
   p = subprocess.Popen("%s/velveth | grep CATEGORIES"%(velvetPath), shell=True, stdin=None, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
   (checkStdout, checkStderr) = p.communicate()
   if checkStderr != "":
      print "Warning: Cannot determine Velvet number of supported libraries"
   else:
      mymatch = re.split('\s+=\s+', checkStdout.strip())
      if (len(mymatch) == 2 and mymatch[1] != None):
         CATEGORIES = float(mymatch[1])

   velvetgCommandLine = ""

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

   return velvetgCommandLine
 
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
         if lib.innie:
            velvethCommandLine += " -shortPaired%s "%(currLibString)
         else:
            velvethCommandLine += " -shortMatePaired%s "%(currLibString)
      else:
         velvethCommandLine += " -short%s "%(currLibString)
      velvethCommandLine += "%s/Preprocess/out/lib%d.seq "%(_settings.rundir, lib.id)
      currLibID += 1
      if (currLibID > 1):
         currLibString = "%d"%(currLibID) 

   # now run velveth
   run_process(_settings, "%s"%(velvethCommandLine), "Assemble")

   # now build velvetg command line
   velvetgCommandLine = "%s/velvetg %s/Assemble/out/ "%(velvetPath, _settings.rundir) 
   velvetgCommandLine += getVelvetGCommand(velvetPath)
   velvetgCommandLine += " %s"%(getProgramParams(_settings.METAMOS_UTILS, "%s.spec"%(name), "", "-"))
   velvetgCommandLine += " -read_trkg yes -scaffolding no -amos_file yes";
   run_process(_settings, "%s"%(velvetgCommandLine), "Assemble")

   # make symlinks
   run_process(_settings, "rm %s/Assemble/out/%s.afg"%(_settings.rundir, _settings.PREFIX), "Assemble")
   run_process(_settings, "ln %s/Assemble/out/velvet_asm.afg %s/Assemble/out/%s.afg"%(_settings.rundir, _settings.rundir, _settings.PREFIX),"Assemble")
   run_process(_settings, "rm %s/Assemble/out/%s.asm.contig"%(_settings.rundir, _settings.PREFIX),"Assemble")
   run_process(_settings, "ln %s/Assemble/out/contigs.fa %s/Assemble/out/%s.asm.contig"%(_settings.rundir, _settings.rundir, _settings.PREFIX), "Assemble")

def runSparseAssembler(sparsePath, name):
   if not os.path.exists(sparsePath + os.sep + "SparseAssembler"):
      print "Error: %s not found in %s. Please check your path and try again.\n"%(name, sparsePath)
      raise(JobSignalledBreak)

   sparseLibLine = ""
   libsAdded = 0
   currLibString = "";
   for lib in _readlibs:
      format = lib.format
      if lib.format == "fasta":
         print "Warning: sparse assembler requires fastq files, processing library %d as fastq\n"%(lib.id)
         format = "fastq"

      if format == "fastq":
         libsAdded += 1
         if lib.mated:
            run_process(_settings, "ln %s/Preprocess/out/lib%d.1.fastq %s/Assemble/out/lib%d.1.fastq"%(_settings.rundir,
 lib.id, _settings.rundir, lib.id), "Assemble")
            run_process(_settings, "ln %s/Preprocess/out/lib%d.2.fastq %s/Assemble/out/lib%d.2.fastq"%(_settings.rundir,
 lib.id, _settings.rundir, lib.id), "Assemble")
            sparseLibLine += "p1 lib%d.1.fastq p2 lib%d.2.fastq"%(lib.id, lib.id)
         else:
            run_process(_settings, "ln %s/Preprocess/out/lib%d.seq %s/Assemble/out/lib%d.seq"%(_settings,rundir, lib.id, settings_rundir, lib.id), "Assemble")
            sparseLibLine += "f lib%d.fastq"%(lib.id)

   if libsAdded == 0:
      print "Error: SparseAssembler was selected but no libraries are in fastq format. Cannot run assembly\n"
      raise(JobSignalledBreak)
 
   # now we can run the program, we start with a two-step correction
   run_process(_settings, "%s/ReadsDenoiser %s %s"%(sparsePath, getProgramParams(_settings.METAMOS_UTILS, "%s.spec"%(name), "ReadDenoiser", ""), sparseLibLine), "Assemble")
   sparseLibLine = sparseLibLine.replace("lib", "Denoised_lib");
   run_process(_settings, "%s/ReadsDenoiser %s %s"%(sparsePath, getProgramParams(_settings.METAMOS_UTILS, "%s.spec"%(name), "ReadDenoiserStep2", ""), sparseLibLine), "Assemble")
   sparseLibLine = sparseLibLine.replace("lib", "Denoised_lib");

   # now run the actual assembler
   run_process(_settings, "%s/SparseAssembler %s %s"%(sparsePath, getProgramParams(_settings.METAMOS_UTILS, "%s.spec"%(name), "SparseAssembler", ""), sparseLibLine), "Assemble")

   # create symlinks
   run_process(_settings, "rm %s/Assemble/out/%s.asm.contig"%(_settings.rundir, _settings.PREFIX),"Assemble")
   run_process(_settings, "ln %s/Assemble/out/Contigs.txt %s/Assemble/out/%s.asm.contig"%(_settings.rundir, _settings.rundir, _settings.PREFIX), "Assemble")

def runMetaVelvet(velvetPath, metavelvetPath, name):
   # check for metavelvet
   if not os.path.exists(metavelvetPath + os.sep + "meta-velvetg"):
      print "Error: %s not found in %s. Please check your path and try again.\n"%(name, velvetPath)
      raise(JobSignalledBreak)

   # run velvet first
   runVelvet(velvetPath, name)

   # now run the extras
   velvetgCommandLine = "%s/meta-velvetg %s/Assemble/out/ "%(metavelvetPath, _settings.rundir)
   velvetgCommandLine += getVelvetGCommand(velvetPath)
   velvetgCommandLine += " %s"%(getProgramParams(_settings.METAMOS_UTILS, "%s.spec"%(name), "", "-"))
   velvetgCommandLine += " -read_trkg yes -scaffolding no -amos_file yes"
   run_process(_settings, "%s"%(velvetgCommandLine), "Assemble")
   
   # get coverage peaks
   # metavelvet comes with a peak estimator but it is not packaged in v1.01 only in the git repository, why?
   p = subprocess.Popen("%s/scripts/scriptEstimatedCovMulti.py %s/Assemble/out/meta-velvetg.LastGraph-stats.txt |tail -n 1"%(metavelvetPath, _settings.rundir), shell=True, stdin=None, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
   (checkStdout, checkStderr) = p.communicate()
   if checkStderr != "":
      print "Warning: Cannot determine MetaVelvet coverages"
   else:
      checkStdout.strip()
      # finally re-run velvetg with the peaks
      velvetgCommandLine = "%s/meta-velvetg %s/Assemble/out/ "%(metavelvetPath, _settings.rundir)
      velvetgCommandLine += getVelvetGCommand(velvetPath)
      velvetgCommandLine += " -read_trkg yes -scaffolding no -amos_file yes"
      velvetgCommandLine += " -exp_covs %s"%(checkStdout.strip())
      run_process(_settings, "%s"%(velvetgCommandLine), "Assemble")

   # make symlinks
   run_process(_settings, "rm %s/Assemble/out/%s.afg"%(_settings.rundir, _settings.PREFIX), "Assemble")
   run_process(_settings, "ln %s/Assemble/out/meta-velvetg.asm.afg %s/Assemble/out/%s.afg"%(_settings.rundir, _settings.rundir, _settings.PREFIX),"Assemble")
   run_process(_settings, "rm %s/Assemble/out/%s.asm.contig"%(_settings.rundir, _settings.PREFIX),"Assemble")
   run_process(_settings, "ln %s/Assemble/out/meta-velvetg.contigs.fa %s/Assemble/out/%s.asm.contig"%(_settings.rundir, _settings.rundir, _settings.PREFIX), "Assemble")
       
@posttask(touch_file("%s/Logs/assemble.ok"%(_settings.rundir))) 
@files("%s/Preprocess/out/preprocess.success"%(_settings.rundir),["%s/Logs/assemble.ok"%(_settings.rundir)])
#@posttask(create_symlink,touch_file("completed.flag"))
@follows(Preprocess)
def Assemble(input,output):
   #pick assembler
   if "Assemble" in _skipsteps or "assemble" in _skipsteps:
      run_process(_settings, "touch %s/Logs/assemble.skip"%(_settings.rundir), "Assemble")
      return 0
   if _asm == "none" or _asm == None:
      pass
   elif _asm == "soapdenovo":
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

      if not os.path.exists(_settings.SOAPDENOVO + os.sep + "soap63"):
         print "Error: SOAPdenovo not found in %s. Please check your path and try again.\n"%(_settings.SOAPDENOVO)
         raise(JobSignalledBreak)

      soapOptions = getProgramParams(_settings.METAMOS_UTILS, "soap.spec", "pregraph", "-") 
      soapContigOptions = getProgramParams(_settings.METAMOS_UTILS, "soap.spec", "contig", "-")

      #start stopwatch
      if _settings.kmer > 63:
          
          run_process(_settings, "%s/soap127 pregraph -p %d %d %s -s %s/soapconfig.txt -o %s/Assemble/out/%s.asm"%(_settings.SOAPDENOVO, _settings.threads, _settings.kmer, soapOptions, _settings.rundir,_settings.rundir,_settings.PREFIX),"Assemble")#SOAPdenovo config.txt
          run_process(_settings, "%s/soap127 contig -g %s/Assemble/out/%s.asm %s"%(_settings.SOAPDENOVO,_settings.rundir,_settings.PREFIX, soapContigOptions),"Assemble")#SOAPdenovo config.txt
      else:
          
          run_process(_settings, "%s/SOAPdenovo-63mer pregraph -p %d -d -K %d %s -s %s/soapconfig.txt -o %s/Assemble/out/%s.asm"%(_settings.SOAPDENOVO, _settings.threads, _settings.kmer, soapOptions, _settings.rundir,_settings.rundir,_settings.PREFIX),"Assemble")#SOAPdenovo config.txt
          run_process(_settings, "%s/SOAPdenovo-63mer contig -g %s/Assemble/out/%s.asm -R -M 3"%(_settings.SOAPDENOVO, _settings.rundir,_settings.PREFIX),"Assemble")#SOAPdenovo config.txt
      #if OK, convert output to AMOS
   elif _asm == "metaidba":
      bowtie_mapping = 1
      for lib in _readlibs:
          if lib.format != "fasta"  or (lib.mated and not lib.interleaved):
              print "Warning: meta-IDBA requires reads to be in (interleaved) fasta format, converting library"
          #apparently connect = scaffold? need to convert fastq to interleaved fasta to run, one lib per run??
          #print "%s/metaidba --read %s/Preprocess/out/lib%d.fasta --output  %s/Assemble/out/%s.asm --mink 21 --maxk %d --cover 1 --connect"%(_settings.METAIDBA,_settings.rundir,lib.id,_settings.rundir,_settings.PREFIX,_settings.kmer)

          metaidbaOptions = getProgramParams(_settings.METAMOS_UTILS, "metaidba.spec", "", "--")
          run_process(_settings, "%s/metaidba --read %s/Preprocess/out/lib%d.fasta --output  %s/Assemble/out/%s.asm %s --maxk %d"%(_settings.METAIDBA,_settings.rundir,lib.id,_settings.rundir,_settings.PREFIX,metaidbaOptions,_settings.kmer),"Assemble")
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

      mated = False;
      for lib in _readlibs:
          if lib.mated:
             mated = True;

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
      run_process(_settings, "%s/toAmos -o %s/Assemble/out/%s.mates.afg -m %s/Preprocess/out/all.seq.mates -ace %s/Assemble/out/%s.ace"%(_settings.AMOS,_settings.rundir, _settings.PREFIX, _settings.rundir, _settings.rundir, _settings.PREFIX),"Assemble")
      # get info on EID/IIDs for contigs
      run_process(_settings, "cat %s/Assemble/out/%s.mates.afg | grep -A 3 \"{CTG\" |awk '{if (match($1, \"iid\") != 0) {IID = $1} else if (match($1, \"eid\") != 0) {print $1\" \"IID; } }'|sed s/eid://g |sed s/iid://g > %s/Assemble/out/454eidToIID"%(_settings.rundir, _settings.PREFIX, _settings.rundir),"Assemble")
      run_process(_settings, "java -cp %s convert454GraphToCTL %s/Assemble/out/454eidToIID %s/Assemble/out/assembly/454ContigGraph.txt > %s/Assemble/out/%s.graph.cte"%(_settings.METAMOS_JAVA, _settings.rundir, _settings.rundir, _settings.rundir, _settings.PREFIX),"Assemble")
      run_process(_settings, "cat %s/Assemble/out/%s.mates.afg %s/Assemble/out/%s.graph.cte > %s/Assemble/out/%s.afg"%(_settings.rundir, _settings.PREFIX, _settings.rundir, _settings.PREFIX, _settings.rundir, _settings.PREFIX),"Assemble")
    
      # make symlink for subsequent steps
      run_process(_settings, "rm %s/Assemble/out/%s.asm.contig"%(_settings.rundir, _settings.PREFIX),"Assemble")
      run_process(_settings, "ln %s/Assemble/out/assembly/454AllContigs.fna %s/Assemble/out/%s.asm.contig"%(_settings.rundir, _settings.rundir, _settings.PREFIX),"Assemble")
      if mated == True:
         run_process(_settings, "ln %s/Assemble/out/assembly/454Scaffolds.fna %s/Assemble/out/%s.asm.scafSeq"%(_settings.rundir, _settings.rundir, _settings.PREFIX),"Assemble")
      else:
         run_process(_settings, "ln %s/Assemble/out/assembly/454AllContigs.fna %s/Assemble/out/%s.asm.scafSeq"%(_settings.rundir, _settings.rundir, _settings.PREFIX),"Assemble")

   elif _asm == "amos":
      run_process(_settings, "rm -rf %s/Assemble/in/%s.bnk"%(_settings.rundir, _settings.PREFIX), "Assemble")
      for lib in _readlibs:
         if lib.format == "fasta":
            run_process(_settings, "%s/toAmos_new -s %s/Preprocess/out/lib%d.seq -b %s/Assemble/in/%s.bnk "%(_settings.AMOS,_settings.rundir,lib.id,_settings.rundir, _settings.PREFIX),"Assemble")
         elif lib.format == "fastq":
            run_process(_settings, "%s/toAmos_new -Q %s/Preprocess/out/lib%d.seq -i --libname lib%d --min %d --max %d -b %s/Assemble/in/%s.bnk "%(_settings.AMOS,_settings.rundir,lib.id,lib.id,lib.mean,lib.stdev,_settings.rundir,_settings.PREFIX),"Assemble")
      run_process(_settings, "%s/hash-overlap -B %s/Assemble/in/%s.bnk"%(_settings.AMOS, _settings.rundir, _settings.PREFIX), "Assemble")
      run_process(_settings, "%s/tigger -b %s/Assemble/in/%s.bnk"%(_settings.AMOS, _settings.rundir, _settings.PREFIX), "Assemble")
      run_process(_settings, "%s/make-consensus -B -b %s/Assemble/in/%s.bnk"%(_settings.AMOS, _settings.rundir, _settings.PREFIX), "Assemble")
      run_process(_settings, "%s/bank2fasta -b %s/Assemble/in/%s.bnk > %s.asm.contig"%(_settings.AMOS, _settings.rundir, _settings.PREFIX, _settings.PREFIX), "Assemble")
   elif _asm.lower() == "ca":
      #runCA script
      frglist = ""
      matedString = ""
      for lib in _readlibs:
         if lib.format == "fastq":
            if lib.mated:
               matedString = "-insertsize %d %d -%s -mates"%(lib.mean, lib.stdev, "innie" if lib.innie else "outtie") 
            else:
               matedString = "-reads"
            run_process(_settings, "%s/fastqToCA -libraryname %s -technology illumina %s %s/Preprocess/out/lib%d.seq > %s/Preprocess/out/lib%d.frg"%(_settings.CA, lib.sid, matedString, _settings.rundir, lib.id, _settings.rundir, lib.id),"Assemble")
         elif lib.format == "fasta":
            if lib.mated:
               matedString = "-mean %d -stddev %d -m %s/Preprocess/out/lib%d.seq.mates"%(lib.mean, lib.stdev, _settings.rundir, lib.id)
            run_process(_settings, "%s/convert-fasta-to-v2.pl -l %s %s -s %s/Preprocess/out/lib%d.seq -q %s/Preprocess/out/lib%d.seq.qual > %s/Preprocess/out/lib%d.frg"%(_settings.CA, lib.sid, matedString, _settings.rundir, lib.id, _settings.rundir, lib.id, _settings.rundir, lib.id),"Assemble")
         frglist += "%s/Preprocess/out/lib%d.frg "%(_settings.rundir, lib.id)

      run_process(_settings, "%s/runCA -p %s -d %s/Assemble/out/ -s %s/config/asm.spec %s"%(_settings.CA,_settings.PREFIX,_settings.rundir,_settings.METAMOS_UTILS,frglist),"Assemble")
      #convert CA to AMOS
      run_process(_settings, "%s/gatekeeper -dumpfrg -allreads %s.gkpStore > %s.frg"%(_settings.CA, _settings.PREFIX, _settings.PREFIX),"Assemble")
      run_process(_settings, "%s/terminator -g %s.gkpStore -t %s.tigStore/ 2 -o %s"%(_settings.CA, _settings.PREFIX, _settings.PREFIX, _settings.PREFIX),"Assemble")
      run_process(_settings, "%s/asmOutputFasta -p %s < %s.asm"%(_settings.CA, _settings.PREFIX, _settings.PREFIX), "Assemble")
      run_process(_settings, "ln %s.utg.fasta %s.asm.contig"%(_settings.PREFIX, _settings.PREFIX), "Assemble")
   elif _asm == "velvet":
      runVelvet(_settings.VELVET, "velvet")
   elif _asm == "velvet-sc":
      runVelvet(_settings.VELVET_SC, "velvet-sc")
   elif _asm == "metavelvet":
      runMetaVelvet(_settings.VELVET, _settings.METAVELVET, "metavelvet")
   elif _asm.lower() == "sparseassembler":
      runSparseAssembler(_settings.SPARSEASSEMBLER, "SparseAssembler");
   elif generic.checkIfExists(STEP_NAMES.ASSEMBLE, _asm.lower()):
      generic.execute(STEP_NAMES.ASSEMBLE, _asm.lower())
   else:  
      print "Error: %s is an unknown assembler. No valid assembler specified."%(_asm)
      raise(JobSignalledBreak)

   #if _usecontigs:
   #     map2contig()
   #stop here, for now
   #sys.exit(0)
   #check if sucessfully completed   
