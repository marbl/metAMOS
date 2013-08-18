#!python

import os, sys, string, time, BaseHTTPServer, getopt, re, subprocess, webbrowser
from operator import itemgetter

from utils import *
from findreps import FindRepeats
from annotate import Annotate
from fannotate import FunctionalAnnotation
sys.path.append(INITIAL_UTILS)
from ruffus import *

_readlibs = []
_skipsteps = []
_settings = Settings()
_retainBank = False
_asm = None
_mated = False

def init(reads, skipsteps, retainBank, asm):
   global _readlibs
   global _skipsteps
   global _retainBank
   global _asm
   global _mated

   _readlibs = reads
   _skipsteps = skipsteps
   _retainBank = retainBank
   _asm = asm.lower()

   for lib in _readlibs:
      if lib.mated == True:
         _mated = True
         break

@follows(FunctionalAnnotation)
@posttask(touch_file("%s/Logs/scaffold.ok"%(_settings.rundir)))
@files(["%s/Assemble/out/%s.asm.contig"%(_settings.rundir,_settings.PREFIX)],"%s/Scaffold/out/%s.scaffolds.final"%(_settings.rundir,_settings.PREFIX))
def Scaffold(input,output):
   if "Scaffold" in _skipsteps or "scaffold" in _skipsteps:
      run_process(_settings, "touch %s/Logs/scaffold.skip"%(_settings.rundir), "Scaffold")
      return 0
   global _retainBank

   # check if we need to do scaffolding
   numMates = 0

   # check if we need to retain the bank
   if not os.path.isdir("%s/Scaffold/in/%s.bnk"%(_settings.rundir, _settings.PREFIX)):
      _retainBank = False

   if not _retainBank:
       run_process(_settings, "rm -rf %s/Scaffold/in/%s.bnk"%(_settings.rundir,_settings.PREFIX),"Scaffold")
       if _asm == "newbler":
          p = subprocess.Popen("cat %s/Assemble/out/%s.graph.cte |grep \"{CTL\" |wc -l"%(_settings.rundir, _settings.PREFIX), stdin=None, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
          (checkStdout, checkStderr) = p.communicate()
          numMates = int(checkStdout.strip())

       if _asm == "newbler":
          run_process(_settings, "rm -rf %s/Scaffold/in/%s.bnk"%(_settings.rundir, _settings.PREFIX),"Scaffold")
          # build the bank for amos
          run_process(_settings, "%s/bank-transact -b %s/Scaffold/in/%s.bnk -c -m %s/Assemble/out/%s.afg"%(_settings.AMOS,_settings.rundir, _settings.PREFIX, _settings.rundir, _settings.PREFIX),"Scaffold")
       #elif _asm == "velvet" or _asm == "velvet-sc" or _asm == "metavelvet":
       #   run_process(_settings, "rm -rf %s/Scaffold/in/%s.bnk"%(_settings.rundir, _settings.PREFIX), "Scaffold")
       #   run_process(_settings, "%s/bank-transact -b %s/Scaffold/in/%s.bnk -c -m %s/Assemble/out/%s.afg"%(_settings.AMOS, _settings.rundir, _settings.PREFIX, _settings.rundir, _settings.PREFIX), "Scaffold")
       elif _asm == "ca" or _asm == "CA":
          run_process(_settings, "%s/toAmos_new -a %s/Assemble/out/%s.asm -f %s/Assemble/out/%s.frg -b %s/Scaffold/in/%s.bnk -U "%(_settings.AMOS, _settings.rundir, _settings.PREFIX, _settings.rundir, _settings.PREFIX, _settings.rundir, _settings.PREFIX),"Scaffold")
       else:
           for lib in _readlibs:
        
               if lib.format == "fasta":
                   run_process(_settings, "%s/toAmos_new -s %s/Preprocess/out/lib%d.seq -m %s/Assemble/out/%s.lib%d.mappedmates -b %s/Scaffold/in/%s.bnk "%(_settings.AMOS,_settings.rundir,lib.id,_settings.rundir, _settings.PREFIX,lib.id,_settings.rundir,_settings.PREFIX),"Scaffold")

               elif lib.format == "fastq":
                   matedStr = ""
                   if lib.mated:
                      matedStr = "-i --min %d --max %d --libname lib%d"%(lib.mmin, lib.mmax, lib.id) 
                   run_process(_settings, "%s/toAmos_new -Q %s/Preprocess/out/lib%d.seq %s -b %s/Scaffold/in/%s.bnk "%(_settings.AMOS,_settings.rundir,lib.id,matedStr,_settings.rundir,_settings.PREFIX),"Scaffold")

           run_process(_settings, "%s/toAmos_new -c %s/Assemble/out/%s.asm.tigr -b %s/Scaffold/in/%s.bnk "%(_settings.AMOS,_settings.rundir,_settings.PREFIX,_settings.rundir,_settings.PREFIX),"Scaffold")
   else:
       run_process(_settings, "perl %s/bank-unlock %s/Scaffold/in/%s.bnk"%(_settings.AMOS,_settings.rundir,_settings.PREFIX),"SCAFFOLD")
       run_process(_settings, "rm %s/Scaffold/in/%s.bnk/CTE.*"%(_settings.rundir,_settings.PREFIX),"SCAFFOLD")
       run_process(_settings, "rm %s/Scaffold/in/%s.bnk/CTL.*"%(_settings.rundir,_settings.PREFIX),"SCAFFOLD")
       run_process(_settings, "rm %s/Scaffold/in/%s.bnk/MTF.*"%(_settings.rundir,_settings.PREFIX),"SCAFFOLD")
       run_process(_settings, "rm %s/Scaffold/in/%s.bnk/SCF.*"%(_settings.rundir,_settings.PREFIX),"SCAFFOLD")

   # after the banks are created, skip the scaffolding when we have no mates
   if _mated == False and numMates == 0:
       print "No mate pair info available for scaffolding, skipping"
       run_process(_settings, "%s/bank2fasta -eid -b %s/Scaffold/in/%s.bnk > %s/Scaffold/out/%s.contigs"%(_settings.AMOS, _settings.rundir, _settings.PREFIX, _settings.rundir, _settings.PREFIX), "Scaffold")
       run_process(_settings, "ln %s/Scaffold/out/%s.contigs %s/Scaffold/out/%s.linearize.scaffolds.final"%(_settings.rundir, _settings.PREFIX, _settings.rundir, _settings.PREFIX), "Scaffold")
       #run_process(_settings, "touch %s/Scaffold/out/%s.linearize.scaffolds.final"%(_settings.rundir, _settings.PREFIX), "Scaffold")
       #_skipsteps.append("FindScaffoldORFS")
       #_skipsteps.append("Propagate")
       return 0

   #use asmQC to update mate insert lens
   run_process(_settings, "%s/asmQC -b %s/Scaffold/in/%s.bnk -scaff -recompute -update -numsd 2 "%(_settings.AMOS,_settings.rundir,_settings.PREFIX),"SCAFFOLD")
   run_process(_settings, "perl %s/bank-unlock %s/Scaffold/in/%s.bnk"%(_settings.AMOS,_settings.rundir,_settings.PREFIX),"SCAFFOLD")
   #calls to Bambus2, goBambus2 script
   # first, parse the parameters
   markRepeatParams = getProgramParams(_settings.METAMOS_UTILS, "bambus.spec", "MarkRepeats", "-")
   orientContigParams = getProgramParams(_settings.METAMOS_UTILS, "bambus.spec", "OrientContigs", "-")

   run_process(_settings, "%s/clk -b %s/Scaffold/in/%s.bnk"%(_settings.BAMBUS2,_settings.rundir,_settings.PREFIX),"Scaffold")
   run_process(_settings, "%s/Bundler -b %s/Scaffold/in/%s.bnk"%(_settings.BAMBUS2,_settings.rundir,_settings.PREFIX),"Scaffold")
   run_process(_settings, "%s/MarkRepeats %s -b %s/Scaffold/in/%s.bnk > %s/Scaffold/in/%s.reps"%(_settings.BAMBUS2,markRepeatParams,_settings.rundir,_settings.PREFIX,_settings.rundir,_settings.PREFIX),"Scaffold")
   run_process(_settings, "%s/OrientContigs %s -b %s/Scaffold/in/%s.bnk -repeats %s/Scaffold/in/%s.reps "%(_settings.BAMBUS2,orientContigParams,_settings.rundir,_settings.PREFIX, _settings.rundir, _settings.PREFIX),"Scaffold")

   # output results
   run_process(_settings, "%s/bank2fasta -eid  -b %s/Scaffold/in/%s.bnk > %s/Scaffold/out/%s.contigs"%(_settings.AMOS,_settings.rundir,_settings.PREFIX,_settings.rundir,_settings.PREFIX),"Scaffold")
   run_process(_settings, "%s/OutputMotifs -b %s/Scaffold/in/%s.bnk > %s/Scaffold/out/%s.motifs"%(_settings.BAMBUS2,_settings.rundir,_settings.PREFIX,_settings.rundir,_settings.PREFIX),"Scaffold")
   run_process(_settings, "%s/OutputResults -b %s/Scaffold/in/%s.bnk -p %s/Scaffold/out/%s "%(_settings.BAMBUS2,_settings.rundir,_settings.PREFIX,_settings.rundir,_settings.PREFIX),"Scaffold")
   run_process(_settings, "%s/OutputScaffolds -b %s/Scaffold/in/%s.bnk > %s/Scaffold/out/%s.scaffolds.final"%(_settings.BAMBUS2,_settings.rundir,_settings.PREFIX,_settings.rundir,_settings.PREFIX),"Scaffold")

   # generate linearize results
   run_process(_settings, "%s/Linearize -b %s/Scaffold/in/%s.bnk"%(_settings.BAMBUS2,_settings.rundir,_settings.PREFIX),"Scaffold")
   run_process(_settings, "%s/OutputResults -b %s/Scaffold/in/%s.bnk -p %s/Scaffold/out/%s.linearize "%(_settings.BAMBUS2,_settings.rundir,_settings.PREFIX,_settings.rundir,_settings.PREFIX),"Scaffold")
   run_process(_settings, "%s/OutputScaffolds -b %s/Scaffold/in/%s.bnk > %s/Scaffold/out/%s.linearize.scaffolds.final"%(_settings.BAMBUS2,_settings.rundir,_settings.PREFIX,_settings.rundir,_settings.PREFIX),"Scaffold")


