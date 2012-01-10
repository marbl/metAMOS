#!python

import os, sys, string, time, BaseHTTPServer, getopt, re, subprocess, webbrowser
from operator import itemgetter

from utils import *
from findreps import FindRepeats

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

   _readlibs = reads
   _skipsteps = skipsteps
   _retainBank = retainBank
   _asm = asm

   for lib in _readlibs:
      if lib.mated:
         _mated = True
         break

@follows(FindRepeats)
@files(["%s/Assemble/out/%s.asm.contig"%(_settings.rundir,_settings.PREFIX)],"%s/Scaffold/out/%s.scaffolds.final"%(_settings.rundir,_settings.PREFIX))
def Scaffold(input,output):
   # check if we need to do scaffolding
   numMates = 0
   if not _retainBank:
       run_process(_settings, "rm -rf %s/Scaffold/in/%s.bnk"%(_settings.rundir,_settings.PREFIX),"Scaffold")
       if _asm == "newbler":
          p = subprocess.Popen("cat %s/Assemble/out/%s.graph.cte |grep \"{CTL\" |wc -l"%(_settings.rundir, _settings.PREFIX), stdin=None, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
          (checkStdout, checkStderr) = p.communicate()
          numMates = int(checkStdout.strip())

       if _mated == False and numMates == 0:

          print "No mate pair info available for scaffolding, skipping"
          run_process(_settings, "touch %s/Scaffold/out/%s.linearize.scaffolds.final"%(_settings.rundir, _settings.PREFIX), "Scaffold")
          _skipsteps.append("FindScaffoldORFS")
          return 0

       if _asm == "soap":
           for lib in _readlibs:
        
               if lib.format == "fasta":
                   run_process(_settings, "%s/toAmos_new -s %s/Preprocess/out/lib%d.seq -m %s/Assemble/out/%s.lib%d.mappedmates -b %s/Scaffold/in/%s.bnk "%(_settings.AMOS,_settings.rundir,lib.id,_settings.rundir, _settings.PREFIX,lib.id,_settings.rundir,_settings.PREFIX),"Scaffold")

               elif format == "fastq":
                   run_process(_settings, "%s/toAmos_new -Q %s/Preprocess/out/lib%d.seq -m %s/Assemble/out/%s.lib%d.mappedmates -b %s/Scaffold/in/%s.bnk "%(_settings.AMOS,_settings.rundir,lib.id,_settings.rundir,_settings.PREFIX, lib.id,_settings.rundir,_settings.PREFIX),"Scaffold")

           run_process(_settings, "%s/toAmos_new -c %s/Assemble/out/%s.asm.tigr -b %s/Scaffold/in/%s.bnk "%(_settings.AMOS,_settings.rundir,_settings.PREFIX,_settings.rundir,_settings.PREFIX),"Scaffold")
       elif _asm == "metaidba":
          for lib in _readlibs:
              run_process(_settings, "%s/toAmos_new -s %s/Preprocess/out/lib%d.seq -m %s/Assemble/out/%s.lib%d.mappedmates -b %s/Scaffold/in/%s.bnk "%(_settings.AMOS,_settings.rundir,lib.id,_settings.rundir, _settings.PREFIX,lib.id,_settings.rundir,_settings.PREFIX),"Scaffold")
          run_process(_settings, "%s/toAmos_new -c %s/Assemble/out/%s.asm.tigr -b %s/Scaffold/in/%s.bnk "%(_settings.AMOS,_settings.rundir,_settings.PREFIX,_settings.rundir,_settings.PREFIX),"Scaffold")
       elif _asm == "newbler":
          run_process(_settings, "rm -rf %s/Scaffold/in/%s.bnk"%(_settings.rundir, _settings.PREFIX),"Scaffold")
          # build the bank for amos
          run_process(_settings, "%s/bank-transact -b %s/Scaffold/in/%s.bnk -c -m %s/Assemble/out/%s.afg"%(_settings.AMOS,_settings.rundir, _settings.PREFIX, _settings.rundir, _settings.PREFIX),"Scaffold")
       elif _asm == "velvet" or _asm == "velvet-sc":
          run_process(_settings, "rm -rf %s/Scaffold/in/%s.bnk"%(_settings.rundir, _settings.PREFIX), "Scaffold")
          run_process(_settings, "%s/bank-transact -b %s/Scaffold/in/%s.bnk -c -m %s/Assemble/out/%s.afg"%(_settings.AMOS, _settings.rundir, _settings.PREFIX, _settings.rundir, _settings.PREFIX), "Scaffold")
       elif _asm == "ca" or _asm == "CA":
          run_process(_settings, "%s/toAmos_new -a %s/Assemble/out/%s.asm -f %s/Assemble/out/%s.frg -b %s/Scaffold/in/%s.bnk -U "%(_settings.AMOS, _settings.rundir, _settings.PREFIX, _settings.rundir, _settings.PREFIX, _settings.rundir, _settings.PREFIX),"Scaffold")


   else:
       run_process(_settings, "%s/bank-unlock %s/Scaffold/in/%s.bnk"%(_settings.AMOS,_settings.rundir,_settings.PREFIX))
       run_process(_settings, "rm %s/Scaffold/in/%s.bnk/CTE.*"%(_settings.rundir,_settings.PREFIX),"SCAFFOLD")
       run_process(_settings, "rm %s/Scaffold/in/%s.bnk/CTL.*"%(_settings.rundir,_settings.PREFIX),"SCAFFOLD")
       run_process(_settings, "rm %s/Scaffold/in/%s.bnk/MTF.*"%(_settings.rundir,_settings.PREFIX),"SCAFFOLD")
       run_process(_settings, "rm %s/Scaffold/in/%s.bnk/SCF.*"%(_settings.rundir,_settings.PREFIX),"SCAFFOLD")

   #calls to Bambus2, goBambus2 script
   # first, parse the parameters
   markRepeatParams = getProgramParams(_settings.METAMOS_UTILS, "bambus.spec", "MarkRepeats", "-")
   orientContigParams = getProgramParams(_settings.METAMOS_UTILS, "bambus.spec", "OrientContigs", "-")

   run_process(_settings, "%s/clk -b %s/Scaffold/in/%s.bnk"%(_settings.AMOS,_settings.rundir,_settings.PREFIX),"Scaffold")
   run_process(_settings, "%s/Bundler -b %s/Scaffold/in/%s.bnk"%(_settings.AMOS,_settings.rundir,_settings.PREFIX),"Scaffold")
   run_process(_settings, "%s/MarkRepeats %s -b %s/Scaffold/in/%s.bnk > %s/Scaffold/in/%s.reps"%(_settings.AMOS,markRepeatParams,_settings.rundir,_settings.PREFIX,_settings.rundir,_settings.PREFIX),"Scaffold")
   run_process(_settings, "%s/OrientContigs %s -b %s/Scaffold/in/%s.bnk -repeats %s/Scaffold/in/%s.reps "%(_settings.AMOS,orientContigParams,_settings.rundir,_settings.PREFIX, _settings.rundir, _settings.PREFIX),"Scaffold")

   # output results
   run_process(_settings, "%s/bank2fasta  -b %s/Scaffold/in/%s.bnk > %s/Scaffold/out/%s.contigs"%(_settings.AMOS,_settings.rundir,_settings.PREFIX,_settings.rundir,_settings.PREFIX),"Scaffold")
   run_process(_settings, "%s/OutputMotifs -b %s/Scaffold/in/%s.bnk > %s/Scaffold/out/%s.motifs"%(_settings.AMOS,_settings.rundir,_settings.PREFIX,_settings.rundir,_settings.PREFIX),"Scaffold")
   run_process(_settings, "%s/OutputResults -b %s/Scaffold/in/%s.bnk -p %s/Scaffold/out/%s "%(_settings.AMOS,_settings.rundir,_settings.PREFIX,_settings.rundir,_settings.PREFIX),"Scaffold")
   run_process(_settings, "%s/OutputScaffolds -b %s/Scaffold/in/%s.bnk > %s/Scaffold/out/%s.scaffolds.final"%(_settings.AMOS,_settings.rundir,_settings.PREFIX,_settings.rundir,_settings.PREFIX),"Scaffold")

   # generate linearize results
   run_process(_settings, "%s/Linearize -b %s/Scaffold/in/%s.bnk"%(_settings.AMOS,_settings.rundir,_settings.PREFIX),"Scaffold")
   run_process(_settings, "%s/OutputResults -b %s/Scaffold/in/%s.bnk -p %s/Scaffold/out/%s.linearize "%(_settings.AMOS,_settings.rundir,_settings.PREFIX,_settings.rundir,_settings.PREFIX),"Scaffold")
   run_process(_settings, "%s/OutputScaffolds -b %s/Scaffold/in/%s.bnk > %s/Scaffold/out/%s.linearize.scaffolds.final"%(_settings.AMOS,_settings.rundir,_settings.PREFIX,_settings.rundir,_settings.PREFIX),"Scaffold")


