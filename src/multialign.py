#!python                                                                                                                                                                                          
import os, sys, string, time, BaseHTTPServer, getopt, re, subprocess, webbrowser
from operator import itemgetter

from utils import *
from assemble import Assemble
from findscforfs import FindScaffoldORFS

sys.path.append(INITIAL_UTILS)
from ruffus import *

_readlibs = []
_skipsteps = []
_forcesteps = []
_aln = None
_settings = Settings()
_refgenomes = ""

def init(reads, skipsteps, forcesteps, aln, refgenomes):
   global _readlibs
   global _skipsteps
   global _forcesteps
   global _refgenomes
   global _aln
   _readlibs = reads
   _skipsteps = skipsteps
   _forcesteps = forcesteps
   _aln = aln
   _refgenomes = refgenomes


@follows(Assemble)
@posttask(touch_file("%s/Logs/multialign.ok"%(_settings.rundir)))
@files("%s/Assemble/out/%s.asm.contig"%(_settings.rundir,_settings.PREFIX),"%s/MultiAlign/out/%s.tree"%(_settings.rundir,_settings.PREFIX))
def MultiAlign(input,output):
   if "MultiAlign" in _skipsteps:
      # can this be done automatically by ruffus pipeline?                                                                                                                                        
      run_process(_settings, "touch %s/Logs/multialign.skip"%(_settings.rundir), "MultiAlign")
      run_process(_settings, "touch %s/MultiAlign/out/%s.tree"%(_settings.rundir, _settings.PREFIX), "MultiAlign")
      return 0

   #run_process(_settings, "%s -p blastp -i %s/FindORFS/out/%s.faa -d %s/Abundance/out/markers.pfasta -m 8 -b 10 -v 10 -a %s -o %s/Abundance/out/%s.blastp"%(blastc, _settings.rundir,_settings.PREFIX,_settings.rundir,_settings.threads,_settings.rundir,_settings.PREFIX),"Abundance")

   pass



