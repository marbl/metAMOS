B77;10003;0c#!python

import os, sys, string, time, BaseHTTPServer, getopt, re, subprocess, webbrowser
from operator import itemgetter

from utils import *
from preprocess import Preprocess
from assemble import Assemble
sys.path.append(INITIAL_UTILS)
from ruffus import *

_readlibs = []
_skipsteps = []
_settings = Settings()
_asm = None
_mapper = "bowtie"
def init(reads, skipsteps, asm,mapper):
   global _readlibs
   global _asm
   global _skipsteps
   _mapper = mapper
   _readlibs = reads
   _skipsteps = skipsteps
   _asm = asm

@files("%s/Assemble/out/%s.bout"%(_settings.rundir,_settings.PREFIX))
#@posttask(create_symlink,touch_file("completed.flag"))
@follows(MapReads)
def CalcDist(input,output):

   if "CalcDist" in _skipsteps or "calcdist" in _skipsteps:
      return 0
   #given read pairs mapped to contigs, calc insert length
   
