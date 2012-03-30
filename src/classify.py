#!python

import os, sys, string, time, BaseHTTPServer, getopt, re, subprocess, webbrowser
from operator import itemgetter

from utils import *
from propagate import Propagate

sys.path.append(INITIAL_UTILS)
from ruffus import *

_readlibs = []
_skipsteps = []
_cls = None
_settings = Settings()

def init(reads, skipsteps, cls):
   global _readlibs
   global _skipsteps
   global _cls

   _readlibs = reads
   _skipsteps = skipsteps
   _cls = cls

@follows(Propagate)
@files("%s/Propagate/out/%s.clusters"%(_settings.rundir,_settings.PREFIX),"%s/Classify/out/sorted.txt"%(_settings.rundir))
def Classify(input,output):
   if _cls == "phylosift" or _cls == "PhyloSift" or _cls == "Phylosift":
       run_process(_settings, "python %s/python/sort_contigs.py %s/Propagate/out/%s.clusters %s/DB/tax_key.tab %s/Classify/out %s/Scaffold/in/%s.bnk %s"%(_settings.METAMOS_UTILS, _settings.rundir, _settings.PREFIX, _settings.METAMOS_UTILS,_settings.rundir, _settings.rundir, _settings.PREFIX,_settings.AMOS),"Classify")

   else:
       run_process(_settings, "python %s/python/sort_contigs.py %s/Propagate/out/%s.clusters %s/DB/class_key.tab %s/Classify/out %s/Scaffold/in/%s.bnk %s"%(_settings.METAMOS_UTILS, _settings.rundir, _settings.PREFIX, _settings.METAMOS_UTILS,_settings.rundir, _settings.rundir, _settings.PREFIX,_settings.AMOS),"Classify")

