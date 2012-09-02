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

if "Propagate" in _skipsteps or _cls == None:
   run_process(_settings, "touch %s/Propagate/out/%s.clusters"%(_settings.rundir, _settings.PREFIX), "Classify")


@follows(Propagate)
@posttask(touch_file("%s/Logs/classify.ok"%(_settings.rundir)))
@files("%s/Propagate/out/%s.clusters"%(_settings.rundir,_settings.PREFIX),"%s/Classify/out/sorted.txt"%(_settings.rundir))
def Classify(input,output):
   if "Classify" in _skipsteps or _cls == None:
       run_process(_settings, "touch %s/Logs/classify.skip"%(_settings.rundir), "Classify")       
       return 0
   if _cls == "FCP" or _cls == "fcp" or _cls == "phylosift" or _cls == "PhyloSift" or _cls == "Phylosift":
       run_process(_settings, "python %s/python/sort_contigs.py %s/Propagate/in/%s.clusters %s/Propagate/out/%s.clusters %s/DB/tax_key.tab %s/Classify/out %s/Scaffold/in/%s.bnk %s"%(_settings.METAMOS_UTILS, _settings.rundir, _settings.PREFIX, _settings.rundir, _settings.PREFIX, _settings.METAMOS_UTILS,_settings.rundir, _settings.rundir, _settings.PREFIX,_settings.AMOS),"Classify")

   else:
       run_process(_settings, "python %s/python/sort_contigs.py %s/Propagate/in/%s.clusters %s/Propagate/out/%s.clusters %s/DB/class_key.tab %s/Classify/out %s/Scaffold/in/%s.bnk %s"%(_settings.METAMOS_UTILS, _settings.rundir, _settings.PREFIX, _settings.rundir, _settings.PREFIX, _settings.METAMOS_UTILS,_settings.rundir, _settings.rundir, _settings.PREFIX,_settings.AMOS),"Classify")

