#!python

import os, sys, string, time, BaseHTTPServer, getopt, re, subprocess, webbrowser
from operator import itemgetter

from utils import *
from propagate import Propagate

sys.path.append(INITIAL_UTILS)
from ruffus import *
from sort_contigs import *
_readlibs = []
_skipsteps = []
_cls = None
_lowmem = False
_settings = Settings()

def init(reads, skipsteps, cls, low):
   global _readlibs
   global _skipsteps
   global _cls
   global _lowmem

   _readlibs = reads
   _skipsteps = skipsteps
   _cls = cls
   _lowmem = low

@follows(Propagate)
@posttask(touch_file("%s/Logs/classify.ok"%(_settings.rundir)))
@files("%s/Propagate/out/%s.clusters"%(_settings.rundir,_settings.PREFIX),"%s/Classify/out/sorted.txt"%(_settings.rundir))
def Classify(input,output):
   if "Classify" in _skipsteps or _cls == None:
       run_process(_settings, "touch %s/Propagate/out/%s.clusters"%(_settings.rundir, _settings.PREFIX), "Classify")
       run_process(_settings, "touch %s/Logs/classify.skip"%(_settings.rundir), "Classify")       
       return 0

   if _cls == "FCP" or _cls == "fcp" or _cls == "phylosift" or _cls == "PhyloSift" or _cls == "Phylosift":
       #run_process(_settings, "python %s/python/sort_contigs.py %s/Propagate/in/%s.clusters %s/Propagate/out/%s.clusters %s/Propagate/out/%s.reads.clusters %s/tax_key.tab %s/FindORFS/out/%s.fna.bnk %s/FindORFS/out/%s.faa.bnk %s/FindORFS/out/%s.gene.map %s/Classify/out %s/Scaffold/in/%s.bnk %s"%(_settings.METAMOS_UTILS, _settings.rundir, _settings.PREFIX, _settings.rundir, _settings.PREFIX, _settings.rundir, _settings.PREFIX, _settings.DB_DIR, _settings.rundir, _settings.PREFIX, _settings.rundir, _settings.PREFIX, _settings.rundir, _settings.PREFIX, _settings.rundir, _settings.rundir, _settings.PREFIX,_settings.AMOS),"Classify")
       sort_contigs("%s/Propagate/in/%s.clusters"%(_settings.rundir,_settings.PREFIX),"%s/Propagate/out/%s.clusters"%(_settings.rundir,_settings.PREFIX),"%s/Propagate/out/%s.reads.clusters"%(_settings.rundir,_settings.PREFIX),"%s/tax_key.tab"%(_settings.DB_DIR),"%s/FindORFS/out/%s.fna.bnk"%(_settings.rundir,_settings.PREFIX),"%s/FindORFS/out/%s.faa.bnk"%(_settings.rundir,_settings.PREFIX),"%s/FindORFS/out/%s.gene.map"%(_settings.rundir,_settings.PREFIX),"%s/Classify/out"%(_settings.rundir),"%s/Scaffold/in/%s.bnk"%(_settings.rundir,_settings.PREFIX),"%s"%(_settings.AMOS), _lowmem)
   else:
       #run_process(_settings, "python %s/python/sort_contigs.py %s/Propagate/in/%s.clusters %s/Propagate/out/%s.clusters %s/Propagate/out/%s.reads.clusters %s/class_key.tab %s/FindORFS/out/%s.fna.bnk %s/FindORFS/out/%s.faa.bnk %s/FindORFS/out/%s.gene.map %s/Classify/out %s/Scaffold/in/%s.bnk %s"%(_settings.METAMOS_UTILS, _settings.rundir, _settings.PREFIX, _settings.rundir, _settings.PREFIX, _settings.rundir, _settings.PREFIX, _settings.DB_DIR, _settings.rundir, _settings.PREFIX, _settings.rundir, _settings.PREFIX, _settings.rundir, _settings.PREFIX, _settings.rundir, _settings.rundir, _settings.PREFIX,_settings.AMOS),"Classify")
       sort_contigs("%s/Propagate/in/%s.clusters"%(_settings.rundir,_settings.PREFIX),"%s/Propagate/out/%s.clusters"%(_settings.rundir,_settings.PREFIX),"%s/Propagate/out/%s.reads.clusters"%(_settings.rundir,_settings.PREFIX),"%s/class_key.tab"%(_settings.DB_DIR),"%s/FindORFS/out/%s.fna.bnk"%(_settings.rundir,_settings.PREFIX),"%s/FindORFS/out/%s.faa.bnk"%(_settings.rundir,_settings.PREFIX),"%s/FindORFS/out/%s.gene.map"%(_settings.rundir,_settings.PREFIX),"%s/Classify/out"%(_settings.rundir),"%s/Scaffold/in/%s.bnk"%(_settings.rundir,_settings.PREFIX),"%s"%(_settings.AMOS), _lowmem)
