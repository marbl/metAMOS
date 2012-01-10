#!python

import os, sys, string, time, BaseHTTPServer, getopt, re, subprocess, webbrowser
from operator import itemgetter

from utils import *
from annotate import Annotate
from findscforfs import FindScaffoldORFS

sys.path.append(INITIAL_UTILS)
from ruffus import *

_readlibs = []
_skipsteps = []
_cls = None
_settings = Settings()

def init(reads, skipsteps, cls):
   global _readlibs
   global _skipsteps

   _readlibs = reads
   _skipsteps = skipsteps
   _cls = cls

@follows(FindScaffoldORFS, Annotate)
@files("%s/DB/class_key.tab"%(_settings.METAMOS_UTILS),"%s/Propagate/out/%s.clusters"%(_settings.rundir,_settings.PREFIX))
def Propagate(input,output):
   #run propogate java script
   # create s12.annots from Metaphyler output
   if _cls == "metaphyler":
       run_process(_settings, "python %s/python/create_mapping.py %s/DB/class_key.tab %s/Abundance/out/%s.classify.txt %s/Propagate/in/%s.annots"%(_settings.METAMOS_UTILS,_settings.METAMOS_UTILS,_settings.rundir,_settings.PREFIX,_settings.rundir,_settings.PREFIX),"Propagate")
   if _cls == "amphora2" or _cls == "Amphora2" or _cls == "amphora":
       run_process(_settings, "cp %s/Assemble/out/%s.annots %s/Propagate/in/%s.annots"%(_settings.rundir,_settings.PREFIX,_settings.rundir,_settings.PREFIX),"Propagate")
   # strip headers from file and contig name prefix
   
   run_process(_settings, "cat %s/Propagate/in/%s.annots |sed s/contig_//g |grep -v contigID > %s/Propagate/in/%s.clusters"%(_settings.rundir,_settings.PREFIX,_settings.rundir,_settings.PREFIX),"Propagate")
   run_process(_settings, "%s/FilterEdgesByCluster -b %s/Scaffold/in/%s.bnk -clusters %s/Propagate/in/%s.clusters -noRemoveEdges > %s/Propagate/out/%s.clusters"%(_settings.AMOS,_settings.rundir,_settings.PREFIX,_settings.rundir,_settings.PREFIX,_settings.rundir,_settings.PREFIX),"Propagate")

