#!python

import os, sys, string, time, BaseHTTPServer, getopt, re, subprocess, webbrowser
from operator import itemgetter

from utils import *
from annotate import Annotate
from scaffold import Scaffold
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
   global _cls

   _readlibs = reads
   _skipsteps = skipsteps
   _cls = cls

@follows(Scaffold)
@files("%s/Annotate/out/%s.annots"%(_settings.rundir, _settings.PREFIX),"%s/Propagate/out/%s.clusters"%(_settings.rundir,_settings.PREFIX))
def Propagate(input,output):
   #run propogate java script
   # create s12.annots from Metaphyler output

   if _cls == "metaphyler":
       run_process(_settings, "python %s/python/create_mapping.py %s/DB/class_key.tab %s/Abundance/out/%s.classify.txt %s/Propagate/in/%s.annots"%(_settings.METAMOS_UTILS,_settings.METAMOS_UTILS,_settings.rundir,_settings.PREFIX,_settings.rundir,_settings.PREFIX),"Propagate")
   if _cls == "phylosift" or _cls == "PhyloSift" or _cls == "Phylosift" or _cls == "FCP" or _cls == "fcp":
       run_process(_settings, "ln -s %s/Annotate/out/%s.annots %s/Propagate/in/%s.annots"%(_settings.rundir,_settings.PREFIX,_settings.rundir,_settings.PREFIX),"Propagate")

   # strip headers from file and contig name prefix
   
   # some output from the classifiers (for example PhyloSift) outputs multiple contigs with the same classification on one line
   # the line looks like ctg1","ctg2 class so we don't know which is right and we skip it in the classification below
   run_process(_settings, "cat %s/Propagate/in/%s.annots | grep -v \"\\\"\" | grep -v contigID |sed s/utg//g |sed s/ctg//g > %s/Propagate/in/%s.clusters"%(_settings.rundir,_settings.PREFIX,_settings.rundir,_settings.PREFIX),"Propagate")

   if "Propagate" in _skipsteps or "propagate" in _skipsteps:
      run_process(_settings, "ln -s %s/Propagate/in/%s.clusters %s/Propagate/out/%s.clusters"%(_settings.rundir, _settings.PREFIX, _settings.rundir, _settings.PREFIX), "Propagate")
   else:
      run_process(_settings, "%s/FilterEdgesByCluster -b %s/Scaffold/in/%s.bnk -clusters %s/Propagate/in/%s.clusters -noRemoveEdges > %s/Propagate/out/%s.clusters"%(_settings.AMOS,_settings.rundir,_settings.PREFIX,_settings.rundir,_settings.PREFIX,_settings.rundir,_settings.PREFIX),"Propagate")
