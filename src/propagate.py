#!python

import os, sys, string, time, BaseHTTPServer, getopt, re, subprocess, webbrowser
from operator import itemgetter

from utils import *
from annotate import Annotate
from scaffold import Scaffold
from findscforfs import FindScaffoldORFS
from abundance import Abundance

sys.path.append(INITIAL_UTILS)
from ruffus import *
from create_mapping import *
_readlibs = []
_skipsteps = []
_cls = None
_mated = False
_settings = Settings()

def init(reads, skipsteps, cls):
   global _readlibs
   global _skipsteps
   global _cls
   global _mated

   _readlibs = reads
   _skipsteps = skipsteps
   _cls = cls
   for lib in _readlibs:
      if lib.mated == True:
         _mated = True
         break

@follows(Abundance)
@posttask(touch_file("%s/Logs/propagate.ok"%(_settings.rundir)))
@files("%s/Annotate/out/%s.annots"%(_settings.rundir, _settings.PREFIX),"%s/Logs/propagate.ok"%(_settings.rundir))
def Propagate(input,output):
   if _cls == "metaphyler":
       #run_process(_settings, "python %s/python/create_mapping.py %s/class_key.tab %s/Abundance/out/%s.classify.txt %s/Propagate/in/%s.annots"%(_settings.METAMOS_UTILS,_settings.DB_DIR,_settings.rundir,_settings.PREFIX,_settings.rundir,_settings.PREFIX),"Propagate")
       create_mapping("%s/class_key.tab"%(_settings.DB_DIR),"%s/Abundance/out/%s.classify.txt"%(_settings.rundir,_settings.PREFIX),"%s/Propagate/in/%s.annots"%(_settings.rundir,_settings.PREFIX))
   else:
       run_process(_settings, "ln %s/Annotate/out/%s.annots %s/Propagate/in/%s.annots"%(_settings.rundir,_settings.PREFIX,_settings.rundir,_settings.PREFIX),"Propagate")

   # some output from the classifiers (for example PhyloSift) outputs multiple contigs with the same classification on one line
   # the line looks like ctg1","ctg2 class so we don't know which is right and we skip it in the classification below
   run_process(_settings, "cat %s/Propagate/in/%s.annots | grep -v \"\\\"\" | grep -v contigID |sed s/utg//g |sed s/ctg//g > %s/Propagate/in/%s.clusters"%(_settings.rundir,_settings.PREFIX,_settings.rundir,_settings.PREFIX),"Propagate")

   numMates = 0
   if os.path.exists("%s/Assemble/out/%s.graph.cte"%(_settings.rundir, _settings.PREFIX)): 
      p = subprocess.Popen("cat %s/Assemble/out/%s.graph.cte |grep \"{CTL\" |wc -l"%(_settings.rundir, _settings.PREFIX), stdin=None, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      (checkStdout, checkStderr) = p.communicate()
      numMates = int(checkStdout.strip())

   if "Propagate" in _skipsteps or "propagate" in _skipsteps or _cls == None or (_mated == False and numMates == 0):
       run_process(_settings, "touch %s/Logs/propagate.skip"%(_settings.rundir), "Propagate")
       run_process(_settings, "cp %s/Propagate/in/%s.clusters %s/Propagate/out/%s.clusters"%(_settings.rundir, _settings.PREFIX, _settings.rundir, _settings.PREFIX), "Propagate")
   else:
      run_process(_settings, "%s/FilterEdgesByCluster -b %s/Scaffold/in/%s.bnk -clusters %s/Propagate/in/%s.clusters -noRemoveEdges > %s/Propagate/out/%s.clusters"%(_settings.AMOS,_settings.rundir,_settings.PREFIX,_settings.rundir,_settings.PREFIX,_settings.rundir,_settings.PREFIX),"Propagate")

   # here we also propagate to the reads within contigs
   readctg_dict = {}
   for lib in _readlibs:
     ctgfile = open("%s/Assemble/out/%s.lib%dcontig.reads"%(_settings.rundir, _settings.PREFIX, lib.id), 'r')
     for line in ctgfile.xreadlines():
        line = line.replace("\n","")
        read, ctg = line.split()
        if ctg in readctg_dict:
           readctg_dict[ctg].append(read)
        else:
           readctg_dict[ctg] = [read,]
     ctgfile.close()

   read_annots = {}
   known_annots = {}
   annotsfile = open("%s/Propagate/out/%s.clusters"%(_settings.rundir, _settings.PREFIX), 'r')
   for line in annotsfile.xreadlines():
      line = line.replace("\n", "")
      ctg, annot = line.split()
      known_annots[ctg] = annot
   annotsfile.close()

   annotsfile = open("%s/Propagate/in/%s.annots"%(_settings.rundir, _settings.PREFIX), "r")
   for line in annotsfile.xreadlines():
     line = line.replace("\n", "")
     ctg, annot = line.split()
     if ctg not in readctg_dict.keys() and ctg not in known_annots.keys():
        read_annots[ctg] = annot
   annotsfile.close()

   if "Propagate" not in _skipsteps:
      annotsfile = open("%s/Propagate/out/%s.clusters"%(_settings.rundir, _settings.PREFIX), 'a')
      for ctg in read_annots:
          annotsfile.write("%s\t%s\n"%(ctg, read_annots[ctg]))
      annotsfile.close()

   annotsfile = open("%s/Propagate/out/%s.clusters"%(_settings.rundir, _settings.PREFIX), 'r')
   annotreads = open("%s/Propagate/out/%s.reads.clusters"%(_settings.rundir, _settings.PREFIX), 'w')
   for line in annotsfile.xreadlines():
     line = line.replace("\n", "")
     ctg, annot = line.split()
     if ctg in readctg_dict:
        for x in readctg_dict[ctg]:
           annotreads.write("%s\t%s\n"%(x, annot))
     else:
        annotreads.write("%s\t%s\n"%(ctg, annot))
   annotsfile.close()
   annotreads.close()
   readctg_dict.clear()
