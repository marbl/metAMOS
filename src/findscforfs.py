#!python

import os, sys, string, time, BaseHTTPServer, getopt, re, subprocess, webbrowser
from operator import itemgetter

from utils import *
from scaffold import Scaffold
from findorfs import parse_genemarkout

sys.path.append(INITIAL_UTILS)
from ruffus import *

_readlibs = []
_skipsteps = []
_settings = Settings()

def init(reads, skipsteps):
   global _readlibs
   global _skipsteps

   _readlibs = reads
   _skipsteps = skipsteps


@follows(Scaffold)
@files("%s/Scaffold/out/%s.linearize.scaffolds.final"%(_settings.rundir,_settings.PREFIX),"%s/FindScaffoldORFS/out/%s.scaffolds.orfs"%(_settings.rundir,_settings.PREFIX))
def FindScaffoldORFS(input,output):
   if "FindScaffoldORFS" in _skipsteps:
      run_process(_settings, "touch %s/FindScaffoldORFS/out/%s.scaffolds.faa"%(_settings.rundir, _settings.PREFIX))
      return 0

   run_process(_settings, "%s/gmhmmp -o %s/FindScaffoldORFS/out/%s.scaffolds.orfs -m %s/config/MetaGeneMark_v1.mod -d -a %s/Scaffold/out/%s.linearize.scaffolds.final"%(_settings.GMHMMP,_settings.rundir,_settings.PREFIX,_settings.METAMOS_UTILS,_settings.rundir,_settings.PREFIX),"FindScaffoldORFS")
   parse_genemarkout("%s/FindScaffoldORFS/out/%s.scaffolds.orfs"%(_settings.rundir,_settings.PREFIX),1, "FindScaffoldORFS")
   #run_process(_settings, "unlink %s/FindORFS/in/%s.scaffolds.faa"%(_settings.rundir,_settings.PREFIX))
   #run_process(_settings, "ln -t %s/Annotate/in/ -s %s/FindORFS/out/%s.scaffolds.faa"%(_settings.rundir,_settings.rundir,_settings.PREFIX))
