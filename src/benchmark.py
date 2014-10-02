#!python

import os, sys, string, time, BaseHTTPServer, getopt, re, subprocess, webbrowser
from operator import itemgetter

from utils import *
from findorfs import FindORFS,setRunFast
from mapreads import getMeanSD
sys.path.append(INITIAL_UTILS)
from ruffus import *

import generic

_readlibs = []
_skipsteps = []
_settings = Settings()


def init(reads, skipsteps,rulers):
   global _readlibs
   global _skipsteps

   _rulers = rulers
   _readlibs = reads
   _skipsteps = skipsteps

def getRulers:
   #emulate get workflows? in runPipeline, not here
   pass

def getClassifierName(classifier):
   (clsname, citation) = getProgramCitations(_settings, classifier.lower())
   if citation == "UNKNOWN":
      asmName = assembler.upper()

   return clsname


#
#generate test sets for a given process
#outline
#need to skip scaffold, assemble
#this should be ran as follows:
#initPipeline -w measur_sensor
# --this defines input data and test cases?
# --initPipeline generates sim data, wgsim, multimode 
# for each sim dataset 9 total
#  -- runPipeline -c kraken
#  -- runPipeline -c lmat

#runPipeline -c kraken,lmat,siann
# --this would be better as multi-metamos, multiple metamos runs with a single classifier?
#then, only preprocess,classification,benchmark,postprocess 
#OR, preprocess,classification,postprocess
#    preprocess,classification,postprocess
#    then benchmark?
#TBD
#1. define test cases, read in from spec file?
#2. define methods to be tested (user input?)
#3. run all methods
#4. compare output to expected output
#5. create summary table of result

@posttask(touch_file("%s/Logs/benchmark.ok"%(_settings.rundir)))
@merge(FindORFS, ["%s/Logs/benchmark.ok"%(_settings.rundir)])
def Benchmark (input_file_names, output_file_name):
   benchmarkScores = dict()

   if "Benchmark" in _skipsteps or "benchmark" in _skipsteps:
      run_process(_settings, "touch %s/Logs/benchmark.skip"%(_settings.rundir), "Benchmark")
      if len(input_file_names) > 0:
         bestAssembler = getAsmName(input_file_names[0])
         bestAssembly = "%s/Assemble/out/%s.asm.contig"%(_settings.rundir, bestAssembler)
      else:
         bestAssembler = "none"
         bestAssembly = "%s/Assemble/out/%s.asm.contig"%(_settings.rundir, bestAssembler)
         run_process(_settings, "touch %s/Assemble/out/%s.asm.contig"%(_settings.rundir, bestAssembler), "Benchmark")
         run_process(_settings, "touch %s/Assemble/out/%s.contig.cvg"%(_settings.rundir, bestAssembler), "Benchmark")
         for lib in _readlibs:
            run_process(_settings, "touch %s/Assemble/out/%s.lib%dcontig.reads"%(_settings.rundir, bestAssembler, lib.id), "Benchmark")

   #how to match step to ruler type? supported steps for each ruler type? 
   #rulers i have in mind: sensor, scale, survey, screen, signature

   #for all classifiers, measure with all rulers, store in ./Benchmark/out/results.txt
   #then, for multi-runs, collate */Benchmark/out/results.txt into fig/table

   #format for output file: 
   #Ruler Classifier Measurement#1 Measurement#2 .. Measurement#N \n
   #Sensor Kraken 0.9 0.7 .. 0.8 \n
   #Sensor FCP 0.9 0.7 .. 0.8 \n
   #Sensor Phylosift 0.9 0.7 .. 0.8 \n
   #Survey Kraken 0.9 0.7 .. 0.8 \n
   #Survey FCP 0.9 0.7 .. 0.8 \n
   #Survey Phylosift 0.9 0.7 .. 0.8 \n

   for ruler in rulers:
       for c in classifiers:
           #truth table part of ruler? defines a dataset and method to evaluate
           #returns recall & precision
 
           measurements = ruler.run("%s/ClassifyReads/out/results.txt")
           print c.name #ie kraken
           for m in measurements:
               print m.name #ie precision
               print m.value #ie 0.9

 


     


