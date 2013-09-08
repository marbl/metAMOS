#!python

import os, sys, string, time, BaseHTTPServer, getopt, re, subprocess, webbrowser
from operator import itemgetter

from utils import *
from mapreads import MapReads 
sys.path.append(INITIAL_UTILS)
from ruffus import *

import generic

_readlibs = []
_skipsteps = []
_settings = Settings()

def init(reads, skipsteps):
   global _readlibs
   global _skipsteps
   global _usecontigs

   _readlibs = reads
   _skipsteps = skipsteps

def getAsmName(input_file_name):
   asmName = input_file_name.replace("%s/Assemble/out/"%(_settings.rundir), "")
   return asmName.replace(".contig.cvg", "")

@posttask(touch_file("%s/Logs/validate.ok"%(_settings.rundir)))
@merge(MapReads, ["%s/Validate/out/%s.lap"%(_settings.rundir, _settings.PREFIX)])
def Validate (input_file_names, output_file_name):
   bestScore = 0
   bestAssembler = ""
   bestAssembly = ""
   selectedAsm = open("%s/Validate/out/%s.asm.selected"%(_settings.rundir, _settings.PREFIX), 'w')
   lapfile = open("%s/Validate/out/%s.lap"%(_settings.rundir,_settings.PREFIX),'w')

   if "Validate" in _skipsteps or "validate" in _skipsteps:
      run_process(_settings, "touch %s/Logs/validate.skip"%(_settings.rundir), "Validate")
      bestAssembler = getAsmName(input_file_names[0])
      bestAssembly = input_file_names[0].replace(".contig.cvg", ".asm.contig") 
   elif len(input_file_names) == 1:
      run_process(_settings, "touch %s/Logs/validate.skip"%(_settings.rundir), "Validate")
      bestAssembler = getAsmName(input_file_names[0])
      bestAssembly = input_file_names[0].replace(".contig.cvg", ".asm.contig") 
   elif not os.path.exists("%s/bowtie2"%(_settings.BOWTIE2)) or not os.path.exists("%s/aligner/calc_prob.py"%(_settings.LAP)):
      run_process(_settings, "touch %s/Logs/validate.skip"%(_settings.rundir), "Validate")
      print "Warning! LAP is not available, cannot select best assembly, chosing first available: %s!"%(getAsmName(input_file_names[0]))
      bestAssembler = getAsmName(input_file_names[0])
      bestAssembly = input_file_names[0].replace(".contig.cvg", ".asm.contig") 
   else:
      # build string of files to use for validation
      os.environ["BT2_HOME"]=_settings.BOWTIE2

      pairedReads = ""
      unpairedReads = ""
      for lib in _readlibs:
         if lib.mated and pairedReads == "":
            pairedReads="-1 %s/Preprocess/out/lib%d.1.fastq -2 %s/Preprocess/out/lib%d.2.fastq -m %d -t %d -I %d -X %d"%(_settings.rundir, lib.id, _settings.rundir, lib.id, lib.mean, lib.stdev, (lib.mean-5*lib.stdev), (lib.mean+5*lib.stdev))
         elif not lib.paired and unpairedReads == "":
            unpaired=" -i %s/Preprocess/out/lib%d.fastq"%(_settings.rundir, lib.id)

      for input_file_name in input_file_names:
         assembler = getAsmName(input_file_name)
         assembly = input_file_name.replace(".contig.cvg", ".asm.contig")
         abundanceFile = ""
         if os.path.exists("%s/Assemble/out/%s.contig.cvg"%(_settings.rundir, assembler)):
            abundanceFile = "-n %s/Assemble/out/%s.contig.cvg"%(_settings.rundir, assembler)
         run_process(_settings, "python %s/aligner/calc_prob.py -p %d -a %s %s %s > %s/Validate/out/%s.prob"%(_settings.LAP, _settings.threads, assembly, abundanceFile, pairedReads if pairedReads != "" else unpairedReads, _settings.rundir, assembler), "Validate")
         score = getCommandOutput("python %s/aligner/sum_prob.py -i %s/Validate/out/%s.prob"%(_settings.LAP, _settings.rundir, assembler), True).split()[0]
         if _settings.VERBOSE:
            print "*** metAMOS assembler %s score %s."%(assembler, score)
         lapfile.write("%s\t%s\n"%(assembler, score))
         lapfile.flush()

         if bestAssembler == "" or float(score) > bestScore:
            bestScore = float(score)
            bestAssembler = assembler
            bestAssembly = assembly

   if _settings.VERBOSE:
      print "*** metAMOS assembler %s selected."%(bestAssembler)  
   run_process(_settings, "ln %s %s/Assemble/out/%s.asm.contig"%(bestAssembly, _settings.rundir, _settings.PREFIX), "Validate")
   run_process(_settings, "ln %s/Assemble/out/%s.asm.tigr %s/Assemble/out/%s.asm.tigr"%(_settings.rundir, bestAssembler, _settings.rundir, _settings.PREFIX), "Validate")
   run_process(_settings, "ln %s/Assemble/out/%s.contig.cnt %s/Assemble/out/%s.contig.cnt"%(_settings.rundir, bestAssembler, _settings.rundir, _settings.PREFIX), "Validate")
   run_process(_settings, "ln %s/Assemble/out/%s.contig.cvg %s/Assemble/out/%s.contig.cvg"%(_settings.rundir, bestAssembler, _settings.rundir, _settings.PREFIX), "Validate")
   if os.path.exists("%s/Assemble/out/%s.afg"%(_settings.rundir, bestAssembler)):
      run_process(_settings, "ln %s/Assemble/out/%s.afg %s/Assemble/out/%s.afg"%(_settings.rundir, bestAssembler, _settings.rundir, _settings.PREFIX), "Validate")

   for lib in _readlibs: 
      run_process(_settings, "ln %s/Assemble/out/%s.lib%d.badmades %s/Assemble/out/%s.lib%d.badmates"%(_settings.rundir, bestAssembler, lib.id, _settings.rundir, _settings.PREFIX, lib.id), "Validate")
      run_process(_settings, "ln %s/Assemble/out/%s.lib%d.hdr %s/Assemble/out/%s.lib%d.hdr"%(_settings.rundir, bestAssembler, lib.id, _settings.rundir, _settings.PREFIX, lib.id), "Validate")
      run_process(_settings, "ln %s/Assemble/out/%s.lib%d.mappedmates %s/Assemble/out/%s.lib%d.mappedmates"%(_settings.rundir, bestAssembler, lib.id, _settings.rundir, _settings.PREFIX, lib.id), "Validate")
      run_process(_settings, "ln %s/Assemble/out/%s.lib%d.mates_in_diff_contigs %s/Assemble/out/%s.lib%d.mates_in_diff_contigs"%(_settings.rundir, bestAssembler, lib.id, _settings.rundir, _settings.PREFIX, lib.id), "Validate")
      run_process(_settings, "ln %s/Assemble/out/%s.lib%dcontig.reads %s/Assemble/out/%s.lib%dcontig.reads"%(_settings.rundir, bestAssembler, lib.id, _settings.rundir, _settings.PREFIX, lib.id), "Validate")
      run_process(_settings, "ln %s/Assemble/out/%s.lib%d.unaligned.fasta %s/Assemble/out/lib%d.unaligned.fasta"%(_settings.rundir, bestAssembler, lib.id, _settings.rundir, lib.id), "Validate")
      if os.path.exists("%s/Assemble/out/%s.lib%d.unaligned.fastq"%(_settings.rundir, bestAssembler, lib.id)):
         run_process(_settings, "ln %s/Assemble/out/%s.lib%d.unaligned.fastq %s/Assemble/out/lib%d.unaligned.fastq"%(_settings.rundir, bestAssembler, lib.id, _settings.rundir, lib.id), "Validate")

   selectedAsm.write("%s"%(bestAssembler))
   selectedAsm.close()
   lapfile.close()

