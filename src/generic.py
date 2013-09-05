#!python

import os, sys, string, time, BaseHTTPServer, getopt, re, subprocess, webbrowser
from datetime import date
from datetime import time
from datetime import datetime
from operator import itemgetter
from utils import *
from task import JobSignalledBreak

_readlibs = []
_skipsteps = []
_settings = ""
def init(skipsteps, readlibs):
   global _skipsteps
   global _readlibs
   global _settings
   _skipsteps = skipsteps
   _readlibs = readlibs
   _settings = Settings()

# currently lots of special keywords, get away from this
# current keywords
#   - MEM - max memory limit
#   - LIB - library identifier
#   - INPUT - replace with input to the program
#   - MACHINE - replaced with linux-64bit, etc
#   - FIRST - replaced with left mates in mated read or interleaved or unpaired reads otherwise
#   - SECOND - replaced with right mates, on in paired non-interleaved libs
#   - THREADS - replaced with thread parameter specified
#   - KMER - the kmer requested
#   - OFFSET - the phred offset (33/64)
#   - PREFIX - the desied output name
#   - DB - the location of dbs
class GenericProgram:
   # for now we only support new assemblers
   stepName = STEP_NAMES.ASSEMBLE 
   inputType = INPUT_TYPE.FASTQ 
   name = ""
   output = ""
   location = ""
   threads = ""
   paired = ""
   paired_interleaved = ""
   unpaired = ""
   allowPartition = False
   commandList = []
   isValid = False

   def __init__(self, name = "", step=STEP_NAMES.ASSEMBLE):
      self.name = name.lower()
      self.stepName = step
      self.output = ""
      self.location = ""
      self.threads = ""
      self.paired = ""
      self.paired_interleaved = ""
      self.unpaired = ""
      self.allowPartition = False 
      self.commandList = []
      self.isValid = False 

   def read(self):
      # populate from file
      self.isValid = True

      cfg = getProgramParams(_settings.METAMOS_UTILS, "%s.spec"%(self.name), "CONFIG", "\n", "#", " ")
      for line in cfg.split("\n"):
         (defn, sp, value) = line.partition(' ')         
         value = value.strip()
         if (defn == "name"):
            if (not value.lower() == self.name):
               print "Error: name %s does not match expected %s\n"%(value, self.name)
               self.isValid = False
            self.name = value
         elif (defn == "input"):
            try:
               self.inputType = INPUT_TYPE.mapping[value.upper()]
            except KeyError:
               print "Error: unknown input type specfied %s. Must be one of %s\n"%(value.upper(), INPUT_TYPE.mapping.keys())
               self.isValid = False
         elif (defn == "output"):
            self.output = value
         elif (defn == "location"):
           # get the location path
            self.location = value
            self.location = self.location.replace("[MACHINE]", "%s-%s"%(_settings.OSTYPE, _settings.MACHINETYPE))
            if not os.path.isabs(self.location):
               self.location = os.path.abspath("%s%s%s"%(_settings.METAMOS_UTILS, os.sep, self.location))

         elif (defn == "threads"):
            self.threads = value
         elif (defn == "paired"):
            self.paired = value  
         elif (defn == "paired_interleaved"):
            self.paired_interleaved = value
         elif (defn == "unpaired"):
            self.unpaired = value
         elif (defn == "partition"):
            self.allowPartition = str2bool(value.strip())
         elif (defn == "commands"):
            self.commandList.append(value)

      # check integrity
      if (len(self.commandList) == 0):
         print "Error no commands to run provided\n"
         self.isValid = False
      if self.inputType == INPUT_TYPE.FASTA or self.inputType == INPUT_TYPE.FASTQ:
         if ((self.paired == "" and self.paired_interleaved == "") or self.unpaired == ""):
            print "Error: must specify handler for both paired and unpaired sequences"
            self.isValid = False

      for command in self.commandList:
         (commandName, sp, junk) = command.partition(' ')
         theCommand = "%s%s%s"%(self.location, os.sep, commandName)
         if (not os.path.exists(theCommand)):
            # try to find in path
            self.location = getFromPath(commandName, self.name)

   def execute(self):
      print "Starting execute"

      if not self.isValid:
         print "Error: cannot execute %s. Invalid configuration provided\n"%(self.name)
         raise(JobSignalledBreak)

      avram = getAvailableMemory(_settings)
      params = []
      outputs = []
      listOfInput = ""
      offset = ""

      if self.inputType == INPUT_TYPE.FASTA or self.inputType == INPUT_TYPE.FASTQ:
         params.append("")
         suffix = ""
         if self.inputType == INPUT_TYPE.FASTA:
            suffix = "fasta"
         else:
            suffix = "fastq"

         for lib in _readlibs:
            for read in lib.reads:
               if offset == "":
                  offset = read.qformat
               elif not offset == read.qformat:
                  print "Error: inconsistent PHRED offsets in libraries. Previous library had %s and current library %d has %s\n"%(offset, lib.id, read.qformat)
                  raise(JobSignalledBreak)
               
            if lib.mated:
               if self.paired == "":
                  params[0] = params[0] + " " + self.paired_interleaved.replace("[LIB]", "%d"%(lib.id)).replace("[FIRST]", "%s/Preprocess/out/lib%d.%s"%(_settings.rundir, lib.id, suffix))
               else: 
                  params[0] = params[0] + " " + self.paired.replace("[LIB]", "%d"%(lib.id)).replace("[FIRST]", "%s/Preprocess/out/lib%d.1.%s"%(_settings.rundir, lib.id, suffix)).replace("[SECOND]", "%s/Preprocess/out/lib%d.2.%s"%(_settings.rundir, lib.id, suffix))
            elif not lib.mated:
                  params[0] = params[0] + " " + self.unpaired.replace("[LIB]", "%d"%(lib.id)).replace("[FIRST]", "%s/Preprocess/out/lib%d.%s"%(_settings.rundir, lib.id, suffix))
      elif self.inputType == INPUT_TYPE.CONTIGS:
         params.append("%s/Annotate/in/%s.asm.contig"%(_settings.rundir, _settings.PREFIX))
         listOfInput += ":%s/Annotate/in/%s.asm.contig"%(_settings.rundir, _settings.PREFIX)
         if _settings.annotate_unmapped:
            outputs.append("%s/Annotate/out/%s.intermediate.ctg"%(_settings.rundir,_settings.PREFIX))
            for lib in _readlibs:
               listOfInput += ":%s/Assemble/out/lib%d.unaligned.fasta"%(_settings.rundir, lib.id)
               run_process(_settings, "ln %s/Assemble/out/lib%d.unaligned.fasta %s/Annotate/in/lib%d.unaligned.fasta"%(_settings.rundir, lib.id, _settings.rundir, lib.id), STEP_NAMES.reverse_mapping[self.stepName].title())
               params.append("%s/Annotate/in/lib%d.unaligned.fasta"%(_settings.rundir, lib.id))
               outputs.append("%s/Annotate/out/%s.intermediate.lib%d"%(_settings.rundir, _settings.PREFIX, lib.id))
      elif self.inputType == INPUT_TYPE.ORF_FA:
         params.append("%s/Annotate/in/%s.fna"%(_settings.rundir, _settings.PREFIX))
         listOfInput += ":%s/Annotate/in/%s.fna"%(_settings.rundir, _settings.PREFIX)
      elif self.inputType == INPUT_TYPE.ORF_AA:
         params.append("%s/Annotate/in/%s.faa"%(_settings.rundir, _settings.PREFIX))
         listOfInput += ":%s/Annotate/in/%s.faa"%(_settings.rundir, _settings.PREFIX)
         print "Parameters are %s\n"%(params)

      # get the thread parameters
      # when thread params not available, should partition the data
      threadParams = ""
      if not self.threads == "":
         threadParams =  "%s %d"%(self.threads,_settings.threads)

      i = 0
      for param in params:
         for command in self.commandList:
            (commandName, sp, junk) = command.partition(' ')
            theCommand = "%s%s%s"%(self.location, os.sep, commandName)
            if (not os.path.exists(theCommand)): 
               print "Error: requested to run %s but not available in specified location %s. Please check your specification and try again"%(commandName, self.location)
               raise(JobSignalledBreak)

            param = param + " " + getProgramParams(_settings.METAMOS_UTILS, "%s.spec"%(self.name.lower()), commandName, "-").replace("[KMER]", "%d"%(_settings.kmer))
            command = self.location + os.sep + command.replace("[INPUT]", param).replace("[DB]", _settings.DB_DIR).replace("[MEM]", "%d"%(avram)).replace("[THREADS]", threadParams).replace("[OFFSET]", "33" if offset.lower() == "sanger" else "64").replace("[OUTPUT]", self.output.replace("[PREFIX]", outputs[i]) if len(outputs) > i else "%s%s%s%sout%s%s"%(_settings.rundir, os.sep, STEP_NAMES.reverse_mapping[self.stepName].title(), os.sep, os.sep, self.output.replace("[PREFIX]", _settings.PREFIX)))
            run_process(_settings, command, STEP_NAMES.reverse_mapping[self.stepName].title())
         i+=1

      i = 0
      for output in outputs:
         if i == 0:
            run_process(_settings, "unlink %s/%s/out/%s"%(_settings.rundir, STEP_NAMES.reverse_mapping[self.stepName].title(), self.output.replace("[PREFIX]", _settings.PREFIX)), STEP_NAMES.reverse_mapping[self.stepName].title())
         run_process(_settings, "cat %s >> %s/%s/out/%s"%(self.output.replace("[PREFIX]", output), _settings.rundir, STEP_NAMES.reverse_mapping[self.stepName].title(), self.output.replace("[PREFIX]", _settings.PREFIX)), STEP_NAMES.reverse_mapping[self.stepName].title())
         i = 1

      # symlink results
      symlinkCmd = "ln %s/%s/out/%s %s/%s/out/%s%s"%(_settings.rundir, STEP_NAMES.reverse_mapping[self.stepName].title(), self.output.replace("[PREFIX]", _settings.PREFIX), _settings.rundir, STEP_NAMES.reverse_mapping[self.stepName].title(), _settings.PREFIX, STEP_OUTPUTS.reverse_mapping[self.stepName])
      run_process(_settings, symlinkCmd, STEP_NAMES.reverse_mapping[self.stepName].title())
      return listOfInput

def getSupportedList(path, step):
   if (step < STEP_NAMES.ASSEMBLE or step > STEP_NAMES.ANNOTATE):
      return []
   programs = getProgramParams(path, "%s.generic"%(STEP_NAMES.reverse_mapping[step]))
   return programs.strip().split()
   
def checkIfExists(step, programName):
   if (step < STEP_NAMES.ASSEMBLE or step > STEP_NAMES.ANNOTATE):
      return False

   programs = getSupportedList(_settings.METAMOS_UTILS, step)
   return programName.lower() in programs

def getLocation(step, programName):
   if (step < STEP_NAMES.ASSEMBLE or step > STEP_NAMES.ANNOTATE):
      return "UNKNOWN"
   if not checkIfExists(step, programName):
      return "UNKNOWN" 
   dispatch = GenericProgram(programName, step)
   dispatch.read()

   return dispatch.location

def getName(step, programName):
   if (step < STEP_NAMES.ASSEMBLE or step > STEP_NAMES.ANNOTATE):
      return "UNKNOWN"
   if not checkIfExists(step, programName):
      return "UNKNOWN"
   dispatch = GenericProgram(programName, step)
   dispatch.read()

   return dispatch.name

def execute(step, programName):
   if (step < STEP_NAMES.ASSEMBLE or step > STEP_NAMES.ANNOTATE):
      print "Error: unsupported step #%d\n"%(step)
      raise(JobSignalledBreak)

   if not checkIfExists(step, programName):
      print "Error: attempted to execute program %s which is not supported\n"%(programName)
      raise(JobSignalledBreak)

   dispatch = GenericProgram(programName, step)
   dispatch.read()
   return dispatch.execute()
