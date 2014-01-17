#!python

import os, sys, string, time, BaseHTTPServer, getopt, re, subprocess, webbrowser
from datetime import date
from datetime import time
from datetime import datetime
from operator import itemgetter
from utils import *
from task import JobSignalledBreak

LIBRARY_TYPES = enum("PAIRED", "MATED")
TECHNOLOGY_TYPES = enum("SOLEXA", "SANGER", "PACBIO", "454")
SYSTEM_COMMANDS = [ "mkdir", "mv", "bash", "ln", "rm", "cp", "ls" ]

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
#   - ORIENTATION - replaced with either innie or outtie
#   - ORIENTATION_FIGURE - replaced with ---> <--- or <--- ---> for pe and mp, respectively
#   - MEAN - replaced with mean
#   - SD - replaced with standard dev
#   - MIN - replaced with min
#   - MAX - replaced with max
#   - THREADS - replaced with thread parameter specified
#   - KMER - the kmer requested
#   - OFFSET - the phred offset (33/64)
#   - PREFIX - the desied output name
#   - DB - the location of dbs
#   - RUNDIR - the location where the program is running
#   - LOCATION - the location where the program executable lives
#   - TECHNOLOGY - the type of sequencing data
class GenericProgram:
   # for now we only support new assemblers
   stepName = STEP_NAMES.ASSEMBLE 
   inputType = INPUT_TYPE.FASTQ 
   name = ""
   output = ""
   scaffoldOutput = ""
   backupOutput = ""
   location = ""
   threads = ""
   threadsSupportedOn = []
   paired = ""
   paired_interleaved = ""
   mated = ""
   mated_interleaved = ""
   unpaired = ""
   requiredLibs = dict() 
   maxLibs = 0
   maxLibIDLen = 0
   config = None
   allowPartition = False
   commandList = []
   isValid = False
   technologyParams = dict()

   def __init__(self, name = "", step=STEP_NAMES.ASSEMBLE):
      self.name = name.lower()
      self.stepName = step
      self.output = ""
      self.scaffoldOutput = ""
      self.backupOutput = ""
      self.location = ""
      self.threads = ""
      self.threadsSupportedOn = []
      self.paired = ""
      self.paired_interleaved = ""
      self.mated = ""
      self.mated_interleaved = ""
      self.unpaired = ""
      self.requiredLibs = dict()
      self.maxLibs = 1
      self.allowPartition = False 
      self.commandList = []
      self.isValid = False 
      self.technologyParams = dict()

   def getCommandBin(self, command):
      # skip any leading special flags
      reading = False
      count = 0
      for c in command:
        if c == "[":
           reading=True

        if c != " " and not reading:
           break
        count += 1

        if c == "]":
           reading=False
      (commandName, sp, junk) = command[count:].partition(' ')
      if commandName in SYSTEM_COMMANDS:
         theCommand = getFromPath(commandName, self.name, False) + os.sep + commandName
      else:
         theCommand = "%s%s%s"%(self.location, os.sep, commandName)
      return (commandName, theCommand)

   def read(self):
      # populate from file
      self.isValid = True

      cfg = getProgramParams(_settings.METAMOS_UTILS, "%s.spec"%(self.name), "CONFIG", "\n", "#", " ")
      for line in cfg.split("\n"):
         (defn, sp, value) = line.partition(' ')         
         value = value.strip()
         if (defn == "name"):
            if (not value.lower() == self.name.lower()):
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
         elif (defn == "scaffoldOutput"):
            self.scaffoldOutput = value
         elif (defn == "backupOutput"):
            self.backupOutput = value
         elif (defn == "location"):
           # get the location path
            self.location = value
            self.location = self.location.replace("[MACHINE]", "%s-%s"%(_settings.OSTYPE, _settings.MACHINETYPE)).replace("[MACHINE-CA]", "%s-%s"%(_settings.OSTYPE, _settings.MACHINETYPE.replace("x86_64", "amd64")))
            if not os.path.isabs(self.location):
               self.location = os.path.abspath("%s%s%s"%(_settings.METAMOS_UTILS, os.sep, self.location))
         elif (defn == "threads"):
            self.threads = value
         elif (defn == "threadsSupport"):
            self.threadsSupportedOn = value.strip().split()
         elif (defn == "paired"):
            self.paired = value  
         elif (defn == "paired_interleaved"):
            self.paired_interleaved = value
         elif (defn == "mated"):
            self.mated = value
         elif (defn == "mated_interleaved"):
            self.mated_interleaved = value
         elif (defn == "unpaired"):
            self.unpaired = value
         elif (defn == "maxlibs"):
            self.maxLibs = int(value)
         elif (defn == "maxlibid"):
            self.maxLibIDLen = int(value)
         elif (defn == "config"):
            self.config = value
         elif (defn == "required"):
            splitLine = value.strip().split(",")
            for val in splitLine:
               val = val.replace("[", "").replace("]", "")
               if val.upper() in LIBRARY_TYPES.mapping:
                  self.requiredLibs[LIBRARY_TYPES.mapping[val.upper()]] = True
               else:
                 print "Warning: unknown library requirement %s specified"%(val)
         elif (defn == "partition"):
            self.allowPartition = str2bool(value.strip())
         elif (defn == "commands"):
            self.commandList.extend(value.split("&&"))
         elif(defn.upper() in TECHNOLOGY_TYPES.mapping):
            self.technologyParams[TECHNOLOGY_TYPES.mapping[defn.upper()]] = value

      # check integrity
      if (len(self.commandList) == 0):
         print "Error no commands to run provided\n"
         self.isValid = False
      if self.inputType == INPUT_TYPE.FASTA or self.inputType == INPUT_TYPE.FASTQ:
         if ((self.paired == "" and self.paired_interleaved == "") and self.unpaired == ""):
            print "Error: must at least specify handler for paired or unpaired sequences"
            self.isValid = False
         for reqLib in self.requiredLibs:
            if self.requiredLibs[reqLib] == True:
               if reqLib == LIBRARY_TYPES.PAIRED:
                  if (self.paired == "" and self.paired_interleaved == ""):
                     print "Error: paired library is required but no handler specified"
                     self.isValid = False
               elif reqLib == LIBRARY_TYPES.MATED:
                  if (self.mated == "" and self.mated_interleaved == ""):
                     print "Error: mated library is required but no handler specified"
                     self.isValid = False
               else:
                  print "Error: unknown library requirement specified"
      for command in self.commandList:
         (commandName, theCommand) = self.getCommandBin(command)
         if (not os.path.exists(theCommand)):
            # try to find in path
            self.location = getFromPath(commandName, self.name, False)

   def getThreadParams(self):
      threadParams = ""
      if (len(self.threadsSupportedOn) == 0 or _settings.OSTYPE.upper() in self.threadsSupportedOn):
         if self.threads != "":
            if self.threads.endswith("="):
               threadParams =  "%s%d"%(self.threads,_settings.threads)
            else:
               threadParams =  "%s %d"%(self.threads,_settings.threads)
         elif self.config != None and len(self.config) != 0:
            threadParams = "%d"%(_settings.threads)

      return threadParams.strip()

   def getLibInput(self, separator=" "):
         techParams = ""
         techAdded = []

         suffix = ""
         if self.inputType == INPUT_TYPE.FASTA:
            suffix = "fasta"
         else:
            suffix = "fastq"

         input = ""
         offset = ""
         havePE = haveMP = False
         pelist = ""
         mplist = ""

         if self.maxLibs == 0:
            self.maxLibs = len(_readlibs)
         libCount = 0

         for lib in _readlibs:
            if libCount > self.maxLibs:
               print "Warning: selected assembler %s only supports a maximum of %d libs, not using all libraries"%(self.name, self.maxLibs)
               break

            for read in lib.reads:
               if offset == "":
                  offset = read.qformat
               elif not offset == read.qformat:
                  print "Error: inconsistent PHRED offsets in libraries. Previous library had %s and current library %d has %s\n"%(offset, lib.id, read.qformat)
                  raise(JobSignalledBreak)

            orientation = "outtie"
            orientationFig = "<--- --->"
            if lib.innie:
               orientation = "innie"
               orientationFig = "---> <---"
            if lib.format == "sff":
               pyroType = TECHNOLOGY_TYPES.mapping["454"]
               technology = TECHNOLOGY_TYPES.reverse_mapping[pyroType].lower()
               if pyroType in self.technologyParams and technology not in techAdded:
                  techParams = "%s_SETTINGS %s"%(technology.upper(), self.technologyParams[pyroType])
                  techAdded.append(technology)
            else:
               technology = TECHNOLOGY_TYPES.reverse_mapping[TECHNOLOGY_TYPES.SOLEXA].lower()
               if TECHNOLOGY_TYPES.SOLEXA in self.technologyParams and technology not in techAdded:
                  techParams = "%s_SETTINGS %s"%(technology.upper(), self.technologyParams[TECHNOLOGY_TYPES.SOLEXA])
                  techAdded.append(technology)

            # eventually check others

            if lib.mated:
               if lib.innie:
                  havePE = True
                  pelist = "%s lib%s"%(pelist, lib.id)
                  if self.paired == "":
                     input = input + separator + self.paired_interleaved.replace("[LIB]", "%d"%(lib.id)).replace("[FIRST]", "%s/Preprocess/out/lib%d.%s"%(_settings.rundir, lib.id, suffix)).replace("[MEAN]", "%d"%(lib.mean)).replace("[SD]", "%d"%(lib.stdev)).replace("[MIN]", "%d"%(lib.mmin)).replace("[MAX]", "%d"%(lib.mmax)).replace("[TECHNOLOGY]", technology).replace("[ORIENTATION]", orientation).replace("[ORIENTATION_FIGURE]", orientationFig)
                  else:
                     input = input + separator + self.paired.replace("[LIB]", "%d"%(lib.id)).replace("[FIRST]", "%s/Preprocess/out/lib%d.1.%s"%(_settings.rundir, lib.id, suffix)).replace("[SECOND]", "%s/Preprocess/out/lib%d.2.%s"%(_settings.rundir, lib.id, suffix)).replace("[MEAN]", "%d"%(lib.mean)).replace("[SD]", "%d"%(lib.stdev)).replace("[MIN]", "%d"%(lib.mmin)).replace("[MAX]", "%d"%(lib.mmax)).replace("[TECHNOLOGY]", technology).replace("[ORIENTATION]", orientation).replace("[ORIENTATION_FIGURE]", orientationFig)
               else:
                  haveMP = True
                  mplist = "%s lib%s"%(mplist, lib.id) 
                  if self.mated == "":
                     input = input + separator + self.mated_interleaved.replace("[LIB]", "%d"%(lib.id)).replace("[FIRST]", "%s/Preprocess/out/lib%d.%s"%(_settings.rundir, lib.id, suffix)).replace("[MEAN]", "%d"%(lib.mean)).replace("[SD]", "%d"%(lib.stdev)).replace("[MIN]", "%d"%(lib.mmin)).replace("[MAX]", "%d"%(lib.mmax)).replace("[TECHNOLOGY]", technology).replace("[ORIENTATION]", orientation).replace("[ORIENTATION_FIGURE]", orientationFig)
                  else:
                     input = input + separator + self.mated.replace("[LIB]", "%d"%(lib.id)).replace("[FIRST]", "%s/Preprocess/out/lib%d.1.%s"%(_settings.rundir, lib.id, suffix)).replace("[SECOND]", "%s/Preprocess/out/lib%d.2.%s"%(_settings.rundir, lib.id, suffix)).replace("[MEAN]", "%d"%(lib.mean)).replace("[SD]", "%d"%(lib.stdev)).replace("[MIN]", "%d"%(lib.mmin)).replace("[MAX]", "%d"%(lib.mmax)).replace("[TECHNOLOGY]", technology).replace("[ORIENTATION]", orientation).replace("[ORIENTATION_FIGURE]", orientationFig)
            elif not lib.mated:
                  input = input + separator + self.unpaired.replace("[LIB]", "%d"%(lib.id)).replace("[FIRST]", "%s/Preprocess/out/lib%d.%s"%(_settings.rundir, lib.id, suffix)).replace("[TECHNOLOGY]", technology)

            return (havePE, haveMP, offset, input, pelist.strip(), mplist.strip(), techParams.strip())

   def processConfig(self, libs, offset, pelist, mplist, techParams, avram):
      configFile = ""

      if self.config != None and len(self.config) != 0:
         try:
            configFile = "%s/%s/out/%s.config"%(_settings.rundir, STEP_NAMES.reverse_mapping[self.stepName].title(), _settings.PREFIX)
            template = open("%s/%s"%(_settings.METAMOS_UTILS, self.config), 'r')
            spec = open(configFile, 'w')

            for line in template.xreadlines():
               line = line.strip()
               if "[INPUT]" in line:
                  for l in libs.strip().split("\\n"):
                     spec.write(l.strip() + "\n")
               else:
                  line = line.replace("[TECHNOLOGY_PARAMETERS]", techParams).replace("[MACHINE]", "%s-%s"%(_settings.OSTYPE, _settings.MACHINETYPE)).replace("[MACHINE-CA]", "%s-%s"%(_settings.OSTYPE, _settings.MACHINETYPE.replace("x86_64", "amd64"))).replace("[MPI]", _settings.MPI).replace("[PELIST]", pelist).replace("[MPLIST]", mplist).replace("[PREFIX]", _settings.PREFIX).replace("[DB]", _settings.DB_DIR).replace("[MEM]", "%d"%(avram)).replace("[THREADS]",self.getThreadParams()).replace("[OFFSET]", "33" if offset.lower() == "sanger" else "64").replace("[OUTPUT]", "%s%s%s%sout%s%s"%(_settings.rundir, os.sep, STEP_NAMES.reverse_mapping[self.stepName].title(), os.sep, os.sep, self.output.replace("[PREFIX]", _settings.PREFIX))).replace("[RUNDIR]", "%s%s%s%sout"%(_settings.rundir, os.sep, STEP_NAMES.reverse_mapping[self.stepName].title(), os.sep)).replace("[KMER]", "%s"%(_settings.kmer)).replace("[LOCATION]", "%s%s"%(self.location,os.sep))
                  spec.write(line + "\n")
            template.close()
            spec.close()
         except IOError:
            print "Error: could not load config file for %s."%(self.name)
            raise(JobSignalledBreak)
      return configFile

   def execute(self):
      if not self.isValid:
         print "Error: cannot execute %s. Invalid configuration provided\n"%(self.name)
         raise(JobSignalledBreak)

      avram = getAvailableMemory(_settings)
      params = []
      outputs = []
      techParams = ""
      listOfInput = ""
      offset = pelist = mplist = ""

      if self.inputType == INPUT_TYPE.FASTA or self.inputType == INPUT_TYPE.FASTQ:
         if self.config != None and len(self.config) != 0:
            (havePE, haveMP, offset, input, pelist, mplist, techParams) = self.getLibInput("\n")
            configFile = self.processConfig(input, offset, pelist, mplist, techParams, avram)
            params.append("%s"%(configFile))
         else:
            (havePE, haveMP, offset, input, pelist, mplist, techParams) = self.getLibInput()
            params.append(input)

         # check for required libs
         for reqLib in self.requiredLibs.keys():
            if self.requiredLibs[reqLib] == True:
               if reqLib == LIBRARY_TYPES.PAIRED:
                  if havePE == False:
                     print "Error: cannot run %s. Required PE library but not specified"%(self.name)
                     return
               elif reqLib == LIBRARY_TYPES.MATED:
                  if haveMP == False:
                     print "Error: cannot run %s. Required PE library but not specified"%(self.name)
                     return
               else:
                  print "Error: unknown library type %s"%(reqLib)
                  raise(JobSignalledBreak)
      elif self.inputType == INPUT_TYPE.CONTIGS:
         params.append("%s/Annotate/in/%s.asm.contig"%(_settings.rundir, _settings.PREFIX))
         listOfInput += ":%s/Annotate/in/%s.asm.contig"%(_settings.rundir, _settings.PREFIX)
         if _settings.annotate_unmapped:
            outputs.append("%s/Annotate/out/%s.intermediate.ctg"%(_settings.rundir,_settings.PREFIX))
            for lib in _readlibs:
               if os.path.exists("%s/Assemble/out/lib%d.unaligned.fasta"%(_settings.rundir, lib.id)) and os.path.getsize("%s/Assemble/out/lib%d.unaligned.fasta"%(_settings.rundir, lib.id)) > 0: 
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

      # lots of tools expect to be able to get themselves from path so add run dir to path
      oldLDPath = ""
      oldPath = os.environ["PATH"]
      if "LD_LIBRARY_PATH" in os.environ:
         oldLDPath = os.environ["LD_LIBRARY_PATH"]
      os.environ["PATH"] = self.location + os.pathsep + oldPath
      os.environ["LD_LIBRARY_PATH"] = self.location + os.pathsep + oldLDPath

      i = 0
      for param in params:
         for command in self.commandList:
            (commandName, theCommand) = self.getCommandBin(command)
            if (not os.path.exists(theCommand)): 
               print "Error: requested to run %s (%s) but not available in specified location %s. Please check your specification and try again"%(self.name, commandName, self.location)
               return
            if ("[MPI]" in command and not os.path.exists(_settings.MPI)):
               print "Error: requested to run %s but required MPI is not available"%(commandName)
               command = command.replace("[MPI]", "").replace("[THREADS]", "")

            if ("[FIRST]" in command or "[SECOND]" in command):
               if self.maxLibs != 1:
                  print "Error: program %s requested library file as part of command but more than one library is available"%(self.name)
                  return
               else:
                  suffix = ""
                  if self.inputType == INPUT_TYPE.FASTA:
                     suffix = "fasta"
                  else:
                     suffix = "fastq"
                  command = command.replace("[FIRST]", "%s/Preprocess/out/lib%d.1.%s"%(_settings.rundir, 1, suffix)).replace("[SECOND]", "%s/Preprocess/out/lib%d.2.%s"%(_settings.rundir, 1, suffix))

            param = param + " " + getProgramParams(_settings.METAMOS_UTILS, "%s.spec"%(self.name.lower()), commandName, "-").replace("[KMER]", "%s"%(_settings.kmer))
            command = command.replace("[TECHNOLOGY_PARAMETERS]", techParams).replace(commandName, theCommand, 1).replace("[MPI]", _settings.MPI).replace("[INPUT]", param).replace("[PELIST]", pelist).replace("[MPLIST]", mplist).replace("[PREFIX]", _settings.PREFIX).replace("[DB]", _settings.DB_DIR).replace("[MEM]", "%d"%(avram)).replace("[THREADS]", self.getThreadParams()).replace("[OFFSET]", "33" if offset.lower() == "sanger" else "64").replace("[OUTPUT]", self.output.replace("[PREFIX]", outputs[i]) if len(outputs) > i else "%s%s%s%sout%s%s"%(_settings.rundir, os.sep, STEP_NAMES.reverse_mapping[self.stepName].title(), os.sep, os.sep, self.output.replace("[PREFIX]", _settings.PREFIX))).replace("[RUNDIR]", "%s%s%s%sout"%(_settings.rundir, os.sep, STEP_NAMES.reverse_mapping[self.stepName].title(), os.sep)).replace("[KMER]", "%s"%(_settings.kmer))
            run_process(_settings, command, STEP_NAMES.reverse_mapping[self.stepName].title())
         i+=1

      i = 0
      for output in outputs:
         if i == 0:
            run_process(_settings, "unlink %s/%s/out/%s"%(_settings.rundir, STEP_NAMES.reverse_mapping[self.stepName].title(), self.output.replace("[PREFIX]", _settings.PREFIX)), STEP_NAMES.reverse_mapping[self.stepName].title())
         run_process(_settings, "cat %s >> %s/%s/out/%s"%(self.output.replace("[PREFIX]", output), _settings.rundir, STEP_NAMES.reverse_mapping[self.stepName].title(), self.output.replace("[PREFIX]", _settings.PREFIX)), STEP_NAMES.reverse_mapping[self.stepName].title())
         i = 1

      # symlink results
      programOut = "%s/%s/out/%s"%(_settings.rundir, STEP_NAMES.reverse_mapping[self.stepName].title(), self.output.replace("[PREFIX]", _settings.PREFIX))
      stepOut = "%s/%s/out/%s%s"%(_settings.rundir, STEP_NAMES.reverse_mapping[self.stepName].title(), _settings.PREFIX, STEP_OUTPUTS.reverse_mapping[self.stepName])

      if programOut != stepOut:
         symlinkCmd = "ln %s %s"%(programOut, stepOut) 
         run_process(_settings, "unlink %s"%(stepOut), STEP_NAMES.reverse_mapping[self.stepName])
         if not os.path.exists(programOut) and self.backupOutput != "" and self.stepName == STEP_NAMES.ASSEMBLE:
            unpairedOut = "%s/%s/out/%s"%(_settings.rundir, STEP_NAMES.reverse_mapping[self.stepName].title(), self.backupOutput.replace("[PREFIX]", _settings.PREFIX))
            symlinkCmd = "ln %s %s"%(unpairedOut, stepOut)

         run_process(_settings, symlinkCmd, STEP_NAMES.reverse_mapping[self.stepName].title())

      if _settings.doscaffolding and self.scaffoldOutput != "" and self.stepName == STEP_NAMES.ASSEMBLE:
         scfOut = "%s/%s/out/%s"%(_settings.rundir, STEP_NAMES.reverse_mapping[self.stepName].title(), self.scaffoldOutput.replace("[PREFIX]", _settings.PREFIX))
         if os.path.exists(scfOut):
            run_process(_settings, "ln %s %s/%s/out/%s%s"%(scfOut, _settings.rundir, STEP_NAMES.reverse_mapping[self.stepName].title(), _settings.PREFIX, STEP_OUTPUTS.reverse_mapping[STEP_NAMES.mapping["SCAFFOLD"]]), STEP_NAMES.reverse_mapping[self.stepName].title()) 
         else:
            run_process(_settings, "ln %s %s/%s/out/%s%s"%(programOut, _settings.rundir, STEP_NAMES.reverse_mapping[self.stepName].title(), _settings.PREFIX, STEP_OUTPUTS.reverse_mapping[STEP_NAMES.mapping["SCAFFOLD"]]), STEP_NAMES.reverse_mapping[self.stepName].title())    

      # restore path
      os.environ["PATH"] = oldPath
      os.environ["LD_LIBRARY_PATH"] = oldLDPath
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

def execute(step, programName, settings):
   if (step < STEP_NAMES.ASSEMBLE or step > STEP_NAMES.ANNOTATE):
      print "Error: unsupported step #%d\n"%(step)
      raise(JobSignalledBreak)

   if not checkIfExists(step, programName):
      print "Error: attempted to execute program %s which is not supported\n"%(programName)
      raise(JobSignalledBreak)

   global _settings
   _settings = settings

   dispatch = GenericProgram(programName, step)
   dispatch.read()
   return dispatch.execute()
