#!python

import os, sys, string, time, BaseHTTPServer, getopt, re, subprocess, webbrowser
from datetime import date
from datetime import time
from datetime import datetime
from operator import itemgetter
from utils import *

_workflows = dict()
_verbose = False

class Workflow:
   name = ""
   programList = set()
   super = []
   readlibs = []
   commandList = []
   path = None

   def __init__(self, name = "", path = sys.path[0]):
      self.name = name.lower()
      self.super = []
      self.programList = set()
      self.readlibs = []
      self.asmcontigs = []
      self.commandList = []
      self.isValid = False 
      self.path = path

   def getName(self):
      nameList = self.super
      nameList.append(self.name)
      return nameList

   def read(self):
      # populate from file
      self.isValid = True

      workflow_file = open("%s%s%s.ini"%(self.path, os.sep, self.name.lower()), 'r')

      lastPos = 0
      while True:
         lastPos = workflow_file.tell()
         line = workflow_file.readline()

         if "programs:" in line:
            splitLine = line.replace("programs:", "").strip().split()
            for prog in splitLine:
               self.programList.add(prog.lower())
         elif "inherit:" in line:
            splitLine = line.replace("inherit:", "").strip().split()
            for master in splitLine:
               self.super.append(master)
               masterFlow = None
               if master in _workflows:
                  masterFlow = _workflows[master]
               else:
                  masterFlow = Workflow(line.replace("inherit:", "").strip(), self.path)
                  masterFlow.read()
               self.programList.update(masterFlow.programList)
               self.readlibs.extend(masterFlow.readlibs)
               self.commandList.extend(masterFlow.asmcontigs)
               self.asmcontigs.extend(masterFlow.asmcontigs)
         else:
            break
      workflow_file.seek(lastPos)

      # now read the library info
      (readasm, readl) = readConfigInfo(workflow_file)
      self.readlibs.extend(readl)
      self.asmcontigs.extend(readasm)
      workflow_file.close()

      if _verbose:
         print "Read workflow %s"%(self.name)
         print "Inherits from %s"%(",".join(self.super))
         print "Supported programs %s"%(",".join(self.programList))
         print "Asm contigs %s"%(",".join(self.asmcontigs))
         for lib in self.readlibs:
            print "Lib %s %s %s %s %s %s %s %s %s"%(lib.id, lib.format, lib.mean, lib.stdev, lib.mmin, lib.mmax, lib.mated, lib.interleaved, lib.innie) 


def getSupportedWorkflows(path = sys.path[0]):
   enabled = []

   if os.path.exists("%s%sUtilities%sworkflows%sworkflows.ini"%(path, os.sep, os.sep, os.sep)):
      workflow_file = open("%s%sUtilities%sworkflows%senabled.ini"%(path, os.sep, os.sep, os.sep), 'r')
      for line in workflow_file.xreadlines():
         enabled.append(line.strip())

      workflow_file.close()
   else:
      print "Error: unable to find any workflows"
      raise JobSignalledBreak

   return enabled
   
def updateSupportedWorkflows(workflowList, path = sys.path[0]):
   workflow_file = open("%s%sUtilities%sworkflows%senabled.ini"%(path, os.sep, os.sep, os.sep), 'w')
   workflow_file.write("\n".join(workflowList))
   workflow_file.close()

def getAllWorkflows(path = sys.argv[0]):
   global _workflows

   if (len(_workflows) == 0):
     if os.path.exists("%s%sUtilities%sworkflows%sworkflows.ini"%(path, os.sep, os.sep, os.sep)):
        workflow_file = open("%s%sUtilities%sworkflows%sworkflows.ini"%(path, os.sep, os.sep, os.sep), 'r') 
        for line in workflow_file.xreadlines():
           workflowName = line.strip()
           workflow = Workflow(workflowName, "%s%sUtilities%sworkflows"%(path, os.sep, os.sep))
           workflow.read()
           _workflows[workflowName] = workflow 
        workflow_file.close()
     else:
        print "Error: unable to find workflows, defaulting to core only"
        _workflows["core"] = Workflow("core")
        plist = ["setuptools", "pysam", "psutil", "cython", "numpy", "matplotlib", "AMOS", "kraken", "LAP", "KronaTools", "fastqc", "uniprot"]
        _workflows["core"].programList = set(plist)
        _workflows["core"].isValid = True

   return _workflows
