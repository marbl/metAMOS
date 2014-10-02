#!python

import os, sys, string, time, BaseHTTPServer, getopt, re, subprocess, webbrowser
from datetime import date
from datetime import time
from datetime import datetime
from operator import itemgetter
from utils import readConfigInfo
from utils import str2bool

_workflows = dict()
_verbose = False

class Workflow:
   name = ""
   programList = set()
   super = []
   stepList = []
   readlibs = []
   commandList = ""
   path = set()

   def __init__(self, name = "", path = "%s%sUtilities%sworkflows"%(sys.path[0], os.sep, os.sep)):
      self.name = name.lower()
      self.md5 = None
      self.super = []
      self.programList = set()
      self.readlibs = []
      self.asmcontigs = []
      self.stepList = ["Preprocess","Assemble","MapReads","FindORFS","Classify","Abundance","Scaffold","Propagate","Bin","Postprocess"]
      self.commandList = ""
      self.isValid = False
      self.isModifiable = False
      if isinstance(path, set):
         self.path.union(path)
      else:
         self.path.add(path)
      self.path.add("%s%sUtilities%sworkflows"%(sys.path[0], os.sep, os.sep))

   def getDerivedName(self):
      nameList = self.super
      nameList.append(self.name)
      return nameList

   def getMD5(self):
      md5 = hashlib.md5()
      md5.update(self.isValid)
      return md5.hexdigest()

   def canModify(self):
      return len(self.commandList) == 0 or self.isModifiable == True

   def getRequiresPrograms(self):
      workflow_file = None
      for p in self.path:
         if os.path.exists("%s%s%s.ini"%(p, os.sep, self.name.lower())):
            workflow_file = open("%s%s%s.ini"%(p, os.sep, self.name.lower()), 'r')
      if workflow_file == None:
         print "Error: cannot open workflow %s"%(self.name)
         raise(JobSignalledBreak)

      lastPos = 0
      for line in workflow_file.xreadlines():
         if "programs:" in line:
            splitLine = line.replace("programs:", "").strip().split()
            for prog in splitLine:
               self.programList.add(prog.lower())
            break
      workflow_file.close()
      return len(self.programList) > 0

   def read(self, comment = "#"):
      # populate from file
      self.isValid = True

      workflow_file = None
      for p in self.path:
         if os.path.exists("%s%s%s.ini"%(p, os.sep, self.name.lower())):
            workflow_file = open("%s%s%s.ini"%(p, os.sep, self.name.lower()), 'r')
      if workflow_file == None:
         print "Error: cannot open workflow %s"%(self.name)
         raise(JobSignalledBreak)

      lastPos = 0
      while True:
         lastPos = workflow_file.tell()
         line = workflow_file.readline()
         (line, sep, commentLine) = line.partition(comment)
         if len(commentLine) != 0:
            continue

         if "programs:" in line:
            splitLine = line.replace("programs:", "").strip().split()
            for prog in splitLine:
               self.programList.add(prog.lower())

         if "steps:" in line:
            self.stepList = []
            splitLine = line.replace("steps:", "").strip().split()
            for step in splitLine:
               if step not in stepList:
                   self.stepList.append(step.lower())

         elif "md5:" in line:
            self.md5 = line.replace("md5:", "").strip()
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
               self.commandList = (masterFlow.commandList + " " + self.commandList).strip()
               self.asmcontigs.extend(masterFlow.asmcontigs)
               self.isModifiable = masterFlow.isModifiable or self.isModifiable
         elif "modify:" in line:
            splitLine = line.replace("modify:", "").strip()
            self.isModifiable = self.isModifiable or str2bool(splitLine)
         elif "command:" in line:
            splitLine = line.replace("command:", "").strip()
            self.commandList = (self.commandList + " " + splitLine).strip()
         else:
            break
      workflow_file.seek(lastPos)

      # now read the library info
      (readasm, readl, ignore) = readConfigInfo(workflow_file)
      self.readlibs.extend(readl)

      self.asmcontigs.extend(readasm)
      workflow_file.close()

      if self.md5 != None:
        if self.md5 != getMD5():
           self.isValid = False
 
      if _verbose:
         print "Read workflow %s"%(self.name)
         print "Inherits from %s"%(",".join(self.super))
         print "Supported programs %s"%(",".join(self.programList))
         print "Asm contigs %s"%(",".join(self.asmcontigs))
         for lib in self.readlibs:
            print "Lib %s %s %s %s %s %s %s %s %s"%(lib.id, lib.format, lib.mean, lib.stdev, lib.mmin, lib.mmax, lib.mated, lib.interleaved, lib.innie) 

def getSupportedWorkflows(path = "%s%sUtilities%sworkflows"%(sys.path[0], os.sep, os.sep), includeRuntime = False):
   enabled = [ ]

   if os.path.exists("%s%sworkflows.ini"%(path, os.sep)):
      if os.path.exists("%s%senabled.ini"%(path, os.sep)):
         workflow_file = open("%s%senabled.ini"%(path, os.sep), 'r')
         for line in workflow_file.xreadlines():
            workflowName = line.strip()
            workflow = Workflow(workflowName, path)
            workflow.read()
            enabled.append(workflow)
         workflow_file.close()
      else:
         workflow = Workflow("core", path)
         workflow.read()
         enabled.append(workflow)

   # finally grab any workflows that don't have any program requirements
   if includeRuntime:
      for file in os.listdir("%s%s"%(path, os.sep)):
         if not file.endswith(".ini"):
            continue

         prefix = os.path.splitext(os.path.basename(file))[0]
         if prefix == "workflow" or prefix == "enabled":
            continue
         else:
            workflow = Workflow(prefix, path)
            if workflow.getRequiresPrograms() == False:
               workflow.read()
               enabled.append(workflow)

   return enabled
   
def getSupportedWorkflowNames(path = "%s%sUtilities%sworkflows"%(sys.path[0], os.sep, os.sep), includeRuntime = False):
   wfs = getSupportedWorkflows(path, includeRuntime)

   names = set()
   for wf in wfs:
      names.update(wf.getDerivedName())

   return names

def updateSupportedWorkflows(workflowList, path = "%s%sUtilities%sworkflows"%(sys.path[0], os.sep, os.sep)):
   workflow_file = open("%s%senabled.ini"%(path, os.sep), 'w')
   workflow_file.write("\n".join(workflowList))
   workflow_file.close()

def getAllWorkflows(path = "%s%sUtilities%sworkflows"%(sys.path[0], os.sep, os.sep)):
   global _workflows

   if (len(_workflows) == 0):
     if os.path.exists("%s%sworkflows.ini"%(path, os.sep)):
        workflow_file = open("%s%sworkflows.ini"%(path, os.sep), 'r') 
        for line in workflow_file.xreadlines():
           workflowName = line.strip().lower()
           workflow = Workflow(workflowName, path)
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

if __name__ == "__main__":
   if os.path.exists("%s/Utilities%sworkflows%s%s.ini"%(sys.path[0], os.sep, os.sep, sys.argv[1])):
      workflow = Workflow(sys.argv[1], "%s%sUtilities%sworkflows"%(sys.path[0], os.sep, os.sep))
      workflow.read()

