import os, sys, string, subprocess, distutils.util, site, glob, multiprocessing

def addEnvironmentVar(varName, newValue, sep = " "):
   oldVal = ""
   if varName in os.environ:
      oldVal = os.environ[varName]
      os.environ[varName] = newValue + sep + oldVal
   else:
      os.environ[varName] = newValue
   return oldVal

def updateMakeFileForDarwin(fileName, addedCFlags, addedLDFlags, addFlagsToCompile=False):
   if OSTYPE == "Darwin":
      os.system("cp %s %s.orig"%(fileName, fileName))
      numCF=utils.getCommandOutput("grep -c \"CFLAGS*=\" %s.orig"%(fileName), False).strip()
      numCX=utils.getCommandOutput("grep -c \"CXXFLAGS*=\" %s.orig"%(fileName), False).strip()
      numLD=utils.getCommandOutput("grep -c \"LDFLAGS*=\" %s.orig"%(fileName), False).strip()
      numD=utils.getCommandOutput("grep -c \^DFLAGS*=\" %s.orig"%(fileName), False).strip()

      addCF = False
      addLD = False
      if ((numCF == "" or int(numCF) == 0) and (numCX == "" or int(numCX) == 0)):
         addCF = True
      if ((numCF == "" or int(numCF) == 0) and (numD == "" or int(numD) == 0)):
         addLD = True

      os.system("cat %s.orig |awk '{if (match($0, \"^CFLAGS.*=\")) { print $0\" %s\"; } else if (match($0, \"^CXXFLAGS.*=\")) { print $0\" %s\"; } else if (match($0, \"^LDFLAGS.*=\")) { print $0\" %s\" } else if (match($0, \"^DFLAGS =\")) { print $0\" %s\"; } else { print $0; } }' >%s"%(fileName, addedCFlags, addedCFlags, addedLDFlags, addedLDFlags, fileName))
      if addCF:
         os.system("cp %s %s.orig"%(fileName, fileName))
         os.system("cat %s.orig |awk '{if (NR == 1) { print \"CFLAGS=%s\\nCXXFLAGS=%s\\n\"$0; } else { print $0; } }' > %s"%(fileName, addedCFlags, addedCFlags, fileName))
      if addLD:
         os.system("cp %s %s.orig"%(fileName, fileName))
         os.system("cat %s.orig |awk '{if (NR == 1) { print \"LDFLAGS=%s\\n\"$0; } else { print $0; } }' > %s"%(fileName, addedLDFlags, fileName))

      if addFlagsToCompile:
         os.system("cp %s %s.orig"%(fileName, fileName))
         os.system("cat %s.orig |awk '{if (match($1, \"g++\")) { sub(/g\\+\\+/, \"g++ \\$(CXXFLAGS) \\$(LDFLAGS)\", $0) } print $0; }' > %s"%(fileName, fileName))

def copyPerlLib(pathToCopy, dest):
   if pathToCopy != "":
      pathsToCopy = pathToCopy.strip().split("\n")
      for path in pathsToCopy:
         pathToCopy = os.path.dirname(path)
         os.system("mkdir -p %s"%(dest))
         # copy one at a time in case of conflicts
         for file in os.listdir("%s%s"%(pathToCopy, os.sep)):
            toCopy = file
            file = "%s%s%s"%(pathToCopy, os.sep, toCopy)
            if os.path.exists("%s/%s"%(dest, toCopy)):
               os.system("mv %s/* %s/%s/"%(file, dest, toCopy))
            else:
               os.system("mv %s %s/"%(file, dest))


user_home = os.environ["HOME"]
print "<<Welcome to metAMOS install>>"

#check for python version
if (sys.version_info[0] < 2) or (sys.version_info[0] == 2 and sys.version_info[1] < 6):
  print "Python version is %s. metAMOS requires at least 2.6"%(sys.version)
  sys.exit(1)

#add access to utils.py, for utils dir
METAMOS_ROOT  = os.getcwd().strip()
INITIAL_SRC   = "%s%ssrc"%(METAMOS_ROOT, os.sep)
sys.path.append(INITIAL_SRC)
import utils
import workflow
sys.path.append(utils.INITIAL_UTILS)

shellv = os.environ["SHELL"]

#add site dir
site.addsitedir(utils.INITIAL_UTILS+os.sep+"python"+os.sep+"lib"+os.sep+"python")
site.addsitedir(utils.INITIAL_UTILS+os.sep+"python"+os.sep+"lib64"+os.sep+"python")

if "PYTHONPATH" not in os.environ:
   os.environ["PYTHONPATH"] = ""
os.environ["PYTHONPATH"]+=utils.INITIAL_UTILS+os.sep+"python"+os.pathsep
os.environ["PYTHONPATH"] += utils.INITIAL_UTILS+os.sep+"python"+os.sep+"lib"+os.pathsep
os.environ["PYTHONPATH"] += utils.INITIAL_UTILS+os.sep+"python"+os.sep+"lib"+os.sep+"python"+os.pathsep
os.environ["PYTHONPATH"] += utils.INITIAL_UTILS+os.sep+"python"+os.sep+"lib64"+os.pathsep
os.environ["PYTHONPATH"] += utils.INITIAL_UTILS+os.sep+"python"+os.sep+"lib64"+os.sep+"python"+os.pathsep
sys.path.append(utils.INITIAL_UTILS+os.sep+"python")
sys.path.append(utils.INITIAL_UTILS+os.sep+"python" + os.sep+"lib"+ os.sep+"python")
sys.path.append(utils.INITIAL_UTILS+os.sep+"python" + os.sep+"lib64"+ os.sep+"python")

if 'bash' in shellv or utils.cmdExists('export'):
   os.system("export PYTHONPATH=%s:$PYTHONPATH"%(utils.INITIAL_UTILS+os.sep+"python"))
   os.system("export PYTHONPATH=%s:$PYTHONPATH"%(utils.INITIAL_UTILS+os.sep+"python"+os.sep+"lib"+os.sep+"python"))
   os.system("export PYTHONPATH=%s:$PYTHONPATH"%(utils.INITIAL_UTILS+os.sep+"python"+os.sep+"lib64"+os.sep+"python"))
elif utils.cmdExists('setenv'):
   os.system("setenv PYTHONPATH %s:$PYTHONPATH"%(utils.INITIAL_UTILS+os.sep+"python"))
   os.system("setenv PYTHONPATH %s:$PYTHONPATH"%(utils.INITIAL_UTILS+os.sep+"python"+os.sep+"lib"+os.sep+"python"))
   os.system("setenv PYTHONPATH %s:$PYTHONPATH"%(utils.INITIAL_UTILS+os.sep+"python"+os.sep+"lib64"+os.sep+"python"))
else:
   print "Cannot set PYTHONPATH variable, unknown shell %s\n"%(shellv)

if not os.path.exists("%s"%utils.INITIAL_UTILS+os.sep+"python"+os.sep+"lib"):
    os.system("mkdir %s"%utils.INITIAL_UTILS+os.sep+"python"+os.sep+"lib")
if not os.path.exists("%s"%utils.INITIAL_UTILS+os.sep+"python"+os.sep+"lib64"):
    os.system("mkdir %s"%utils.INITIAL_UTILS+os.sep+"python"+os.sep+"lib64")
if not os.path.exists("%s"%utils.INITIAL_UTILS+os.sep+"python"+os.sep+"lib"+os.sep+"python"):
    os.system("mkdir %s"%utils.INITIAL_UTILS+os.sep+"python"+os.sep+"lib"+os.sep+"python")
if not os.path.exists("%s"%utils.INITIAL_UTILS+os.sep+"python"+os.sep+"lib64"+os.sep+"python"):
    os.system("mkdir %s"%utils.INITIAL_UTILS+os.sep+"python"+os.sep+"lib64"+os.sep+"python")

ALLOW_FAST=True
HAVE_GCC42=False
HAVE_RT=False
HAVE_QUIET_HEAD=False
GCC_VERSION=0.0
try:
   GCC_VERSION=float(utils.getCommandOutput("gcc --version|grep gcc|awk '{print $NF}' |awk -F \".\" '{print $1\".\"$2}'", False))
except:
   try:
      GCC_VERSION=float(utils.getCommandOutput("gcc --version|grep gcc|awk '{print $3}' |awk -F \".\" '{print $1\".\"$2}'", False))
   except:
      print "Warning: cannot determine GCC version"

OSTYPE="Linux"
OSVERSION="1"
MACHINETYPE="x86_64"
kronaTools = "KronaTools-2.2"

#identify machine type
p = subprocess.Popen("echo `uname`", shell=True, stdin=None, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
(checkStdout, checkStderr) = p.communicate()
if checkStderr != "":
   print "Warning: Cannot determine OS, defaulting to %s"%(OSTYPE)
else:
    OSTYPE = checkStdout.strip()

p = subprocess.Popen("echo `uname -r`", shell=True, stdin=None, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
(checkStdout, checkStderr) = p.communicate()
if checkStderr != "":
   print "Warning: Cannot determine OS version, defaulting to %s"%(OSVERSION)
else:
   OSVERSION = checkStdout.strip()

p = subprocess.Popen("echo `uname -m`", shell=True, stdin=None, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
(checkStdout, checkStderr) = p.communicate()
if checkStderr != "":
   print "Warning: Cannot determine system type, defaulting to %s"%(MACHINETYPE)
else:
   MACHINETYPE = checkStdout.strip()

addedCFlags=""
addedLDFlags=""
oldCFlags = ""
oldCPPFlags = ""
oldCXXFlags = ""
oldLDFlags = ""

if OSTYPE == "Darwin":
   p = subprocess.Popen("echo `gcc --version`", shell=True, stdin=None, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
   (checkStdout, checkStderr) = p.communicate()
   if "Apple" not in checkStdout:
      ALLOW_FAST=False
   gcc42 = utils.getCommandOutput("which g++-4.2", False)
   if gcc42 == "":
      HAVE_GCC42=False
   else:
      HAVE_GCC42=True

   # global vars for building
   libPath=""
   clib=utils.getCommandOutput("g++ -print-file-name=libgcc.a", False)
   if clib != "":
      libPath="%s %s"%(libPath, clib)
   cpplib=utils.getCommandOutput("g++ -print-file-name=libstdc++.a", False)
   if cpplib != "":
      libPath="%s %s"%(libPath, cpplib)
   omplib=utils.getCommandOutput("g++ -print-file-name=libgomp.a", False)
   if omplib != "":
      libPath="%s %s"%(libPath, omplib)

   commonFlags="-mmacosx-version-min=10.6 -static-libgcc -static-libstdc++ "
   oldCFlags = addEnvironmentVar("CFLAGS", " %s "%(commonFlags))
   oldCPPFlags = addEnvironmentVar("CPPFLAGS", " %s "%(commonFlags))
   oldCXXFlags = addEnvironmentVar("CXXFLAGS", " %s "%(commonFlags))
   oldLDFlags = addEnvironmentVar("LDFLAGS", " %s "%(libPath))
   addedCFlags="%s %s"%(commonFlags, libPath)
   addedLDFlags="-static-libgcc -static-libstdc++ %s"%(libPath)

libPaths = [ "/usr/lib", "/usr/lib64", "/usr/local/lib/", "/usr/local/lib64/", "/opt/local/lib/", "/opt/local/lib64/"] 
for libPath in libPaths:
   if os.path.exists(libPath + os.sep + "librt.a") or os.path.exists(libPath + os.sep + "librt.so"):
      HAVE_RT=True
      break

p = utils.getCommandOutput("head --help |grep \"\\-q\" |wc -l", False)
if int(p) >= 1:
   HAVE_QUIET_HEAD=True

# get list of supported workflows
enabledWorkflows = set()
packagesToInstall = set()
knownPackages = set()
workflows = workflow.getAllWorkflows("%s/Utilities/workflows"%(METAMOS_ROOT))
for flow in workflows:
   knownPackages.update(workflows[flow].programList)

manual = False
fail = False
nodbs = False

availableWf = workflow.getSupportedWorkflows("%s/Utilities/workflows"%(METAMOS_ROOT), False)
for wf in availableWf:
   enabledWorkflows.update(wf.getDerivedName())
   packagesToInstall.update(wf.programList)

if (len(sys.argv) > 1):
  # should support tool list as well added
  for i in range(1, len(sys.argv)):
      arg = sys.argv[i]
      if arg.lower() in workflows.keys():
         packagesToInstall.update(workflows[arg.lower()].programList)
         enabledWorkflows.update(workflows[arg.lower()].getDerivedName())
      elif arg.lower() == "full":
         for flow in workflows:
             packagesToInstall.update(workflows[flow].programList)
             enabledWorkflows.update(workflows[flow].getDerivedName())
         print "Installing all available workflows"
      elif arg.lower() == "manual":
        manual = True
        for flow in workflows:
           enabledWorkflows.update(workflows[flow].getDerivedName())
      elif arg.lower() == "nodbs":
        nodbs = True

      elif arg.lower() in knownPackages:
         packagesToInstall.add(arg.lower())
         for flow in workflows:
            if arg.lower() in workflows[flow].programList:
               enabledWorkflows.update(workflows[flow].getDerivedName())
               break
      else:
         if arg != "help":
            print "Unknown program or workflow %s specified."%(arg)
         fail = True

if fail or help in sys.argv:
   print "Available workflows: %s"%(" ".join(workflows.keys()))
   print "Available packages: %s"%("\n\t".join(knownPackages))
   exit(1)
    
if manual:
    packagesToInstall = set()

for workflowName in enabledWorkflows:
    print "Selected to install workflowName: %s."%(workflowName.upper())

print "Will automatically install:"
for p in packagesToInstall:
    print "\t%s"%(p.title())

if not os.path.exists("./Utilities/config/usage.ok") and not os.path.exists("./Utilities/config/usage.no"):
    print "MetAMOS would like to record anonymous usage statistics, is this ok ? "
    dl = raw_input("Enter Y/N: ")
    if dl == 'y' or dl == 'Y':
        os.system("echo ok > ./Utilities/config/usage.ok")
    else:
        os.system("echo no > ./Utilities/config/usage.no")

# first the needed python packages
# make sure we have setuptools available
if 1:
   fail = 0
   try:
       import setuptools
   except ImportError:
       fail = 1
   if "setuptools" in packagesToInstall:
       dl = 'y'
   elif fail:
       print "setuptools not found, required for install, download now?"
       dl = raw_input("Enter Y/N: ")
   if fail and (dl == 'y' or dl == "Y"):
       os.system("curl -L https://bitbucket.org/pypa/setuptools/raw/0.7.4/ez_setup.py -o ez_setup.py")
       os.system("python ez_setup.py --user")
       
if 1:
   fail = 0
   try:
       import psutil
   except ImportError:
       fail = 1
   if "psutil" in packagesToInstall:
       dl = 'y'
   elif fail:
       print "psutil not found, required for memory usage estimation, download now?"
       dl = raw_input("Enter Y/N: ")
   if fail and (dl == 'y' or dl == "Y"):
       os.system("curl -L http://psutil.googlecode.com/files/psutil-0.6.1.tar.gz -o ./psutil.tar.gz")
       os.system("tar -C ./Utilities/python -xvf psutil.tar.gz")
       os.system("mv ./Utilities/python/psutil-0.6.1 ./Utilities/python/psutil")
       os.chdir("./Utilities/python/psutil")
       os.system("python setup.py install --home=%spython"%(utils.INITIAL_UTILS+os.sep))
       os.chdir("%s"%(METAMOS_ROOT))
       os.system("rm -rf psutil.tar.gz")
if 1:
   fail = 0
   try:
       import cython
   except ImportError:
       fail = 1
   if "cython" in packagesToInstall:
       dl = 'y'
   elif fail:
       print "cython modules not found, necessary for c-compiling python code, download now?"
       dl = raw_input("Enter Y/N: ")
   if fail and (dl == 'y' or dl == "Y"):
       os.system("curl -L https://github.com/cython/cython/archive/master.zip -o ./cython.zip")
       os.system("unzip ./cython.zip")
       os.system("mv ./cython-master ./Utilities/python/cython")
       os.chdir("./Utilities/python/cython")
       os.system("python setup.py install --home=%spython"%(utils.INITIAL_UTILS+os.sep))
       os.chdir(METAMOS_ROOT)
       os.system("rm -rf cython.zip")

if 1:
   fail = 0
   try:
       import pysam
   except ImportError:
       fail = 1

   if "pysam" in packagesToInstall:
       dl = 'y'
   elif fail:
       print "pysam python modules not found, necessary for bowtie2 alignments, download now?"
       dl = raw_input("Enter Y/N: ")

   if fail and (dl == 'y' or dl == "Y"):
       os.system("curl -L http://pysam.googlecode.com/files/pysam-0.6.tar.gz -o ./pysam.tar.gz")
       os.system("tar -C ./Utilities/python -xvf pysam.tar.gz")
       os.system("mv ./Utilities/python/pysam-0.6 ./Utilities/python/pysam")
       doInstall = True
       #for root install
       #os.system("sudo python ./Utilities/python/pysam/setup.py install")
       os.chdir("./Utilities/python/pysam")
       if OSTYPE == "Darwin":
          if utils.getFromPath("llvm-gcc-4.2", "LLVM GCC"):
             os.system("export CC=llvm-gcc-4.2")
             os.system("export CXX=llvm-g++-4.2")
          else:
             print "Warning: Cannot install pysam on your system. Please install LLVM compiler first."
             doInstall=False
       if doInstall:
          os.system("python setup.py build_ext --inplace")
          os.system("python setup.py build")
          os.system("python setup.py install --home=%spython"%(utils.INITIAL_UTILS+os.sep))
       os.chdir(METAMOS_ROOT)
       os.system("rm -rf pysam.tar.gz")

#WARNING: matplotlib causes install issues for multiple users
   fail = 0
   try:
       import numpy
   except ImportError:
       fail = 1

   if "numpy" in packagesToInstall:
       dl = 'y'
   elif fail:
       print "numpy python modules not found, necessary for html report, download now?"
       dl = raw_input("Enter Y/N: ")
   if fail and (dl == 'y' or dl == "Y"):
       os.system("curl -L http://downloads.sourceforge.net/project/numpy/NumPy/1.7.1/numpy-1.7.1.tar.gz -o ./numpy.tar.gz")
       os.system("tar -C ./Utilities/python -xvf numpy.tar.gz")
       os.system("mv ./Utilities/python/numpy-1.7.1 ./Utilities/python/numpy")
       os.chdir("./Utilities/python/numpy")
       os.system("python setup.py install --home=%spython"%(utils.INITIAL_UTILS+os.sep))
       os.chdir(METAMOS_ROOT)
       os.system("rm -rf numpy.tar.gz")

if 1:
   fail = 0
   try:
       import matplotlib
       if (matplotlib.__version__ < "1.1.0"):
          fail = 1
   except ImportError:
       fail = 1

   if "matplotlib" in packagesToInstall:
       dl = 'y'
   elif fail:
       print "Matplot lib version %s is incompatible with metAMOS. Need version 1.1.0+, download now?"%(matplotlib.__version__) 
       dl = raw_input("Enter Y/N: ")
   if fail and (dl == 'y' or dl == "Y"):
       os.system("curl -L http://downloads.sourceforge.net/project/matplotlib/matplotlib/matplotlib-1.1.0/matplotlib-1.1.0.tar.gz -o ./matplotlib.tar.gz")
       os.system("tar -C ./Utilities/python -xvf matplotlib.tar.gz")
       os.system("mv ./Utilities/python/matplotlib-1.1.0 ./Utilities/python/matplotlib")
       os.chdir("./Utilities/python/matplotlib")
       os.system("python setup.py install --home=%spython"%(utils.INITIAL_UTILS+os.sep))
       os.chdir(METAMOS_ROOT)
       os.system("rm -rf matplotlib.tar.gz")

# now software
if not os.path.exists("./AMOS") or 0:
   if "amos" in packagesToInstall:
       dl = 'y'
   else:
       print "AMOS binaries not found, needed for all steps, download now?"
       dl = raw_input("Enter Y/N: ")
       
   if dl == 'y' or dl == 'Y':
        os.system("curl -L ftp://ftp.cbcb.umd.edu/pub/data/metamos/amos-3.2-BETA-%s-%s.binaries.tar.gz -o ./amos-binaries.tar.gz"%(OSTYPE, MACHINETYPE))
        os.system("tar -xvf amos-binaries.tar.gz")
        os.system("rm -rf amos-binaries.tar.gz")

        # descriptive perl module
        stat = utils.getCommandOutput("perl -MStatistics::Descriptive -e 0 && echo $?", True)
        if stat == "":
           os.system("curl -L ftp://cbcb.umd.edu/pub/data/metamos/Statistics-Descriptive-3.0203.tar.gz -o stat.tar.gz")
           os.system("tar -xvzf stat.tar.gz")
           os.chdir("Statistics-Descriptive-3.0203")
           os.system("perl Makefile.PL PREFIX=`pwd`/build")
           os.system("make install")
           os.chdir("%s"%(METAMOS_ROOT))
           pathToCopy = utils.getCommandOutput("find Statistics-Descriptive-3.0203/build -type d -name \"Statistics\" |grep -v auto", False)
           copyPerlLib(pathToCopy, "AMOS%s%s-%s%slib"%(os.sep, OSTYPE, MACHINETYPE, os.sep))
           os.system("rm -rf stat.tar.gz")
           os.system("rm -rf Statistics-Descriptive-3.0203")

if not os.path.exists("./Utilities/cpp%s%s-%s%skraken"%(os.sep, OSTYPE, MACHINETYPE, os.sep)):
    if "kraken" in packagesToInstall:
       dl = 'y'
    else:
       print "Kraken not found, optional for Annotate step, download now?"
       dl = raw_input("Enter Y/N: ")
    if dl == 'y' or dl == 'Y':
        archive = "kraken.tar.gz"
        os.system("curl -L ftp://ftp.cbcb.umd.edu/pub/data/metamos/kraken-0.10.3-beta.tgz -o %s"%(archive))
        os.system("rm -rf ./Utilities/cpp%s%s-%s%skraken"%(os.sep, OSTYPE, MACHINETYPE, os.sep))
        os.system("tar -xvzf %s"%(archive))
        os.system("mv kraken-0.10.3-beta ./Utilities/cpp/%s%s-%s%skraken"%(os.sep, OSTYPE, MACHINETYPE, os.sep))
        os.chdir("./Utilities/cpp/%s%s-%s%skraken"%(os.sep, OSTYPE, MACHINETYPE, os.sep))
        os.system("./install_kraken.sh `pwd`/bin")
        os.chdir("%s"%(METAMOS_ROOT))
        os.system("rm %s"%archive)

if not os.path.exists("./Utilities/DB/kraken"):
    if "kraken" in packagesToInstall:
       dl = 'y'
    else:
       print "Kraken DB not found, required for Kraken, download now?"
       dl = raw_input("Enter Y/N: ")
    if dl == 'y' or dl == 'Y':
       settings = utils.Settings(1, 1, "", "")
       settings.OSTYPE = OSTYPE
       mem = utils.getAvailableMemory(settings)

       if (mem < 100) and not nodbs:
          print "Insufficient memory to build full Kraken database. Requires at least 100GB of memory, using mini DB"
          archive = "minikraken.tgz"
          os.system("curl -L ftp://ftp.cbcb.umd.edu/pub/data/metamos/%s -o %s"%(archive, archive))
          os.system("tar xvzf %s"%(archive))
          os.system("mv minikraken_* ./Utilities/DB/kraken")
          os.system("rm %s"%(archive))
       elif not nodbs:
          # first we need jellyfish which is used to build DB
          # kraken needs jellyfish, if we don't find it build it and add to path
          jellyfish = utils.getFromPath("jellyfish", "Jellyfish", False)
          if jellyfish == "":
             archive = "jellyfish.tar.gz"
             os.system("curl -L http://www.cbcb.umd.edu/software/jellyfish/jellyfish-1.1.11.tar.gz -o %s"%(archive))
             os.system("tar xvzf %s"%(archive))
             os.system("mv jellyfish-1.1.11 ./Utilities/cpp%s%s-%s%s/jellyfish"%(os.sep, OSTYPE, MACHINETYPE, os.sep))
             os.chdir("./Utilities/cpp%s%s-%s%s/jellyfish"%(os.sep, OSTYPE, MACHINETYPE, os.sep))
             os.system("./configure --prefix=`pwd`")
             os.system("make")
             os.system("make install")
             os.chdir("%s"%(METAMOS_ROOT))

             pathUpdate = "%s/Utilities/cpp%s%s-%s%sjellyfish/bin/"%(METAMOS_ROOT, os.sep, OSTYPE, MACHINETYPE, os.sep)
             if "PATH" in os.environ:
                pathUpdate = "%s%s%s"%(os.environ["PATH"], os.pathsep, pathUpdate)
             os.environ["PATH"]=pathUpdate

          os.chdir("./Utilities/cpp/%s%s-%s%skraken"%(os.sep, OSTYPE, MACHINETYPE, os.sep))
          os.system("./bin/kraken-build --standard --threads %d --db %s/Utilities/DB/kraken"%(multiprocessing.cpu_count() - 1, METAMOS_ROOT)) 
          os.chdir("%s"%(METAMOS_ROOT))

if not os.path.exists("./LAP"):
    if "lap" in packagesToInstall:
       dl = 'y'
    else:
       print "LAP tool not found, needed for multiple assembly pipeline, download now?"
       dl = raw_input("Enter Y/N: ")
    if dl == 'y' or dl == 'Y':
        os.system("curl -L http://www.cbcb.umd.edu/~cmhill/files/lap_release_1.1.zip -o lap_release_1.1.zip")
        os.system("unzip lap_release_1.1.zip")
        os.system("mv ./lap_release_1.1 ./LAP")
        os.system("rm -rf lap_release_1.1.zip")

if not os.path.exists("KronaTools") or 0:
    if "kronatools" in packagesToInstall:
       dl = 'y'
    else:
       print "KronaTools not found, needed for Postprocess, download now?"
       dl = raw_input("Enter Y/N: ")
    if dl == 'y' or dl == 'Y':
        # TODO: KronaTools should be on the FTP site for robustness to URL changes
        os.system("curl -L 'ftp://ftp.cbcb.umd.edu/pub/data/metamos/" + kronaTools + ".tar' -o %s.tar"%(kronaTools))
        os.system("tar -xvf %s.tar"%(kronaTools))
        os.system("rm -rf %s.tar"%(kronaTools))
        os.system("mv %s KronaTools"%(kronaTools))
        os.system("cd KronaTools && ./install.pl --prefix=.")

if not os.path.exists("KronaTools/taxonomy/taxonomy.tab") or 0:
    if "kronatools" in packagesToInstall:
       dl = 'y'
    else:
       print "KronaTools taxonomy data not found, needed for Postprocess, download now (will take around 20 minutes)?"
       dl = raw_input("Enter Y/N: ")
    if (dl == 'y' or dl == 'Y') and not nodbs:
        os.system("cd KronaTools && ./updateTaxonomy.sh")
        os.chdir("%s"%(METAMOS_ROOT))
        os.system("cat KronaTools/taxonomy/taxonomy.tab |awk -F \"\\t\" '{print $1\"\\t\"$NF}' > ./Utilities/DB/tax_key.tab")

if not os.path.exists("./FastQC"):
    if "fastqc" in packagesToInstall:
        dl = 'y'
    else:
       print "FastQC not found, optional for Preprocess, download now?"
       dl = raw_input("Enter Y/N: ")
    if dl == 'y' or dl == 'Y':
        archive = "fastqc_v0.10.0.zip"
        os.system("curl -L http://www.bioinformatics.babraham.ac.uk/projects/fastqc/%s -o %s" % (archive,archive))
        os.system("unzip %s" % archive)
        os.system("rm %s" % archive)
        os.system("chmod u+x FastQC/fastqc")
        
if not os.path.exists("./Utilities/DB/uniprot_sprot.fasta"):
    if "uniprot" in packagesToInstall:
       dl = 'y'
    else:
       print "Uniprot/Swissprot DB not found, optional for Functional Annotation, download now?"
       dl = raw_input("Enter Y/N: ")
    if (dl == 'y' or dl == 'Y') and not nodbs:
        archive = "uniprot.tar.gz"
        os.system("curl -L ftp://ftp.cbcb.umd.edu/pub/data/metamos/%s -o %s" %(archive, archive))
        os.system("tar -C ./Utilities/DB/ -xvf %s" % archive)
        os.system("rm %s"%archive)

# velvet
if not os.path.exists("./Utilities/cpp%s%s-%s%svelvet"%(os.sep, OSTYPE, MACHINETYPE, os.sep)):
    if "velvet" in packagesToInstall:
       dl = 'y'
    else:
       print "Velvet not found, optional for Assemble step, download now?"
       dl = raw_input("Enter Y/N: ")
    if dl == 'y' or dl == 'Y':
        archive = "velvet.tar.gz"
        os.system("curl -L ftp://ftp.cbcb.umd.edu/pub/data/metamos/velvet_1.2.10.tgz -o %s"%(archive))
        os.system("rm -rf ./Utilities/cpp%s%s-%s%svelvet"%(os.sep, OSTYPE, MACHINETYPE, os.sep))
        os.system("tar -xvzf %s"%(archive))
        os.system("mv velvet_1.2.10 ./Utilities/cpp/%s%s-%s%svelvet"%(os.sep, OSTYPE, MACHINETYPE, os.sep))
        os.chdir("./Utilities/cpp/%s%s-%s%svelvet"%(os.sep, OSTYPE, MACHINETYPE, os.sep))
        updateMakeFileForDarwin("Makefile", addedCFlags, addedLDFlags)
        os.system("make clean")
        os.system("make CATEGORIES=16 MAXKMERLENGTH=127 OPENMP=1 BUNDLEDZLIB=1")
        os.chdir("%s"%(METAMOS_ROOT))
        os.system("rm %s"%archive)

# velvet-sc
if not os.path.exists("./Utilities/cpp%s%s-%s%svelvet-sc"%(os.sep, OSTYPE, MACHINETYPE, os.sep)):
    if "velvet-sc" in packagesToInstall:
       dl = 'y'
    else:
       print "Velvet-SC not found, optional for Assemble step, download now?"
       dl = raw_input("Enter Y/N: ")
    if dl == 'y' or dl == 'Y':
        archive = "velvet-sc.tar.gz"
        os.system("curl -L ftp://ftp.cbcb.umd.edu/pub/data/metamos/velvet-sc.tar.gz -o %s"%(archive))
        os.system("rm -rf ./Utilities/cpp%s%s-%s%svelvet-sc"%(os.sep, OSTYPE, MACHINETYPE, os.sep))
        os.system("tar -xvzf %s"%(archive))
        os.system("mv velvet-sc ./Utilities/cpp/%s%s-%s%svelvet-sc"%(os.sep, OSTYPE, MACHINETYPE, os.sep))
        os.chdir("./Utilities/cpp/%s%s-%s%svelvet-sc"%(os.sep, OSTYPE, MACHINETYPE, os.sep))
        updateMakeFileForDarwin("Makefile", addedCFlags, addedLDFlags)
        os.system("make clean")
        os.system("make CATEGORIES=16 MAXKMERLENGTH=127 OPENMP=1")
        os.chdir("%s"%(METAMOS_ROOT))
        os.system("rm %s"%archive)

# metavelvet
if not os.path.exists("./Utilities/cpp%s%s-%s%sMetaVelvet"%(os.sep, OSTYPE, MACHINETYPE, os.sep)):
    if "metavelvet" in packagesToInstall:
       dl = 'y'
    else:
       print "MetaVelvet not found, optional for Assemble step, download now?"
       dl = raw_input("Enter Y/N: ")
    if dl == 'y' or dl == 'Y':
        archive = "MetaVelvet-1.2.02.tgz"
        os.system("curl -L ftp://ftp.cbcb.umd.edu/pub/data/metamos/MetaVelvet-1.2.02.tgz -o %s"%(archive))
        os.system("rm -rf ./Utilities/cpp%s%s-%s%sMetaVelvet"%(os.sep, OSTYPE, MACHINETYPE, os.sep))
        os.system("tar -xvzf %s"%(archive))
        os.system("mv MetaVelvet-1.2.02 ./Utilities/cpp/%s%s-%s%sMetaVelvet"%(os.sep, OSTYPE, MACHINETYPE, os.sep))
        os.chdir("./Utilities/cpp/%s%s-%s%sMetaVelvet"%(os.sep, OSTYPE, MACHINETYPE, os.sep))
        if OSTYPE == "Darwin":
           os.system("cp Utils/Utils.hh Utils/Utils.hh.orig")
           os.system("cat Utils/Utils.hh.orig |awk '{if (match($0, \"#define MAX_STRING_LENGTH\")) { print \"#include <sys/types.h>\\n\"$0; } else { print $0; }}' > Utils/Utils.hh")
           updateMakeFileForDarwin("Makefile", addedCFlags, addedLDFlags)
        os.system("make clean")
        os.system("make CATEGORIES=16 MAXKMERLENGTH=127")
        os.chdir("%s"%(METAMOS_ROOT))
        os.system("rm %s"%archive)
        
if "viritas" in enabledWorkflows or manual:
    if not os.path.exists("./Utilities/cpp%s%s-%s%strnascan"%(os.sep, OSTYPE, MACHINETYPE, os.sep)):
        if "trnascan" in packagesToInstall:
           dl = 'y'
        else:
           print "tRNAscan not found, optional for Annotate step, download now?"
           dl = raw_input("Enter Y/N: ")
        if dl == 'y' or dl == 'Y':
            os.system("curl -L http://lowelab.ucsc.edu/software/tRNAscan-SE.tar.gz -o trnascan.tar")
            os.system("tar xvf trnascan.tar")
            os.system("mv tRNAscan-SE-1.3.1 ./Utilities/cpp/%s%s-%s%strnascan"%(os.sep, OSTYPE, MACHINETYPE, os.sep))
            os.chdir("./Utilities/cpp/%s%s-%s%strnascan"%(os.sep, OSTYPE, MACHINETYPE, os.sep))
            updateMakeFileForDarwin("Makefile", addedCFlags, addedLDFlags)
            os.system("make")
            os.chdir("%s"%(METAMOS_ROOT))
            os.system("rm -rf trnascan.tar")

# now workflow specific tools
if "optional" in enabledWorkflows or manual:
    if not os.path.exists("./Utilities/cpp/%s-%s/metaphylerClassify"%(OSTYPE, MACHINETYPE)) or not os.path.exists("./Utilities/perl/metaphyler/markers/markers.protein") or not os.path.exists("./Utilities/perl/metaphyler/markers/markers.dna"):
        if "metaphyler" in packagesToInstall:
           dl = 'y'
        else:
           print "Metaphyler (latest version) not found, optional for Annotate, download now?"
           dl = raw_input("Enter Y/N: ")
        if dl == 'y' or dl == 'Y':
            os.system("curl -L http://metaphyler.cbcb.umd.edu/MetaPhylerV1.25.tar.gz -o metaphyler.tar.gz")
            os.system("tar -C ./Utilities/perl/ -xvf metaphyler.tar.gz")
            os.system("mv ./Utilities/perl/MetaPhylerV1.25 ./Utilities/perl/metaphyler")
            os.system("mv ./Utilities/perl/metaphyler/installMetaphyler.pl ./Utilities/perl/metaphyler/installMetaphylerFORMATDB.pl");
            os.system("cat ./Utilities/perl/metaphyler/installMetaphylerFORMATDB.pl  |sed 's/formatdb/\.\/Utilities\/cpp\/%s-%s\/formatdb/g' > ./Utilities/perl/metaphyler/installMetaphyler.pl"%(OSTYPE, MACHINETYPE));
            os.system("perl ./Utilities/perl/metaphyler/installMetaphyler.pl")
            os.system("cp ./Utilities/perl/metaphyler/metaphylerClassify ./Utilities/cpp/%s-%s/metaphylerClassify"%(OSTYPE, MACHINETYPE))

    if not os.path.exists("./Utilities/models") or not os.path.exists("./Utilities/DB/blast_data"):
        if "fcp" in packagesToInstall:
           dl = 'y'
        else:
           print "Genome models not found, optional for FCP/NB, download now?"
           dl = raw_input("Enter Y/N: ")
        if (dl == 'y' or dl == 'Y') and not nodbs:
            archive = "fcp_models.tar.gz"
            os.system("curl -L ftp://ftp.cbcb.umd.edu/pub/data/metamos/%s -o %s" %(archive, archive))
            os.system("rm -rf ./Utilities/DB/blast_data")
            os.system("rm -rf ./Utilities/models")
            os.system("tar -C ./Utilities/ -xvf %s" % archive)
            os.system("rm %s"%archive)

    if not os.path.exists("./phylosift") or not os.path.exists("./phylosift/lib/version.pm") or not os.path.exists("./phylosift/lib/Params"):
       if "phylosift" in packagesToInstall:
          dl = 'y'
       else:
          print "PhyloSift binaries not found, optional for Annotate step, download now?"
          dl = raw_input("Enter Y/N: ")
       if dl == 'y' or dl == 'Y':
          if not os.path.exists("./phylosift"):
             #phylosift OSX binaries included inside Linux X86_64 tarball..
             os.system("curl -L http://edhar.genomecenter.ucdavis.edu/~koadman/phylosift/devel/phylosift_20130829.tar.bz2 -o ./phylosift.tar.bz2")
             os.system("tar -xvjf phylosift.tar.bz2")
             os.system("rm -rf phylosift.tar.bz2")
             os.system("mv phylosift_20130829 phylosift")

          if not os.path.exists("./phylosift/legacy/version.pm"):
             #phylosift needs version but doesn't include it
             os.system("curl -L http://www.cpan.org/authors/id/J/JP/JPEACOCK/version-0.9903.tar.gz -o version.tar.gz")
             os.system("tar xvzf version.tar.gz")
             os.chdir("./version-0.9903/")
             os.system("perl Makefile.PL")
             os.system("make")
             os.system("cp -r blib/lib/* ../phylosift/lib")
             os.chdir(METAMOS_ROOT)
             os.system("rm -rf version.tar.gz")
             os.system("rm -rf version-0.9903")
          if not os.path.exists("./phylosift/lib/Params"):
             os.system("curl -L ftp://ftp.cbcb.umd.edu/pub/data/metamos/params-validate.tar.gz -o ./params-validate.tar.gz")
             os.system("tar xvzf params-validate.tar.gz")
             os.system("rm -rf params-validate.tar.gz")

          # download markers dbs
          if not os.path.exists("./phylosift/share"):
             markerUrl = utils.getCommandOutput("cat phylosift/phylosiftrc |grep marker_base |awk '{print $NF}' |sed s/\;//g", False)
             ncbiUrl = utils.getCommandOutput("cat phylosift/phylosiftrc |grep ncbi_url |awk '{print $NF}' |sed s/\;//g", False)
             os.system("mkdir -p ./phylosift/share/phylosift")
             os.chdir("./phylosift/share/phylosift")
             os.system("curl -L %s/markers.tgz -o marker.tgz"%(markerUrl))
             os.system("tar xvzf marker.tgz")
             os.system("rm marker.tgz")
             os.system("curl -L %s -o ncbi.tgz"%(ncbiUrl))
             os.system("tar xvzf ncbi.tgz") 
             os.system("rm ncbi.tgz")
             os.chdir(METAMOS_ROOT)

    # check the number of files the DB currently is and see if we have the expected number
    dbResult = ""
    if not nodbs:
        dbResult = utils.getCommandOutput("perl ./Utilities/perl/update_blastdb.pl refseq_protein --numpartitions", False)
    if not nodbs and dbResult == "":
       print "Error: could not connect to NCBI, will not be installing refseq protein DB"
    elif not nodbs:
       (dbName, numPartitions) = dbResult.split("\t", 1) 
       print "Checking whether %s is complete. Expecting %d partitions.\n"%(dbName, int(numPartitions))
       numPartitions = int(numPartitions) - 1
    
       if not os.path.exists("./Utilities/DB/refseq_protein.pal") or not os.path.exists("./Utilities/DB/refseq_protein.%02d.psq"%(int(numPartitions))) or not os.path.exists("./Utilities/DB/allprots.faa"):
           if "phmmer" in packagesToInstall:
               dl = 'y'
           else:
              print "refseq protein DB not found or incomplete, needed for Annotate step, download now?"
              dl = raw_input("Enter Y/N: ")
           if dl == 'y' or dl == 'Y':
               print "Download and install refseq protein DB.."
               os.system("perl ./Utilities/perl/update_blastdb.pl refseq_protein")
               os.system("mv refseq_protein.*.tar.gz ./Utilities/DB/")
           
               fileList = glob.glob("./Utilities/DB/refseq_protein.*.tar.gz") 
               for file in fileList:
                  os.system("tar -C ./Utilities/DB/ -xvf %s"%(file))
               print "    running fastacmd (might take a few min)..."
               os.system(".%sUtilities%scpp%s%s-%s%sfastacmd -d ./Utilities/DB/refseq_protein -p T -a T -D 1 -o ./Utilities/DB/allprots.faa"%(os.sep, os.sep, os.sep, OSTYPE, MACHINETYPE, os.sep))

# sra toolkit
if not os.path.exists("./Utilities/cpp%s%s-%s%ssra"%(os.sep, OSTYPE, MACHINETYPE, os.sep)):
    sra = utils.getFromPath("srapath", "SRA PATH", False)
    if sra == "":
       if "sra" in packagesToInstall:
          dl = 'y'
       else:
          print "SRA binaries not found, optional for initPipeline step, download now?"
          dl = raw_input("Enter Y/N: ")
       if dl == 'y' or dl == 'Y':
           if OSTYPE == 'Linux' and MACHINETYPE == "x86_64":
              os.system("curl -L http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.3.3-3/sratoolkit.2.3.3-3-centos_linux64.tar.gz -o sra.tar.gz")
           elif  OSTYPE == "Darwin" and MACHINETYPE == "x86_64":
               os.system("curl -L http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.3.3-3/sratoolkit.2.3.3-3-mac64.tar.gz -o sra.tar.gz")
           os.system("tar xvzf sra.tar.gz")
           os.system("mv sratoolkit.2.3.3-3-* ./Utilities/cpp%s%s-%s%ssra"%(os.sep, OSTYPE, MACHINETYPE, os.sep)) 
           os.system("rm -rf sra.tar.gz")

if "isolate" in enabledWorkflows or "imetamos" in enabledWorkflows or manual:
    # check for cmake

    if not os.path.exists("./Utilities/cpp%s%s-%s%scmake"%(os.sep, OSTYPE, MACHINETYPE, os.sep)):
       cmake = utils.getFromPath("cmake", "CMAKE", False)
       if cmake == "":
          if "cmake" in packagesToInstall:
            dl = 'y'
          else:
             print "cmake binaries not found, optional for initPipeline step, download now?"
             dl = raw_input("Enter Y/N: ")
          if dl == 'y' or dl == 'Y':
             os.system("curl -L http://www.cmake.org/files/v2.8/cmake-2.8.12.tar.gz -o cmake.tar.gz")
             os.system("tar xvzf cmake.tar.gz")
             os.system("mv cmake-2.8.12 ./Utilities/cpp%s%s-%s%scmake"%(os.sep, OSTYPE, MACHINETYPE, os.sep))
             os.chdir("./Utilities/cpp%s%s-%s%scmake"%(os.sep, OSTYPE, MACHINETYPE, os.sep))
             os.system("./bootstrap --prefix=`pwd`/build;make;make install")
             os.chdir("%s"%(METAMOS_ROOT))
             os.system("rm cmake.tar.gz")
    if os.path.exists("./Utilities/cpp%s%s-%s%scmake"%(os.sep, OSTYPE, MACHINETYPE, os.sep)):
       pathUpdate = "%s/Utilities/cpp%s%s-%s%scmake/build/bin"%(METAMOS_ROOT, os.sep, OSTYPE, MACHINETYPE, os.sep)

       if "PATH" in os.environ:
          pathUpdate = "%s%s%s"%(os.environ["PATH"], os.pathsep, pathUpdate)
       os.environ["PATH"]=pathUpdate
    os.chdir("%s"%(METAMOS_ROOT))

    if not os.path.exists("./CA"):
      if "ca" in packagesToInstall:
         dl = 'y'
      else:
         print "Celera Assembler binaries not found, optional for Assemble step, download now?"
         dl = raw_input("Enter Y/N: ")
      if dl == 'y' or dl == 'Y':
          os.system("curl -L https://downloads.sourceforge.net/project/wgs-assembler/wgs-assembler/wgs-8.1/wgs-8.1.tar.bz2 -o wgs-8.1.tar.bz2")
          os.system("tar xvjf wgs-8.1.tar.bz2")
          os.system("rm -rf wgs-8.1.tar.bz2")
          os.system("mv wgs-8.1 CA")
          # patch CA to support PacBio sequences and non-apple compilers on OSX
          if not ALLOW_FAST:
             os.system("cd CA/kmer/ && cp configure.sh configure.original")
             os.system("cd CA/kmer/ && cat configure.original |sed s/\-fast//g > configure.sh")
             os.system("cd CA/src/ && cp c_make.as c_make.original")
             os.system("cd CA/src/ && cat c_make.original |sed s/\-fast//g > c_make.as")
          if not HAVE_GCC42:
             os.system("cd CA/src/ && cp c_make.as c_make.original")
             os.system("cd CA/src/ && cat c_make.original |sed s/\-4.2//g > c_make.as")
          if GCC_VERSION >= 4.7:
             os.system("cd CA/src/ && cp c_make.as c_make.original")
             os.system("cd CA/src/ && cat c_make.original |sed s/\-rdynamic//g > c_make.as")
          updateMakeFileForDarwin("CA/kmer/Makefile", addedCFlags, addedLDFlags)
          updateMakeFileForDarwin("CA/samtools/Makefile", addedCFlags, addedLDFlags)
          updateMakeFileForDarwin("CA/src/c_make.as", addedCFlags, addedLDFlags)
          os.system("cd CA/samtools && make")
          os.system("cd CA/kmer && ./configure.sh && gmake install")
          os.system("cd CA/src && gmake")

    if not os.path.exists("./Utilities/cpp%s%s-%s%sRay"%(os.sep, OSTYPE, MACHINETYPE, os.sep)):
       if "ray" in packagesToInstall:
          dl = 'y'
       else:
          print "Ray binaries not found, optional for Assemble step, download now?"
          dl = raw_input("Enter Y/N: ")
       if dl == 'y' or dl == 'Y':
          # check for mpi which is required
          command="mpicxx"
          mpi=utils.getFromPath(command, "MPI", False)
          if not os.path.exists("%s%s%s"%(mpi, os.sep, command)):
             command="openmpicxx"
             mpi=utils.getFromPath(command, "MPI", False)
             if not os.path.exists("%s%s%s"%(mpi, os.sep, command)):
                mpi = command = ""
                print "Error: cannot find MPI, required to build Ray. Please add it to your path."
          if command != "":
             os.system("curl -L http://downloads.sourceforge.net/project/denovoassembler/Ray-v2.2.0.tar.bz2 -o Ray-v2.2.0.tar.bz2")
             os.system("tar xvjf Ray-v2.2.0.tar.bz2")
             os.system("mv Ray-v2.2.0 ./Utilities/cpp/%s%s-%s%sRay"%(os.sep, OSTYPE, MACHINETYPE, os.sep))
             os.chdir("./Utilities/cpp/%s%s-%s%sRay"%(os.sep, OSTYPE, MACHINETYPE, os.sep))
             os.system("make PREFIX=bin MPICXX=%s%s%s MAXKMERLENGTH=128 MPI_IO=y DEBUG=n ASSERT=n EXTRA=\" -march=native\""%(mpi, os.sep, command))
             os.system("make install")
             os.chdir("%s"%(METAMOS_ROOT))
             os.system("rm -rf Ray-v2.2.0.tar.bz2")

    if not os.path.exists("./Utilities/cpp%s%s-%s%skmergenie"%(os.sep, OSTYPE, MACHINETYPE, os.sep)):
        kmerGenie = utils.getFromPath("kmergenie", "Kmer Genie", False)
        if kmerGenie == "":
           if "kmergenie" in packagesToInstall:
              dl = 'y'
           else:
              print "Kmer Genie was not found, optional for Assemble step, download now?"
              dl = raw_input("Enter Y/N: ")
        if dl == 'y' or dl == 'Y':
           os.system("curl -L ftp://ftp.cbcb.umd.edu/pub/data/metamos/kmergenie-1.5692.tar.gz -o kmer.tar.gz")
           os.system("tar xvzf kmer.tar.gz")
           os.system("mv kmergenie-1.5692 ./Utilities/cpp%s%s-%s%skmergenie"%(os.sep, OSTYPE, MACHINETYPE, os.sep))
           os.chdir("./Utilities/cpp%s%s-%s%skmergenie"%(os.sep, OSTYPE, MACHINETYPE, os.sep))
           updateMakeFileForDarwin("makefile", addedCFlags, addedLDFlags)
           os.system("make k=300")
           os.chdir("%s"%(METAMOS_ROOT))
           os.system("rm -rf kmer.tar.gz")

    if not os.path.exists("./Utilities/cpp%s%s-%s%sspades"%(os.sep, OSTYPE, MACHINETYPE, os.sep)):
        spades = utils.getFromPath("spades.py", "SPAdes", False)
        if spades == "":
           if "spades" in packagesToInstall:
              dl = 'y'
           else:
              print "SPAdes was not found, optional for Assemble step, download now?"
              dl = raw_input("Enter Y/N: ")
        if dl == 'y' or dl == 'Y':
           if OSTYPE == "Darwin":
              if GCC_VERSION < 4.7:
                 print "Error: SPAdes requires gcc at least version 4.7, found version %s. Please update and try again"%(GCC_VERSION)
              else:
                 os.system("curl -L http://spades.bioinf.spbau.ru/release3.0.0/SPAdes-3.0.0.tar.gz -o spades.tar.gz")
                 os.system("tar xvzf spades.tar.gz")
                 os.system("mv SPAdes-3.0.0 ./Utilities/cpp%s%s-%s%sspades"%(os.sep, OSTYPE, MACHINETYPE, os.sep))
                 os.chdir("./Utilities/cpp%s%s-%s%sspades"%(os.sep, OSTYPE, MACHINETYPE, os.sep))
                 os.system("export CC=`which gcc` && bash spades_compile.sh")
                 os.chdir("%s"%(METAMOS_ROOT))
           else:
              os.system("curl -L http://spades.bioinf.spbau.ru/release3.0.0/SPAdes-3.0.0-Linux.tar.gz -o spades.tar.gz")
              os.system("tar xvzf spades.tar.gz")
              os.system("mv SPAdes-3.0.0-Linux ./Utilities/cpp%s%s-%s%sspades"%(os.sep, OSTYPE, MACHINETYPE, os.sep))
           os.system("rm -rf spades.tar.gz")

    if not os.path.exists("./Utilities/cpp%s%s-%s%sprokka"%(os.sep, OSTYPE, MACHINETYPE, os.sep)):
       prokaBin = utils.getFromPath("prokka", "Prokka", False)
       dl = 'n'
       if prokaBin == "":
          if "prokka" in packagesToInstall:
             dl = 'y'
          else:
             print "Prokka binaries not found, optional for Assemble step, download now?"
             dl = raw_input("Enter Y/N: ")
       if dl == 'y' or dl == 'Y':
          signalp = utils.getFromPath("signalp", "SignalP", False)
          if signalp == "":
             print "Warning: SignalP is not installed and is required for Prokka's gram option. Please download it and add it to your path."
          os.system("curl -L ftp://ftp.cbcb.umd.edu/pub/data/metamos/prokka-1.7.tar.gz -o prokka-1.7.tar.gz")
          os.system("tar xvzf prokka-1.7.tar.gz")
          os.system("mv prokka-1.7 ./Utilities/cpp%s%s-%s%sprokka"%(os.sep, OSTYPE, MACHINETYPE, os.sep))
          os.system("rm prokka-1.7.tar.gz")

          bioperl = utils.getCommandOutput("perl -MBio::Seq -e 0 && echo $?", True)
          perltime = utils.getCommandOutput("perl -MTime::Piece -e 0 && echo $?", True)
          xmlsimple = utils.getCommandOutput("perl -MXML::Simple -e 0 && echo $?", True)
          storable = utils.getCommandOutput("perl -MStorable -e 0 && echo $?", True) 
          xmlparser = utils.getCommandOutput("perl -MXML:Parser -e 0 && echo $?", True)

          # always install bioperl, otherwise parts may be missing or it may be the wrong version
          # phylosift comes with BioPerl, use it
          os.system("curl -L http://edhar.genomecenter.ucdavis.edu/~koadman/phylosift/devel/phylosift_20130829.tar.bz2 -o ./phylosift.tar.bz2")
          os.system("tar -xvjf phylosift.tar.bz2")
          os.system("rm -rf phylosift.tar.bz2")
          os.system("mv phylosift_20130829/lib ./Utilities/cpp%s%s-%s%sprokka"%(os.sep, OSTYPE, MACHINETYPE, os.sep))
          os.system("rm -rf phylosift_20130829")

          if perltime == "":
             os.system("curl -L http://search.cpan.org/CPAN/authors/id/M/MS/MSERGEANT/Time-Piece-1.08.tar.gz -o time.tar.gz")
             os.system("tar -xvzf time.tar.gz")
             os.chdir("Time-Piece-1.08")
             os.system("perl Makefile.PL PREFIX=`pwd`/build")
             os.system("make install")
             os.chdir("%s"%(METAMOS_ROOT))
             pathToCopy = utils.getCommandOutput("find Time-Piece-1.08/build -type d -name \"Time\" |grep -v auto", False)
             copyPerlLib(pathToCopy, "./Utilities/cpp%s%s-%s%sprokka/lib/"%(os.sep, OSTYPE, MACHINETYPE, os.sep))
             os.system("rm -rf time.tar.gz")
             os.system("rm -rf Time-Piece-1.08") 

          if xmlparser == "":
             os.system("curl -L http://search.cpan.org/CPAN/authors/id/M/MS/MSERGEANT/XML-Parser-2.36.tar.gz -o parse.tar.gz")
             os.system("tar -xvzf parse.tar.gz")
             os.chdir("XML-Parser-2.36")
             os.system("perl Makefile.PL PREFIX=`pwd`/build")
             os.system("make install")
             os.chdir("%s"%(METAMOS_ROOT))
             pathToCopy = utils.getCommandOutput("find XML-Parser-2.36/build -type d -name \"XML\" |grep -v auto", False)
             copyPerlLib(pathToCopy, "./Utilities/cpp%s%s-%s%sprokka/lib/"%(os.sep, OSTYPE, MACHINETYPE, os.sep))
             libUpdate = "%s/Utilities/cpp%s%s-%s%sprokka/lib/"%(METAMOS_ROOT, os.sep, OSTYPE, MACHINETYPE, os.sep)
             if "PERL5LIB" in os.environ:
                libUpdate = "%s%s%s"%(os.environ["PERL5LIB"], os.pathsep, libUpdate)
             os.environ["PERL5LIB"]=libUpdate
             os.system("rm -rf parse.tar.gz")
             os.system("rm -rf XML-Parser-2.36")

          if xmlsimple == "":
             os.system("curl -L http://search.cpan.org/CPAN/authors/id/G/GR/GRANTM/XML-Simple-1.08.tar.gz -o xml.tar.gz")
             os.system("tar -xvzf xml.tar.gz")
             os.chdir("XML-Simple-1.08")
             os.system("perl Makefile.PL PREFIX=`pwd`/build")
             os.system("make install")
             os.chdir("%s"%(METAMOS_ROOT))
             pathToCopy = utils.getCommandOutput("find XML-Simple-1.08/build -type d -name \"XML\" |grep -v auto", False)
             copyPerlLib(pathToCopy, "./Utilities/cpp%s%s-%s%sprokka/lib/"%(os.sep, OSTYPE, MACHINETYPE, os.sep))
             os.system("rm -rf xml.tar.gz")
             os.system("rm -rf XML-Simple-1.08")

          if os.path.exists("./Utilities/cpp%s%s-%s%sprokka/lib"%(os.sep, OSTYPE, MACHINETYPE, os.sep)):
             os.chdir("./Utilities/cpp%s%s-%s%sprokka/bin"%(os.sep, OSTYPE, MACHINETYPE, os.sep))
             os.system("cp prokka prokka.original")
             os.system("cat prokka.original |awk '{if (match($0, \"use strict\")) { print \"use lib \\\"%s/Utilities/cpp%s%s-%s%sprokka/lib\\\";\"; print $0; } else { print $0}}' > prokka"%(METAMOS_ROOT, os.sep, OSTYPE, MACHINETYPE, os.sep))
             os.chdir("%s"%(METAMOS_ROOT))
          # for some reason prokka adds its binaries to the end of path, not beginning so if your path has the wrong version of a program, it will crash. Update
          os.chdir("./Utilities/cpp%s%s-%s%sprokka/bin"%(os.sep, OSTYPE, MACHINETYPE, os.sep))
          os.system("cp prokka prokka.original")
          os.system("cat prokka.original |awk '{if (match($0, \"ENV{PATH}\")) { print \"$ENV{PATH} = $BINDIR . \\\":\\\" . $ENV{PATH};\"; } else { print $0}}' > prokka")
          os.chdir("%s"%(METAMOS_ROOT))

          aragorn = utils.getFromPath("aragorn", "aragorn", False)
          aragornVersion = ""
          if aragorn != "":
             aragornVersion = utils.getCommandOutput("%s/aragorn -h 2>&1 | grep -i '^ARAGORN v' |sed s/v//g |awk '{printf(\"%%2.2f\n\", $2)}'", True)
             if float(aragornVersion) < 1.2:
                aragorn = ""
          if aragorn == "":
             print "Aragorn missing, will install"
             os.system("curl -L http://130.235.46.10/ARAGORN/Downloads/aragorn1.2.36.tgz -o aragorn.tar.gz")
             os.system("tar xvzf aragorn.tar.gz")
             os.chdir("aragorn1.2.36")
             os.system("gcc -O3 -ffast-math -finline-functions %s %s -o aragorn aragorn1.2.36.c"%(addedCFlags, addedLDFlags))
             os.chdir("%s"%(METAMOS_ROOT))
             os.system("mv aragorn1.2.36/aragorn ./Utilities/cpp%s%s-%s%sprokka/binaries%s%s"%(os.sep, OSTYPE, MACHINETYPE, os.sep, os.sep, OSTYPE.lower()))
             os.system("rm -rf aragorn1.2.36")
             os.system("rm aragorn.tar.gz")
          infernal = utils.getFromPath("cmscan", "Infernal", False)
          if infernal == "" and not os.path.exists("./Utilities/cpp%s%s-%s%sprokka/binaries/%s/infernal"%(os.sep, OSTYPE, MACHINETYPE, os.sep, OSTYPE.lower())):
             print "Infernal missing, will install"
             if OSTYPE == "Darwin":
                os.system("curl -L http://selab.janelia.org/software/infernal/infernal-1.1rc4-macosx-intel.tar.gz -o infernal.tar.gz")
             else:
                os.system("curl -L http://selab.janelia.org/software/infernal/infernal-1.1rc4-linux-intel-gcc.tar.gz -o infernal.tar.gz")
             os.system("tar xvzf infernal.tar.gz")
             os.system("mv infernal*/binaries/* ./Utilities/cpp%s%s-%s%sprokka/binaries%s%s"%(os.sep, OSTYPE, MACHINETYPE, os.sep, os.sep, OSTYPE.lower()))
             os.system("rm -rf infernal*")
          barrnap = utils.getFromPath("barrnap", "barrnap", False)
          if barrnap == "" and not os.path.exists("./Utilities/cpp%s%s-%s%sprokka/binaries/%s/barrnap"%(os.sep, OSTYPE, MACHINETYPE, os.sep, OSTYPE.lower())):
             print "Barrnap missing, will install"
             os.system("curl -L http://www.vicbioinformatics.com/barrnap-0.1.tar.gz -o barrnap.tar.gz")
             os.system("tar xvzf barrnap.tar.gz")
             os.system("mv barrnap-0.1/barrnap ./Utilities/cpp%s%s-%s%sprokka/binaries%s%s"%(os.sep, OSTYPE, MACHINETYPE, os.sep, os.sep, OSTYPE.lower()))
             os.system("mv barrnap-0.1/db ./Utilities/cpp%s%s-%s%sprokka/binaries%s%s"%(os.sep, OSTYPE, MACHINETYPE, os.sep, os.sep, OSTYPE.lower()))
             if os.path.exists("./Utilities/cpp%s%s-%s%sprokka/lib"%(os.sep, OSTYPE, MACHINETYPE, os.sep)):
                os.chdir("./Utilities/cpp%s%s-%s%sprokka/binaries%s%s"%(os.sep, OSTYPE, MACHINETYPE, os.sep, os.sep, OSTYPE.lower()))
                os.system("cp barrnap barrnap.original")
                os.system("cat barrnap.original |awk '{if (match($0, \"use strict\")) { print \"use lib \\\"%s/Utilities/cpp%s%s-%s%sprokka/lib\\\";\"; print $0; } else { print $0}}' > barrnap"%(METAMOS_ROOT, os.sep, OSTYPE, MACHINETYPE, os.sep))
                os.chdir("%s"%(METAMOS_ROOT))
             os.system("rm -rf barrnap-0.1")
             os.system("rm barrnap.tar.gz")

          hmmscan = utils.getFromPath("hmmscan", "HMMER3", False)
          hmmscanVersion = ""
          if hmmscan != "":
             hmmscanVersion = utils.getCommandOutput("%s/hmmscan -h | grep '^# HMMER' |awk '{printf(\"%%2.2f\\n\", $3)}'"%(hmmscan), True)
             print "Found HMM SCAN %s %s"%(hmmscan, hmmscanVersion)
             if float(hmmscanVersion) < 3.1:
                hmmscan = ""
          if hmmscan == "" and not os.path.exists("./Utilities/cpp%s%s-%s%sprokka/binaries/%s/hmmscan"%(os.sep, OSTYPE, MACHINETYPE, os.sep, OSTYPE.lower())):
             print "HMMER3 is missing, will install"
             if OSTYPE == "Darwin":
                os.system("curl -L ftp://selab.janelia.org/pub/software/hmmer3/3.1b1/hmmer-3.1b1-macosx-intel.tar.gz -o hmmer.tar.gz")
             elif OSTYPE == "Linux" and MACHINETYPE == "x86_64":
                os.system("curl -L ftp://selab.janelia.org/pub/software/hmmer3/3.1b1/hmmer-3.1b1-linux-intel-x86_64.tar.gz -o hmmer.tar.gz")
             elif OSTYPE == "Linux":
                os.system("curl -L ftp://selab.janelia.org/pub/software/hmmer3/3.1b1/hmmer-3.1b1-linux-intel-ia32.tar.gz -o hmmer.tar.gz") 
             os.system("tar xvzf hmmer.tar.gz")
             os.system("mv hmmer*/binaries/* ./Utilities/cpp%s%s-%s%sprokka/binaries%s%s"%(os.sep, OSTYPE, MACHINETYPE, os.sep, os.sep, OSTYPE.lower()))
             os.system("rm -rf hmmer*")

          gnuparallel = utils.getFromPath("parallel", "GNU Parallel", False)
          if gnuparallel == "" and not os.path.exists("./Utilities/cpp%s%s-%s%sprokka/binaries/%s/parallel"%(os.sep, OSTYPE, MACHINETYPE, os.sep, OSTYPE.lower())):
             print "GNU Parallel is missing, will install"
             os.system("curl -L http://ftp.gnu.org/gnu/parallel/parallel-20130822.tar.bz2 -o parallel.tar.gz")
             os.system("tar xvjf parallel.tar.gz")
             os.chdir("parallel-20130822")
             os.system("./configure --prefix=`pwd`")
             os.system("make install")
             os.chdir("%s"%(METAMOS_ROOT))
             os.system("mv parallel-20130822/bin/parallel ./Utilities/cpp%s%s-%s%sprokka/binaries%s%s"%(os.sep, OSTYPE, MACHINETYPE, os.sep, os.sep, OSTYPE.lower()))
             os.system("rm -rf parallel-20130822")
             os.system("rm parallel.tar.gz")
          
          blastp = utils.getFromPath("blastp", "BLAST+", False)
          if blastp == "" and not os.path.exists("./Utilities/cpp%s%s-%s%sprokka/binaries/%s/blastp"%(os.sep, OSTYPE, MACHINETYPE, os.sep, OSTYPE.lower())):
             os.system("ln %s/Utilities/cpp%s%s-%s%sblastp %s/Utilities/cpp%s%s-%s%sprokka/binaries%s%s%sblastp"%(METAMOS_ROOT, os.sep, OSTYPE, MACHINETYPE, os.sep, METAMOS_ROOT, os.sep, OSTYPE, MACHINETYPE, os.sep, os.sep, OSTYPE.lower(), os.sep))

          prodigal = utils.getFromPath("prodigal", "PRODIGAL", False)
          if prodigal != "":
             prodigalVersion = utils.getCommandOutput("%s/prodigal -v 2>&1 | grep -i '^Prodigal V' |sed s/V//g |awk '{printf(\"%%2.2f\\n\", $2)}'"%(prodigal), True)
             print "Found prodigal %s %s"%(prodigal, prodigalVersion)
             if float(prodigalVersion) < 2.6:
                prodigal = ""

          if prodigal == "":
             os.system("curl -L https://prodigal.googlecode.com/files/Prodigal-2.60.tar.gz -o prodigal.tar.gz")
             os.system("tar xvzf prodigal.tar.gz")
             os.system("mv Prodigal-2.60 prodigal.v2_60") 
             os.system("rm -rf prodigal.tar.gz")
             os.system("curl -L https://prodigal.googlecode.com/files/prodigal_v2.60.bugfix1.tar.gz -o prodigal.tar.gz")
             os.system("tar xvzf prodigal.tar.gz")
             os.system("mv prodigal_v2.60.bugfix1/* prodigal.v2_60/")
             os.chdir("prodigal.v2_60")
             updateMakeFileForDarwin("Makefile", addedCFlags, addedLDFlags)
             os.system("make")
             os.chdir("%s"%(METAMOS_ROOT))
             os.system("mv prodigal.v2_60/prodigal ./Utilities/cpp%s%s-%s%sprokka/binaries%s%s"%(os.sep, OSTYPE, MACHINETYPE, os.sep, os.sep, OSTYPE.lower()))
             os.system("rm -rf prodigal.tar.gz")
             os.system("rm -rf prodigal.v2_60") 
             os.system("rm -rf prodigal_v2.60.bugfix1")

          tbl2asn = utils.getFromPath("tbl2asn", "NCBI Tbl2Asn", False)
          if tbl2asn == "" and not os.path.exists("./Utilities/cpp%s%s-%s%sprokka/binaries/%s/tbl2asn"%(os.sep, OSTYPE, MACHINETYPE, os.sep, OSTYPE.lower())):
             print "NCBI Tbl2Asn is missing, will install"
             if OSTYPE == "Darwin":
                os.system("curl -L ftp://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/tbl2asn/mac.tbl2asn.gz -o tbl2asn.gz")
             elif OSTYPE == "Linux" and MACHINETYPE == "x86_64":
                os.system("curl -L ftp://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/tbl2asn/linux64.tbl2asn.gz -o tbl2asn.gz")
             elif OSTYPE == "Linux":
                os.system("curl -L ftp://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/tbl2asn/linux.tbl2asn.gz -o tbl2asn.gz")
             os.system("gunzip tbl2asn.gz")
             os.system("chmod ug+x tbl2asn")
             os.system("mv tbl2asn ./Utilities/cpp%s%s-%s%sprokka/binaries%s%s"%(os.sep, OSTYPE, MACHINETYPE, os.sep, os.sep, OSTYPE.lower()))
             os.system("rm tbl2asn.gz")

    if not os.path.exists("./Utilities/cpp%s%s-%s%ssoap2"%(os.sep, OSTYPE, MACHINETYPE, os.sep)):
       if "soap2" in packagesToInstall:
          if "soap2" in packagesToInstall:
             dl = 'y'
          else:
             print "SOAPdenovo2 binaries not found, optional for Assemble step, download now?"
             dl = raw_input("Enter Y/N: ")
          if dl == 'y' or dl == 'Y':
             if OSTYPE == "Darwin":
                os.system("curl -L http://sourceforge.net/projects/soapdenovo2/files/SOAPdenovo2/src/r240/SOAPdenovo2-src-r240-mac.tgz -o soap2.tar.gz")
                os.system("tar xvzf soap2.tar.gz")
                os.system("mv r240_noAIOinPregraph ./Utilities/cpp%s%s-%s%ssoap2"%(os.sep, OSTYPE, MACHINETYPE, os.sep))
                os.chdir("./Utilities/cpp%s%s-%s%ssoap2"%(os.sep, OSTYPE, MACHINETYPE, os.sep))
                updateMakeFileForDarwin("Makefile", addedCFlags, addedLDFlags)
                os.system("make clean")
                os.system("make")
                os.system("mkdir bin")
                os.system("mv SOAPdenovo-* bin/")
                os.chdir("%s"%(METAMOS_ROOT))
             else:
                os.system("curl -L http://sourceforge.net/projects/soapdenovo2/files/SOAPdenovo2/bin/r240/SOAPdenovo2-bin-LINUX-generic-r240.tgz -o soap2.tar.gz")
                os.system("tar xvzf soap2.tar.gz")
                os.system("mkdir ./Utilities/cpp%s%s-%s%ssoap2"%(os.sep, OSTYPE, MACHINETYPE, os.sep))
                os.system("mv SOAPdenovo2-bin-LINUX-generic-r240 ./Utilities/cpp%s%s-%s%ssoap2/bin"%(os.sep, OSTYPE, MACHINETYPE, os.sep))
             os.system("curl -L http://sourceforge.net/projects/soapdenovo2/files/GapCloser/src/r6/GapCloser-src-v1.12-r6.tgz -o gapcloser.tar.gz")
             os.system("tar xvzf gapcloser.tar.gz")
             os.chdir("v1.12-r6")
             updateMakeFileForDarwin("Makefile", addedCFlags, addedLDFlags)
             os.system("make")
             os.chdir("%s"%(METAMOS_ROOT))
             os.system("mv v1.12-r6/Release/* ./Utilities/cpp%s%s-%s%ssoap2/bin"%(os.sep, OSTYPE, MACHINETYPE, os.sep))
             os.system("rm -rf soap2.tar.gz")
             os.system("rm -rf v1.12-r6")
             os.system("rm -rf gapcloser.tar.gz")

    if not os.path.exists("./Utilities/cpp%s%s-%s%sMaSuRCA"%(os.sep, OSTYPE, MACHINETYPE, os.sep)):
       masurca = utils.getFromPath("masurca", "MaSuRCA", False)
       if masurca == "":
          if "masurca" in packagesToInstall:
             dl = 'y'
          else:
             print "MaSuRCA binaries not found, optional for Assemble step, download now?"
             dl = raw_input("Enter Y/N: ")
          if dl == 'y' or dl == 'Y':
             if GCC_VERSION < 4.4:
                print "Error: MaSuRCA requires gcc at least version 4.4, found version %s. Please update and try again"%(GCC_VERSION)
             else:
                os.system("curl -L ftp://ftp.cbcb.umd.edu/pub/data/metamos/MaSuRCA-2.2.0.tar.gz -o msrca.tar.gz")
                os.system("tar xvzf msrca.tar.gz")
                os.system("mv ./MaSuRCA-2.2.0 ./Utilities/cpp%s%s-%s%sMaSuRCA"%(os.sep, OSTYPE, MACHINETYPE, os.sep))
                os.chdir("./Utilities/cpp%s%s-%s%sMaSuRCA"%(os.sep, OSTYPE, MACHINETYPE, os.sep))
                os.system("cp install.sh install.orig")
                os.system("cat install.orig |sed s/\-\-prefix/\-\-disable\-shared\ \-\-prefix/g > install.sh")
                # patch CA
                if not ALLOW_FAST:
                   os.system("cd CA/kmer/ && cp configure.sh configure.original")
                   os.system("cd CA/kmer/ && cat configure.original |sed s/\-fast//g > configure.sh")
                if not HAVE_RT:
                   os.system("cd SuperReads/ && cp Makefile.am Makefile.am.original")
                   os.system("cd SuperReads/ && cat Makefile.am.original |sed s/\-lrt//g > Makefile.am")
                   os.system("cd SuperReads/ && cp Makefile.in Makefile.in.original")
                   os.system("cd SuperReads/ && cat Makefile.in.original |sed s/\-lrt//g > Makefile.in")
                if not HAVE_QUIET_HEAD:
                   os.system("cd SuperReads/src && cp masurca masurca.original")
                   os.system("cd SuperReads/src && cat masurca.original |sed s/head\ \-q/head/g > masurca")
                os.system("rm -f CA/kmer/makepath")

                # fix compilation on OSX
                if OSTYPE == "Darwin":
                   os.system("cp SuperReads/include/reallocators.hpp SuperReads/include/reallocators.hpp.orig")
                   testIn = open("SuperReads/include/reallocators.hpp.orig", 'r')
                   testOut = open("SuperReads/include/reallocators.hpp", 'w')
                   for line in testIn.xreadlines():
                      if "T* res = reallocator<T>::operator()(" in line:
                         testOut.write("T* res = reallocator<T>::realloc(ptr, osize, nsize);\n")
                      else:
                         testOut.write(line.strip() + "\n")
                   testIn.close()
                   testOut.close()

                   # dont set static building libs on OSX, sseems to cause compile issues for jellyfish
                   os.environ["CFLAGS"] = oldCFlags  
                   os.environ["CPPFLAGS"] = oldCPPFlags
                   os.environ["CXXFLAGS"] = oldCXXFlags
                   os.environ["LDFLAGS"] = oldLDFlags

                updateMakeFileForDarwin("CA/src/c_make.as", addedCFlags, addedLDFlags)
                os.system("bash install.sh")
                fileOptions = utils.getCommandOutput("file -b --mime-type INSTALL.py", False)
                if fileOptions == "":
                   fileOptions = utils.getCommandOutput("file -b --mime INSTALL.py", False)
                   if fileOptions != "":
                      # fix file command used by MaSuRCA, its not compatible with the system
                      if os.path.exists("bin/expand_fastq"):
                         os.system("cp bin/expand_fastq bin/expand_fastq.orig")
                         testIn = open("bin/expand_fastq.orig", 'r')
                         testOut = open("bin/expand_fastq", 'w')
                         for line in testIn.xreadlines():
                            if "case $(file" in line:
                               testOut.write("case $(file -b --mime \"$FILE\" |awk '{print $1}'|sed s/\\;//g) in\n")
                            else:
                               testOut.write(line.strip() + "\n")
                         testIn.close()
                         testOut.close()
                   else:
                      os.chdir("%s"%(METAMOS_ROOT))
                      os.system("rm -rf ./Utilities/cpp%s%s-%s%sMaSuRCA"%(os.sep, OSTYPE, MACHINETYPE, os.sep))
                # update path to CA which is always hardcoded to Linux-amd64
                os.system("cp bin/masurca bin/masurca.orig")
                os.system("cat bin/masurca.orig | sed s/Linux-amd64/%s-%s/g |sed s/check_exec\\(\\\"jellyfish\\\"/check_exec\\(\\\"jellyfish-2.0\\\"/g > bin/masurca"%(OSTYPE, MACHINETYPE.replace("x86_64", "amd64")))

                if OSTYPE == "Darwin":
                   os.system("cp bin/masurca bin/masurca.orig")
                   os.system("cat bin/masurca.orig | awk '{if (match($0, \"save NUM_SUPER_READS\")) { print $0\"\\n\\tprint FILE \\\"export NUM_SUPER_READS=\\\\$NUM_SUPER_READS\\\\n\\\";\"; } else { print $0}}' | sed  s/\\(\\'..TOTAL_READS\\'/\\(\\\\\\\\\\$ENV{\\'TOTAL_READS\\'}/g| sed s/'<..$NUM_SUPER_READS.'/\"<ENVIRON[\\\\\\\\\\\"NUM_SUPER_READS\\\\\\\\\\\"]\"/g | sed s/'>=..$NUM_SUPER_READS.'/\">=ENVIRON[\\\\\\\\\\\"NUM_SUPER_READS\\\\\\\\\\\"]\"/g > bin/masurca")

                   # reset env variables again
                   addEnvironmentVar("CFLAGS", " %s "%(addedCFlags))
                   addEnvironmentVar("CPPFLAGS", " %s "%(addedCFlags))
                   addEnvironmentVar("CXXFLAGS", " %s "%(addedCFlags))
                   addEnvironmentVar("LDFLAGS", " %s "%(addedLDFlags))

                os.chdir("%s"%(METAMOS_ROOT))
                os.system("rm -rf ./MaSuRCA-2.2.0")
                os.system("rm msrca.tar.gz")

    if not os.path.exists("./Utilities/cpp%s%s-%s%smira"%(os.sep, OSTYPE, MACHINETYPE, os.sep)):
       mira = utils.getFromPath("mira", "MIRA", False)
       if mira == "":
          if "mira" in packagesToInstall:
             dl = 'y'
          else:
             print "MIRA binaries not found, optional for Assemble step, download now?"
             dl = raw_input("Enter Y/N: ")
          if dl == 'y' or dl == 'Y':
             if OSTYPE == "Darwin":
	        os.system("curl -L ftp://ftp.cbcb.umd.edu/pub/data/metamos/mira_4.0rc5_darwin13.0.0_x86_64_static.tar.bz2 -o mira.tar.bz2")
             else:
                os.system("curl -L ftp://ftp.cbcb.umd.edu/pub/data/metamos/mira_4.0rc5_linux-gnu_x86_64_static.tar.bz2 -o mira.tar.bz2")
             os.system("tar xvjf mira.tar.bz2")
             os.system("rm -f mira.tar.bz2")
             os.system("mv `ls -d mira*` ./Utilities/cpp%s%s-%s%smira"%(os.sep, OSTYPE, MACHINETYPE, os.sep))

    if not os.path.exists("./Utilities/cpp%s%s-%s%sidba"%(os.sep, OSTYPE, MACHINETYPE, os.sep)):
       idba = utils.getFromPath("idba", "IDBA-UD", False)
       if idba == "":
          if "idba" in packagesToInstall:
             dl = 'y'
          else:
             print "IDBA-UD binaries not found, optional for Assemble step, download now?"
             dl = raw_input("Enter Y/N: ")
          if dl == 'y' or dl == 'Y':
             os.system("curl -L http://hku-idba.googlecode.com/files/idba-1.1.1.tar.gz -o idba.tar.gz")
             os.system("tar xvzf idba.tar.gz")
             os.system("mv idba-1.1.1 ./Utilities/cpp%s%s-%s%sidba"%(os.sep, OSTYPE, MACHINETYPE, os.sep))
             os.chdir("./Utilities/cpp%s%s-%s%sidba"%(os.sep, OSTYPE, MACHINETYPE, os.sep))
             os.system("mv src/sequence/short_sequence.h src/sequence/short_sequence.orig")
             os.system("cat src/sequence/short_sequence.orig |awk '{if (match($0, \"kMaxShortSequence = 128\")) print \"static const uint32_t kMaxShortSequence = 32768;\"; else print $0}' > src/sequence/short_sequence.h")
             os.system("./configure")
             os.system("make")
             os.chdir("%s"%(METAMOS_ROOT))
             os.system("rm -rf idba.tar.gz")

    if not os.path.exists("./Utilities/cpp%s%s-%s%seautils"%(os.sep, OSTYPE, MACHINETYPE, os.sep)):
       eautils = utils.getFromPath("fastq-mcf", "EA-UTILS", False)
       if eautils == "":
          if "eautils" in packagesToInstall:
             dl = 'y'
          else:
             print "EA-UTILS binaries not found, optional for Assemble step, download now?"
             dl = raw_input("Enter Y/N: ")
          if dl == 'y' or dl == 'Y':
             os.system("curl -L https://ea-utils.googlecode.com/files/ea-utils.1.1.2-537.tar.gz -o eautils.tar.gz")
             os.system("curl -L ftp://ftp.gnu.org/gnu/gsl/gsl-1.16.tar.gz -o gsl.tar.gz")
             os.system("tar xvzf eautils.tar.gz")
             os.system("tar xvzf gsl.tar.gz")
             os.system("mv gsl-1.16 ea-utils.1.1.2-537/gsl")
             os.system("rm ea-utils.1.1.2-537/tidx/utils.cpp")
             os.system("mv ea-utils.1.1.2-537 ./Utilities/cpp%s%s-%s%seautils"%(os.sep, OSTYPE, MACHINETYPE, os.sep))
             os.chdir("./Utilities/cpp%s%s-%s%seautils/gsl"%(os.sep, OSTYPE, MACHINETYPE, os.sep))
             os.system("./configure --prefix=`pwd`/build")
             os.system("make")
             os.system("make install")
             os.chdir("..")
             os.system("mv Makefile Makefile.orig")
             os.system("cat Makefile.orig |sed s/CFLAGS?=/CFLAGS+=/g |sed s/CPPFLAGS?=/CPPFLAGS+=/g > Makefile")
             addEnvironmentVar("CFLAGS", "-L%s/Utilities/cpp%s%s-%s%seautils/gsl/build/lib/"%(METAMOS_ROOT, os.sep, OSTYPE, MACHINETYPE, os.sep))
             addEnvironmentVar("CPPFLAGS", "-L%s/Utilities/cpp%s%s-%s%seautils/gsl/build/lib/"%(METAMOS_ROOT, os.sep, OSTYPE, MACHINETYPE, os.sep))
             os.system("make")
             os.chdir("%s"%(METAMOS_ROOT))
             os.system("rm -rf eautils.tar.gz")
             os.system("rm -rf gsl.tar.gz")

    if not os.path.exists("./Utilities/cpp%s%s-%s%sabyss"%(os.sep, OSTYPE, MACHINETYPE, os.sep)):
       abyss = utils.getFromPath("ABYSS", "ABySS", False)
       if abyss == "":
          if "abyss" in packagesToInstall:
             dl = 'y'
          else:
             print "ABySS binaries not found, optional for Assemble step, download now?"
             dl = raw_input("Enter Y/N: ")
          if dl == 'y' or dl == 'Y':
             os.system("curl -L https://sparsehash.googlecode.com/files/sparsehash-2.0.2.tar.gz -o sparse.tar.gz")
             os.system("tar xvzf sparse.tar.gz")
             os.chdir("sparsehash-2.0.2")
             os.system("./configure --prefix=`pwd`")
             os.system("make install")
             os.chdir("%s"%(METAMOS_ROOT))
             os.system("curl -L http://sourceforge.net/projects/boost/files/boost/1.54.0/boost_1_54_0.tar.gz -o boost.tar.gz")
             os.system("tar xvzf boost.tar.gz")
             os.system("curl -L http://www.bcgsc.ca/platform/bioinfo/software/abyss/releases/1.3.6/abyss-1.3.6.tar.gz -o abyss.tar.gz")
             os.system("tar xvzf abyss.tar.gz")
             os.chdir("abyss-1.3.6")
             os.system("ln -s %s/boost_1_54_0/boost boost"%(METAMOS_ROOT))
             addEnvironmentVar("CFLAGS", "-I%s/sparsehash-2.0.2/include"%(METAMOS_ROOT))
             addEnvironmentVar("CPPFLAGS", "-I%s/sparsehash-2.0.2/include"%(METAMOS_ROOT))
             addEnvironmentVar("CXXFLAGS", "-I%s/sparsehash-2.0.2/include"%(METAMOS_ROOT))

             # sparse hash library has unused variables which cause warnings with gcc 4.8 so disable -Werror
             if GCC_VERSION >= 4.8:
                os.system("mv configure configure.original")
                os.system("cat configure.original |sed s/\-Werror//g > configure")
                os.system("chmod ug+x configure")

             os.system("./configure --enable-maxk=96 --prefix=`pwd`/build")
             os.system("make install")
             os.chdir("%s"%(METAMOS_ROOT))
             os.system("mkdir ./Utilities/cpp%s%s-%s%sabyss"%(os.sep, OSTYPE, MACHINETYPE, os.sep))
             os.system("mv abyss-1.3.6/build/* ./Utilities/cpp%s%s-%s%sabyss/"%(os.sep, OSTYPE, MACHINETYPE, os.sep))

             # update abysss to use installed mpi
             command="mpirun"
             mpi=utils.getFromPath(command, "MPI", False)
             if not os.path.exists("%s%s%s"%(mpi, os.sep, command)):
                command="openmpirun"
                mpi=utils.getFromPath(command, "MPI", False)
                if not os.path.exists("%s%s%s"%(mpi, os.sep, command)):
                   mpi = command = ""
             os.chdir("./Utilities/cpp%s%s-%s%sabyss/bin/"%(os.sep, OSTYPE, MACHINETYPE, os.sep))
             os.system("cp abyss-pe abyss-pe-orig")
             if mpi != "" and os.path.exists("ABYSS-P"):
                testIn = open("abyss-pe-orig", 'r')
                testOut = open("abyss-pe", 'w')
                for line in testIn.xreadlines():
                   if "which mpirun" in line:
                      testOut.write("mpirun?=$(shell which %s)\n"%(command))
                   elif "ifdef np" in line:
                      testOut.write(line)
                      testOut.write("ifneq ($(mpirun),mpirun)\n")
                   elif "ABYSS-P" in line:
                      testOut.write(line)
                      testOut.write("else\n")
                      testOut.write("\tABYSS $(abyssopt) $(ABYSS_OPTIONS) -o $@ $(in) $(se)\n")
                      testOut.write("endif\n")
                   else:
                      testOut.write(line)
                testIn.close()
                testOut.close()
             else:
                print "Error: cannot find MPI in your path. Disabling ABySS threading."
                os.system("cat abyss-pe-orig |awk -v found=0 -v skipping=0 '{if (match($0, \"ifdef np\")) {skipping=1; } if (skipping && match($1, \"ABYSS\")) {print $0; skipping=1; found=1} if (found && match($1, \"endif\")) {skipping=0;found = 0;} else if (skipping == 0) { print $0; } }' > abyss-pe")
             os.chdir("%s"%(METAMOS_ROOT))
             os.system("rm -rf sparsehash-2.0.2")
             os.system("rm -rf sparse.tar.gz")
             os.system("rm -rf abyss-1.3.6")
             os.system("rm -rf abyss.tar.gz")
             os.system("rm -rf boost_1_54_0")
             os.system("rm -rf boost.tar.gz")

    if not os.path.exists("./Utilities/cpp%s%s-%s%ssga"%(os.sep, OSTYPE, MACHINETYPE, os.sep)):
       sga = utils.getFromPath("sga", "SGA", False)
       if sga == "":
          if "sga" in packagesToInstall:
             dl = 'y'
          else:
             print "SGA binaries not found, optional for Assemble step, download now?"
             dl = raw_input("Enter Y/N: ")
          if dl == 'y' or dl == 'Y':
             os.system("curl -L https://sparsehash.googlecode.com/files/sparsehash-2.0.2.tar.gz -o sparse.tar.gz")
             os.system("tar xvzf sparse.tar.gz")
             os.chdir("sparsehash-2.0.2")
             os.system("./configure --prefix=`pwd`")
             updateMakeFileForDarwin("Makefile", addedCFlags, addedLDFlags)
             os.system("make install")
             os.chdir("%s"%(METAMOS_ROOT))
             os.system("curl -L https://github.com/pezmaster31/bamtools/archive/v2.3.0.tar.gz -o bamtools.tar.gz")
             os.system("tar xvzf bamtools.tar.gz")
             os.system("curl -L http://sourceforge.net/projects/bio-bwa/files/bwa-0.7.5a.tar.bz2 -o bwa.tar.bz2")
             os.system("tar xvjf bwa.tar.bz2")
             os.chdir("bwa-0.7.5a")
             os.system("make")
             os.chdir("%s"%(METAMOS_ROOT))
             os.system("curl -L https://github.com/jts/sga/archive/v0.10.10.tar.gz -o sga.tar.gz")
             os.system("tar xvzf sga.tar.gz")
             os.system("mv sga-0.10.10 ./Utilities/cpp%s%s-%s%ssga"%(os.sep, OSTYPE, MACHINETYPE, os.sep))
             os.system("mv bamtools-2.3.0 ./Utilities/cpp%s%s-%s%ssga/bamtools"%(os.sep, OSTYPE, MACHINETYPE, os.sep))
             os.system("mv sparsehash-2.0.2 ./Utilities/cpp%s%s-%s%ssga/sparsehash"%(os.sep, OSTYPE, MACHINETYPE, os.sep))
             os.chdir("./Utilities/cpp%s%s-%s%ssga/bamtools"%(os.sep, OSTYPE, MACHINETYPE, os.sep))
             os.system("mkdir build")
             os.chdir("build")
             os.system("export CC=`which gcc` && cmake ..")
             os.system("make")
             os.chdir("%s"%(METAMOS_ROOT))
             os.chdir("./Utilities/cpp%s%s-%s%ssga/src"%(os.sep, OSTYPE, MACHINETYPE, os.sep))
             # sparse hash library has unused variables which cause warnings with gcc 4.8 so disable -Werror
             if GCC_VERSION >= 4.8:
                os.system("mv configure.ac configure.original")
                os.system("cat configure.original |sed s/\-Werror//g > configure.ac")
             os.system("sh ./autogen.sh")
             os.system("./configure --with-sparsehash=`pwd`/../sparsehash --with-bamtools=`pwd`/../bamtools --prefix=`pwd`/../")
             updateMakeFileForDarwin("Makefile", addedCFlags, addedLDFlags)
             os.system("make install")
             os.chdir("%s"%(METAMOS_ROOT))
             os.system("mv bwa-0.7.5a/bwa ./Utilities/cpp%s%s-%s%ssga/bin/"%(os.sep, OSTYPE, MACHINETYPE, os.sep))
             os.system("cp %s/Utilities/cpp%s%s-%s%ssamtools %s/Utilities/cpp%s%s-%s%ssga/bin%ssamtools"%(METAMOS_ROOT, os.sep, OSTYPE, MACHINETYPE, os.sep, METAMOS_ROOT, os.sep, OSTYPE, MACHINETYPE, os.sep, os.sep))
             os.system("rm -rf sparsehash-2.0.2")
             os.system("rm -rf sparse.tar.gz")
             os.system("rm -rf bamtools-2.3.0")
             os.system("rm -rf bamtools.tar.gz")
             os.system("rm -rf sga-0.10.10")
             os.system("rm -rf sga.tar.gz")
             os.system("rm -rf bwa.tar.bz2")
             os.system("rm -rf bwa-0.7.5a")

    if not os.path.exists("./Utilities/cpp%s%s-%s%sedena"%(os.sep, OSTYPE, MACHINETYPE, os.sep)):
       edena = utils.getFromPath("edena", "EDENA", False)
       if "edena" in packagesToInstall:
          dl = 'y'
       else:
          print "Edena binaries not found, optional for Assemble step, download now?"
          dl = raw_input("Enter Y/N: ")
       if dl == 'y' or dl == 'Y':
          os.system("curl -L ftp://ftp.cbcb.umd.edu/pub/data/metamos/EdenaV3_130110.tar.gz -o edena.tar.gz")
          os.system("tar xvzf edena.tar.gz")
          os.system("mv EdenaV3.130110 ./Utilities/cpp%s%s-%s%sedena"%(os.sep, OSTYPE, MACHINETYPE, os.sep))
          os.chdir("./Utilities/cpp%s%s-%s%sedena"%(os.sep, OSTYPE, MACHINETYPE, os.sep))
          updateMakeFileForDarwin("src/Makefile", addedCFlags, addedLDFlags)
          os.system("make")
          os.chdir("%s"%(METAMOS_ROOT))
          os.system("rm -rf edena.tar.gz")

    if not os.path.exists("./quast"):
        if "quast" in packagesToInstall:
           dl = 'y'
        else:
           print "QUAST tool not found, optional for Validate step, download now?"
           dl = raw_input("Enter Y/N: ")
        if dl == 'y' or dl == 'Y':
            os.system("curl -L http://downloads.sourceforge.net/project/quast/quast-2.2.tar.gz -o quast.tar.gz")
            os.system("tar xvzf quast.tar.gz")
            os.system("mv ./quast-2.2 ./quast")
            os.system("rm -rf quast.tar.gz")

            # since quast requires a reference, also download refseq
            ftpSite = "ftp://ftp.ncbi.nih.gov/genomes/"
            file = "all.fna.tar.gz"
            if not os.path.exists("./Utilities/DB/refseq/") and not nodbs:
                print "Downloading refseq genomes (Bacteria/%s, Viruses/%s)..."%(file,file)
                print "\tThis file is large and may take time to download"
                os.system("curl -L %s/Bacteria/%s -o bacteria.tar.gz"%(ftpSite, file))
                os.system("curl -L %s/Viruses/%s -o viruses.tar.gz"%(ftpSite, file))
                os.system("mkdir -p ./Utilities/DB/refseq/temp")
                os.system("mv bacteria.tar.gz ./Utilities/DB/refseq/temp")
                os.system("mv viruses.tar.gz  ./Utilities/DB/refseq/temp")
                os.chdir("./Utilities/DB/refseq/temp")
                os.system("tar xvzf bacteria.tar.gz")
                os.system("tar xvzf viruses.tar.gz")
                os.chdir("..")
                print "Current directory is %s"%(os.getcwd())
                for file in os.listdir("%s/temp"%(os.getcwd())):
                    file = "%s%stemp%s%s"%(os.getcwd(), os.sep, os.sep, file)
                    if os.path.isdir(file):
                        prefix = os.path.splitext(os.path.basename(file))[0]
                        os.system("cat %s/*.fna > %s.fna"%(file, prefix))
                os.system("rm -rf temp")
                os.chdir("%s"%(METAMOS_ROOT))
   
    if not os.path.exists("./Utilities/cpp%s%s-%s%sfreebayes"%(os.sep, OSTYPE, MACHINETYPE, os.sep)):
        if "freebayes" in packagesToInstall:
           dl = 'y'
        else:
           print "FreeBayes tool not found, optional for Validate step, download now?"
           dl = raw_input("Enter Y/N: ")
        if dl == 'y' or dl == 'Y':
           os.system("git clone --recursive git://github.com/ekg/freebayes.git freebayes")
           os.system("mv ./freebayes ./Utilities/cpp/%s%s-%s%sfreebayes"%(os.sep, OSTYPE, MACHINETYPE, os.sep))
           os.chdir("./Utilities/cpp/%s%s-%s%sfreebayes"%(os.sep, OSTYPE, MACHINETYPE, os.sep))
           updateMakeFileForDarwin("src/makefile", addedCFlags, addedLDFlags)

           # dont set static building libs on OSX, sseems to cause compile issues for jellyfish
           os.environ["CFLAGS"] = oldCFlags
           os.environ["CPPFLAGS"] = oldCPPFlags
           os.environ["CXXFLAGS"] = oldCXXFlags
           os.environ["LDFLAGS"] = oldLDFlags
           os.system("make")
           os.chdir("%s"%(METAMOS_ROOT))
           if OSTYPE == "Darwin":
              # reset env variables again
              addEnvironmentVar("CFLAGS", " %s "%(addedCFlags))
              addEnvironmentVar("CPPFLAGS", " %s "%(addedCFlags))
              addEnvironmentVar("CXXFLAGS", " %s "%(addedCFlags))
              addEnvironmentVar("LDFLAGS", " %s "%(addedLDFlags))
           os.system("make")
           os.chdir("%s"%(METAMOS_ROOT))

    if not os.path.exists("./Utilities/cpp%s%s-%s%scgal"%(os.sep, OSTYPE, MACHINETYPE, os.sep)):
        if "cgal" in packagesToInstall:
           dl = 'y'
        else:
           print "CGAL tool not found, optional for Validate step, download now?"
           dl = raw_input("Enter Y/N: ")
        if dl == 'y' or dl == 'Y':
            os.system("curl -L http://bio.math.berkeley.edu/cgal/cgal-0.9.6-beta.tar -o cgal.tar")
            os.system("tar xvf cgal.tar")
            os.system("mv cgal-0.9.6-beta ./Utilities/cpp/%s%s-%s%scgal"%(os.sep, OSTYPE, MACHINETYPE, os.sep))
            os.chdir("./Utilities/cpp/%s%s-%s%scgal"%(os.sep, OSTYPE, MACHINETYPE, os.sep))
            updateMakeFileForDarwin("makefile", addedCFlags, addedLDFlags, True)
            os.system("make")
            os.chdir("%s"%(METAMOS_ROOT))
            os.system("rm -rf cgal.tar")

    if not os.path.exists("./Utilities/cpp%s%s-%s%sREAPR"%(os.sep, OSTYPE, MACHINETYPE, os.sep)):
       if "reapr" in packagesToInstall:
          dl = 'y'
       else:
           print "REAPR tool not found, optional for Validate step, download now?"
           dl = raw_input("Enter Y/N: ")
       if dl == 'y' or dl == 'Y':
          os.system("curl -L ftp://cbcb.umd.edu/pub/data/metamos/Reapr_1.0.16.tar.gz -o reapr.tar.gz")
          os.system("tar xvzf reapr.tar.gz")
          os.system("mv Reapr_1.0.16 ./Utilities/cpp/%s%s-%s%sREAPR"%(os.sep, OSTYPE, MACHINETYPE, os.sep))

          # find cmake we installed anyway
          if not os.path.exists("./Utilities/cpp%s%s-%s%scmake"%(os.sep, OSTYPE, MACHINETYPE, os.sep)):
             cmake = utils.getFromPath("cmake", "CMAKE", False) + os.sep + "cmake"
          else:
             cmake="%s/Utilities/cpp%s%s-%s%scmake/bin/cmake"%(METAMOS_ROOT, os.sep, OSTYPE, MACHINETYPE, os.sep)

          filespec = utils.getCommandOutput("perl -MFile::Spec::Link -e 0 && echo $?", True)
          if filespec == "":
             os.system("curl -L http://search.cpan.org/CPAN/authors/id/R/RM/RMBARKER/File-Copy-Link-0.113.tar.gz -o file.tar.gz")
             os.system("tar xvzf file.tar.gz")
             os.chdir("File-Copy-Link-0.113")
             os.system("perl Makefile.PL PREFIX=`pwd`/build")
             os.system("make install")
             os.chdir("%s"%(METAMOS_ROOT))
             pathToCopy = utils.getCommandOutput("find File-Copy-Link-0.113/build -type d -name \"File\" |grep -v auto", False)
             copyPerlLib(pathToCopy, "./Utilities/cpp%s%s-%s%sREAPR/lib/"%(os.sep, OSTYPE, MACHINETYPE, os.sep))
             os.system("rm -rf file.tar.gz")
             os.system("rm -rf File-Copy-Link-0.113")
             libUpdate = "%s/Utilities/cpp%s%s-%s%sREAPR/lib/"%(METAMOS_ROOT, os.sep, OSTYPE, MACHINETYPE, os.sep)
             if "PERL5LIB" in os.environ:
                libUpdate = "%s%s%s"%(os.environ["PERL5LIB"], os.pathsep, libUpdate)
             os.environ["PERL5LIB"]=libUpdate

          if OSTYPE == "Darwin":
             os.chdir("./Utilities/cpp/%s%s-%s%sREAPR/third_party/snpomatic/src"%(os.sep, OSTYPE, MACHINETYPE, os.sep))
             os.system("cp snpomatic.h snpomatic.original") 
             os.system("cat snpomatic.original |awk '{if (match($0, \"#include <algorithm>\")) { print $0; print \"#define ulong u_long\"; } else { print $0} }' > snpomatic.h")
             os.chdir("%s"%(METAMOS_ROOT))

             # also need smalt, the reapr distro comes with linux 64 bit only
             os.chdir("./Utilities/cpp/%s%s-%s%sREAPR/third_party"%(os.sep, OSTYPE, MACHINETYPE, os.sep))
             updateMakeFileForDarwin("tabix/Makefile", addedCFlags, addedLDFlags)
             updateMakeFileForDarwin("snpomatic/Makefile", addedCFlags, addedLDFlags)
             os.system("curl -L http://sourceforge.net/projects/smalt/files/smalt-0.7.5.tar.gz -o smalt.tar.gz")
             os.system("tar xvzf smalt.tar.gz")
             os.chdir("./smalt-0.7.5")
             os.system("./configure --prefix=`pwd`/build")
             updateMakeFileForDarwin("Makefile", addedCFlags, addedLDFlags)
             os.system("make install")

             os.chdir("..")
             os.system("rm smalt_x86_64")
             os.system("rm -rf smalt.tar.gz")
             os.system("ln -s smalt-0.7.5/build/bin/smalt smalt_x86_64")
             os.chdir("%s"%(METAMOS_ROOT))
              
          os.chdir("./Utilities/cpp/%s%s-%s%sREAPR"%(os.sep, OSTYPE, MACHINETYPE, os.sep))
          # samtools which reapr includes uses curses lib which is optional so disable it if not found
          os.system("echo \"#include <curses.h>\" > .test.h")
          HAVE_CURSES=utils.getCommandOutput("gcc .test.h && echo $?", True)
          if HAVE_CURSES == "":
             os.chdir("third_party/samtools")
             os.system("mv Makefile Makefile.original")
             os.system("cat Makefile.original | sed s/\-lcurses//g |sed s/\-D_CURSES_LIB=1//g > Makefile")
             os.chdir("../../")

          # reapr comes with its own cmake which has issues building on recent gcc
          # kill it it and use our own
          os.system("cp install.sh install.sh.orig")
          testIn = open("install.sh.orig", 'r')
          testOut = open("install.sh", 'w')
          isSkip = 0;
          for line in testIn.xreadlines():
             if "cmake/bin/cmake" in line:
                testOut.write("%s ..\n"%(cmake))
             elif "cd cmake" in line:
                # skip some lines
                isSkip = 3 
             elif isSkip > 0:
                isSkip -= 1
             else:
                testOut.write(line.strip() + "\n")
          testIn.close()
          testOut.close()

          os.system("export CC=`which gcc` && sh install.sh force")
          os.system("chmod ug+x third_party/smalt_x86_64")
          os.chdir("%s"%(METAMOS_ROOT))

          if os.path.exists("./Utilities/cpp%s%s-%s%sREAPR/lib"%(os.sep, OSTYPE, MACHINETYPE, os.sep)):
             os.chdir("./Utilities/cpp%s%s-%s%sREAPR/"%(os.sep, OSTYPE, MACHINETYPE, os.sep))
             os.system("cp reapr reapr.original")
             os.system("cat reapr.original |awk '{if (match($0, \"use strict\")) { print \"use lib \\\"%s/Utilities/cpp%s%s-%s%sREAPR/lib\\\";\"; print $0; } else { print $0}}' > reapr"%(METAMOS_ROOT, os.sep, OSTYPE, MACHINETYPE, os.sep))
             os.chdir("%s"%(METAMOS_ROOT))

          os.chdir("./Utilities/cpp%s%s-%s%sREAPR/src"%(os.sep, OSTYPE, MACHINETYPE, os.sep))
          # fix samtools link
          os.system("rm samtools")
          os.system("cp ../third_party/samtools/samtools ./") 

          # REAPR has a bug where fasta headers with commas are not properly fixed, patch the bug
          os.system("cp task_facheck.pl task_facheck.pl.original")
          os.system("cat task_facheck.pl.original |awk -v quote=\"'\"  '{if (match($0, \"new_id =~\")) { print \"    $new_id =~ s/[;\"quote\"|:,\\\\+\\\\-\\\\s\\\\(\\\\)\\\\{\\\\}\\\\[\\\\]]/_/g;\"; } else { print $0}}' > task_facheck.pl")
          os.chdir("%s"%(METAMOS_ROOT))
          os.system("rm -rf reapr.tar.gz")
    
    if not os.path.exists("./Utilities/cpp%s%s-%s%sFRCbam"%(os.sep, OSTYPE, MACHINETYPE, os.sep)):
        if "frcbam" in packagesToInstall:
           dl = 'y'
        else:
           print "FRCbam tool not found, optional for Validate step, download now?"
           dl = raw_input("Enter Y/N: ")
        if dl == 'y' or dl == 'Y':
            os.system("curl -L ftp://cbcb.umd.edu/pub/data/metamos/FRC_align-master.zip -o frcbam.zip")
            os.system("unzip frcbam.zip")
            os.system("mv FRC_align-3398ca469b2077d6672b85317eee6fea171b6a27 ./Utilities/cpp/%s%s-%s%sFRCbam"%(os.sep, OSTYPE, MACHINETYPE, os.sep))
            os.chdir("./Utilities/cpp/%s%s-%s%sFRCbam/src/samtools"%(os.sep, OSTYPE, MACHINETYPE, os.sep))
            # samtools which frcbam includes uses curses lib which is optional so disable it if not found
            os.system("echo \"#include <curses.h>\" > .test.h")
            HAVE_CURSES=utils.getCommandOutput("gcc .test.h && echo $?", True)
            if HAVE_CURSES == "":
               os.system("mv Makefile Makefile.original")
               os.system("cat Makefile.original | sed s/\-lcurses//g |sed s/\-D_CURSES_LIB=1//g > Makefile")
            updateMakeFileForDarwin("Makefile", addedCFlags, addedLDFlags)
            os.system("make")
            os.chdir("%s/Utilities/cpp/%s%s-%s%sFRCbam"%(METAMOS_ROOT, os.sep, OSTYPE, MACHINETYPE, os.sep))
            boostFlags = ""
            if os.path.exists("/opt/local/lib/libboost_system-mt.a"):
               os.environ["LDFLAGS"]="-L/opt/local/lib -lboost_system-mt"
            elif os.path.exists("/opt/local/lib/libboost_system.a"):
               os.environ["LDFLAGS"]="-L/opt/local/lib -lboost_system"
            elif os.path.exists("/usr/lib/libboost_system-mt.a"):
               os.environ["LDFLAGS"]="-L/usr/lib -lboost_system-mt"
            elif os.path.exists("/usr/lib/libboost_system.a"):
               os.environ["LDFLAGS"]="-L/usr/lib -lboost_system"
            else:
               # install boost ourselves
               os.system("curl -L http://sourceforge.net/projects/boost/files/boost/1.54.0/boost_1_54_0.tar.gz -o boost.tar.gz")
               os.system("tar xvzf boost.tar.gz")
               os.chdir("boost_1_54_0")
               os.system("sh bootstrap.sh")
               os.system("./b2 install --prefix=`pwd`/build threading=multi")
               ldflags = "-L%s/build/lib -lboost_system"%(os.getcwd())
               if os.path.exists("%s/build/lib/libboost_system-mt.a"%(os.getcwd())):
                  ldflags = "-L%s/build/lib -lboost_system-mt"%(os.getcwd())
               os.environ["LDFLAGS"]=ldflags
               try:
                   os.environ["LD_LIBRARY_PATH"] = os.environ["LD_LIBRARY_PATH"] + os.pathsep + "%s/build/lib"%(os.getcwd())
               except KeyError:
                   os.environ["LD_LIBRARY_PATH"] = "%s/build/lib"%(os.getcwd())
               boostFlags = "--with-boost=%s/build/ --disable-shared --enable-static-boost --enable-static-FRC"%(os.getcwd())
               os.chdir("..")
               os.system("rm -rf boost.tar.gz")
    
            os.system("./configure --prefix=%s/Utilities/cpp/%s%s-%s%sFRCbam/ %s"%(METAMOS_ROOT, os.sep, OSTYPE, MACHINETYPE, os.sep, boostFlags))
            updateMakeFileForDarwin("Makefile", addedCFlags, addedLDFlags)
            os.system("make install")
            if boostFlags != "":
               os.system("cp boost_1_54_0/build/lib/* ./bin")

            os.chdir("%s"%(METAMOS_ROOT))
            os.system("rm -rf frcbam.zip")

    if not os.path.exists("./Utilities/cpp/%s%s-%s%sALE"%(os.sep, OSTYPE, MACHINETYPE, os.sep)):
        if "ale" in packagesToInstall:
           dl = 'y'
        else:
           print "ALE tool not found, optional for Validate step, download now?"
           dl = raw_input("Enter Y/N: ")
        if dl == 'y' or dl == 'Y':
           os.system("curl -L ftp://ftp.cbcb.umd.edu/pub/data/metamos/ale.tar.gz -o ale.tar.gz")
           os.system("tar xvzf ale.tar.gz")
           os.system("mv ALE ./Utilities/cpp/%s%s-%s%sALE"%(os.sep, OSTYPE, MACHINETYPE, os.sep))
           os.chdir("./Utilities/cpp/%s%s-%s%sALE/src"%(os.sep, OSTYPE, MACHINETYPE, os.sep))
           updateMakeFileForDarwin("makefile", addedCFlags, addedLDFlags)
           os.system("make all") 
           os.chdir("%s"%(METAMOS_ROOT))
           os.system("rm -rf ale.tar.gz")

if "deprecated" in enabledWorkflows or manual:
    if not os.path.exists("./Utilities/glimmer-mg"):
        if "glimmer-mg" in packagesToInstall:
           dl = 'y'
        else:
           print "Glimmer-MG not found, optional for FindORFS step. Caution, this will take approx. 24 hours to complete, including Phymm download & install. download & install now?"
           dl = raw_input("Enter Y/N: ")
        if (dl == 'y' or dl == 'Y') or not nodbs:
            archive = "glimmer-mg-0.3.1.tar.gz"
            os.system("curl -L ftp://ftp.cbcb.umd.edu/pub/data/metamos/%s -o %s" %(archive, archive))
            os.system("tar -C ./Utilities/ -xvf %s" % archive)
            os.system("rm %s"%archive)
            os.system("python ./Utilities/glimmer-mg/install_glimmer.py")

# should check for success of installation
workflow.updateSupportedWorkflows(enabledWorkflows)

os.environ["CFLAGS"] = oldCFlags
os.environ["CPPFLAGS"] = oldCPPFlags
os.environ["CXXFLAGS"] = oldCXXFlags
os.environ["LDFLAGS"] = oldLDFlags

sys.path.append(METAMOS_ROOT + os.sep + "Utilities" + os.sep + "python")
from get_setuptools import use_setuptools
use_setuptools()

print "Run setup.py.."
os.system("python setup.py install_scripts --install-dir=`pwd` build_ext")
#print "Compile & optimize"
#distutils.util.byte_compile(['./runPipeline.py'],optimize=2,force=True)
#os.system("chmod a+wrx runPipeline.pyo")
os.system("mv runPipeline.py runPipeline")
os.system("mv initPipeline.py initPipeline")

#remove imports from pth file, if exists                                                                                                                                                          
nf = []
try:
    dir1 = utils.INITIAL_UTILS+os.sep+"python"+os.sep+"lib"+os.sep+"python"
    if not os.path.exists(dir1+os.sep+"easy-install.pth"):
        dir1 = utils.INITIAL_UTILS+os.sep+"python"+os.sep+"lib64"+os.sep+"python"

    nf = open(dir1+os.sep+"easy-install.pth",'r')
    ndata = []
    for line in nf.xreadlines():
        if "import" in line:
            continue
        ndata.append(line)
    nf.close()
    nfo = open(dir1+os.sep+"easy-install.pth",'w')
    for line in ndata:
        nfo.write(line)
    nfo.close()
except IOError:
    pass

validate_install = 0
if validate_install:
    import check_install
    rt = check_install.validate_dir(METAMOS_ROOT,'required_file_list.txt')
    if rt == -1:
        print "MetAMOS not properly installed, please reinstall or contact development team for assistance"
        sys.exit(1)
    
