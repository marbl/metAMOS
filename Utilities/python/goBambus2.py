#pipeline script for assembly + Bambus 2
#contributed by Todd J Treangen

import string, sys, os, subprocess#, spincursor

RED =    "\033[0;31m"
GREEN =  "\033[0;32m"
BLUE =   "\033[0;34m"
NONE =   "\033[0m"

#set a path to AMOS, to local directory, or if already installed in PATH , set to ""
AMOSDIR = ""
LOGFILE = "bambus2.log"

if __name__ == "__main__":
    logfile = open(LOGFILE,'w')
    #input can be CONTIGS, or READS that need assembly into CONTIGS, plus need LINK INFO
    #IF CONTIGS, no need for minimus, but need to call tarchive2amos -c contigfile
    #IF NO CONTIGS, use minimus as default (other also ?)
    #NEED TO ADD IN OTHER USEFUL PARAMETERS, ADD IN PYTHON FILE TO BUILD SYSTEM

    #INPUT: CONTIGS + MATES, READS + TRACE XML, OR AMOS BANK (3 options)
    usage = "\nrun: goBambus2.py <input reads or contigs or amos bank name> <output prefix> [options]\n"
    usage += "eg.: goBambus2.py example.contigs myoutput --all --contigs\n"
    usage += "This script is designed to run the Bambus pipeline and takes either reads or contigs plus XML Trace Archive data as input and outputs scaffolds\n"
    usage += "For further info please contact the Bambus 2 authors: Sergey Koren and Mihai Pop\n"
    if len(sys.argv) < 3:
        print usage
        sys.exit(1)
    if len(sys.argv) > 3:
       options = sys.argv[2:]
       #print "options: ", options
    opts = ["2amos","minimus","clk","bundle","reps","orient","untangle","2fasta","printscaff"]
    xopts = ["new","verbose","all","contigs","reads"]
    #-2amos
    #-minimus
    #-clk
    #-bundle
    #-reps
    #-orient
    #-untangle
    #-2fasta
    #-printScaff
    opt_dict = {}
    xopt_dict = {}
    for opt in opts:
        opt_dict[opt] = [""]

    for opt in xopts:
        xopt_dict[opt] = 0
        
    if len(sys.argv) > 3:
        for opt in opt_dict.keys():
            opt_dict[opt] = []
        for opt in sys.argv[3:]:
            sopt = opt.replace("-","")
            if sopt in opt_dict.keys():
                opt_list = sopt.split(",")
                mopt = opt_list[0]
                xopt = opt_list[1:]
            
                if len(xopt) == 0:
                    xopt = [" "]
                opt_dict[mopt] = xopt
            else:
                try:
                    xopt_dict[sopt] = 1
                except KeyError:
                    pass
            
    prefix =  sys.argv[1].split(".")[0]
    outprefix = sys.argv[2]
    sequence_file = os.getcwd() + "/" + sys.argv[1]

    vtext = subprocess.PIPE
    if xopt_dict["verbose"] == 1:
        vtext = None
    #see if trace archive file present
    xml_file = os.getcwd() + "/" + prefix +".xml"
    exists_xml = True
    if not os.path.exists(xml_file):
        exists_xml = False

    #need to run
    amosbank = outprefix+".bnk"
    if xopt_dict["contigs"] == 1:
        mate_file = os.getcwd() + "/" + prefix +".mates"
        contig_file = os.getcwd() + "/" + prefix +".contig"
        afg_file = os.getcwd() + "/" + outprefix +".afg"
        bank_file = os.getcwd() + "/" + outprefix +".bnk"
        if not os.path.exists(mate_file):
           print "mate pair info not present!"
           sys.exit(1)
        if not os.path.exists(contig_file):
           print "contig info not present!"
           sys.exit(1)
        if not os.path.exists(afg_file) or 1:
            p = subprocess.Popen(AMOSDIR+"toAmos -o %s.afg -s %s -m %s -c %s"%(outprefix,sys.argv[1], mate_file, contig_file), shell=True, stdin=subprocess.PIPE, stdout=vtext, stderr=logfile)

            if xopt_dict["verbose"] == 1:
                print "1) running toAmos"
            
            else:
                print "1) running toAmos...",
            sys.stdout.flush()
            sts = os.waitpid(p.pid,0)[1]

            if xopt_dict["verbose"] != 1:
                 if sts == 0:
                     print "\t\t%s...success%s"%(GREEN,NONE)
                 else:
                     print "\t\t%s...failed%s"%(RED,NONE)
                     sys.exit(1)

        else:
            print "1) running toAmos: AMOS format AFG file already existing, will use existing file"


        if not os.path.exists(bank_file) or 1:
            p = subprocess.Popen(AMOSDIR+"bank-transact -cf -b %s -m %s.afg"%(amosbank,outprefix), shell=True, stdin=subprocess.PIPE, stdout=vtext, stderr=logfile)
            if xopt_dict["verbose"] == 1:
                print "2) running bank-transact"
            
            else:
                print "2) running bank-transact...",
            sys.stdout.flush()
            sts = os.waitpid(p.pid,0)[1]

            if xopt_dict["verbose"] != 1:
                if sts == 0:
                    print "\t\t%s...success%s"%(GREEN,NONE)
                else:
                    print "\t\t%s...failed%s"%(RED,NONE)
                    sys.exit(1)
        else:
            print "2) running bank-transact: BANK file already exists, skipping"

        
    elif xopt_dict["reads"] == 1 and (not exists_xml or xopt_dict["2amos"] ==1):
            mate_file = os.getcwd() + "/" + prefix +".mates"
            afg_file = os.getcwd() + "/" + outprefix +".afg"
            bank_file = os.getcwd() + "/" + outprefix +".bnk"
            if not os.path.exists(afg_file) and not os.path.exists(mate_file):
                print "mate pair info not present!"
                sys.exit(1)

            if not os.path.exists(afg_file) or 0:
                p = subprocess.Popen(AMOSDIR+"toAmos -o %s.afg -s %s -m %s"%(outprefix,sys.argv[1], mate_file), shell=True, stdin=subprocess.PIPE, stdout=vtext, stderr=logfile)
            
                if xopt_dict["verbose"] == 1:
                    print "1) running toAmos"
            
                else:
                    print "1) running toAmos...",

                sys.stdout.flush()
                sts = os.waitpid(p.pid,0)[1]

                if xopt_dict["verbose"] != 1:
                    if sts == 0:
                        print "\t\t%s...success%s"%(GREEN,NONE)
                    else:
                        print "\t\t%s...failed%s"%(RED,NONE)
                        sys.exit(1)

            else:
                print "1) running toAmos: AMOS format AFG file already existing, will use existing file"


            p = subprocess.Popen(AMOSDIR+"minimus %s"%(outprefix), shell=True, stdin=subprocess.PIPE, stdout=vtext, stderr=logfile)

            if xopt_dict["verbose"] == 1:
                print "2) running minimus"
            else:
                print "2) running minimus...",
            sys.stdout.flush()
            sts = os.waitpid(p.pid,0)[1]
            if xopt_dict["verbose"] != 1:
                if sts == 0:
                    print "\t\t%s...success%s"%(GREEN,NONE)
                else:
                    print "\t\t%s...failed%s"%(RED,NONE)
                    sys.exit(1)

        
    elif xopt_dict["reads"] == 1:
        p = subprocess.Popen(AMOSDIR+"tarchive2amos -o %s %s"%(outprefix,sys.argv[1]), shell=True, stdin=subprocess.PIPE, stdout=vtext, stderr=logfile)
        if xopt_dict["verbose"] == 1:
            print "1) running tarchive2amos"
            
        else:
            print "1) running tarchive2amos...",
        sys.stdout.flush()
        sts = os.waitpid(p.pid,0)[1]
        if xopt_dict["verbose"] != 1:
            if sts == 0:
                print "\t\t%s...success%s"%(GREEN,NONE)
            else:
                print "\t\t%s...failed%s"%(RED,NONE)
                sys.exit(1)

        p = subprocess.Popen(AMOSDIR+"minimus %s"%(outprefix), shell=True, stdin=subprocess.PIPE, stdout=vtext, stderr=logfile)

        if xopt_dict["verbose"] == 1:
            print "2) running minimus"
        else:
            print "2) running minimus...",
        sys.stdout.flush()
        sts = os.waitpid(p.pid,0)[1]
        if xopt_dict["verbose"] != 1:
            if sts == 0:
                print "\t\t%s...success%s"%(GREEN,NONE)
            else:
                print "\t\t%s...failed%s"%(RED,NONE)
                sys.exit(1)


    else:
        #assuming existing good AMOS bank
        amosbank = sys.argv[1]
        amosbank_file = os.getcwd() + "/" + amosbank
    
        if not os.path.exists(amosbank_file):
           print "AMOS bank does not exist!"
           sys.exit(1)
    
    if xopt_dict["all"] == 1 or len(opt_dict["clk"]) > 0:
        p = subprocess.Popen(AMOSDIR+"clk -b %s"%(amosbank), shell=True, stdin=subprocess.PIPE, stdout=vtext, stderr=logfile)
        if xopt_dict["verbose"] == 1:
            print "3) running clk"
        else:
            print "3) running clk...",
        sys.stdout.flush()
        sts = os.waitpid(p.pid,0)[1]
        if xopt_dict["verbose"] != 1:
            if sts == 0:
                print "\t\t%s...success%s"%(GREEN,NONE)
            else:
                print "\t\t%s...failed%s"%(RED,NONE)
                sys.exit(1)

    if xopt_dict["all"] == 1 or len(opt_dict["bundle"]) > 0:
        p = subprocess.Popen(AMOSDIR+"Bundler -b %s"%(amosbank), shell=True, stdin=subprocess.PIPE, stdout=vtext, stderr=logfile)
        if xopt_dict["verbose"] == 1:
            print "4) running Bundler"
        else:
            print "4) running Bundler...",
        sys.stdout.flush()
        sts = os.waitpid(p.pid,0)[1]
        if xopt_dict["verbose"] != 1:
            if sts == 0:
                print "\t\t%s...success%s"%(GREEN,NONE)
            else:
                print "\t\t%s...failed%s"%(RED,NONE)
                sys.exit(1)

    if xopt_dict["all"] == 1 or len(opt_dict["reps"]) > 0:
        repfile = open("myreps",'w')
        p = subprocess.Popen(AMOSDIR+"MarkRepeats -b %s"%(amosbank), shell=True, stdin=subprocess.PIPE, stdout=repfile, stderr=logfile)
        if xopt_dict["verbose"] == 1:
            print "5) running MarkRepeats"
        else:
            print "5) running MarkRepeats...",
        sys.stdout.flush()
        sts = os.waitpid(p.pid,0)[1]
        if xopt_dict["verbose"] != 1:
            if sts == 0:
                print "\t\t%s...success%s"%(GREEN,NONE)
            else:
                print "\t\t%s...failed%s"%(RED,NONE)
                sys.exit(1)

        #quality control, check how many repeats file has
        repfile.close()
        repfile = open("myreps",'r')
        nreps = len(repfile.readlines())
        print "\t%s %s repeats found%s"%(BLUE,nreps,NONE)
        
    if xopt_dict["all"] == 1 or len(opt_dict["orient"]) > 0:
        p = subprocess.Popen(AMOSDIR+"OrientContigs -b %s.bnk -repeats myreps -prefix %s  -noreduce -linearize"%(outprefix, prefix+".scaff"), shell=True, stdin=subprocess.PIPE, stdout=vtext, stderr=logfile)

        if xopt_dict["verbose"] == 1:
            print "6) running OrientContigs"
        else:
            print "6) running OrientContigs...",
        sys.stdout.flush()
        sts = os.waitpid(p.pid,0)[1]
        if xopt_dict["verbose"] != 1:
            if sts == 0:
                print "\t\t%s...success%s"%(GREEN,NONE)
            else:
                print "\t\t%s...failed%s"%(RED,NONE)
                sys.exit(1)


    if xopt_dict["all"] == 1 or len(opt_dict["2fasta"]) > 0:
        contigfile = open("%s.contigs.fasta"%(outprefix),'w')
        p = subprocess.Popen(AMOSDIR+"bank2fasta -d -b %s"%(amosbank), shell=True, stdin=subprocess.PIPE, stdout=contigfile, stderr=logfile)
        if xopt_dict["verbose"] == 1:
            print "7) running bank2fasta"
        else:
            print "7) running bank2fasta...",
        sys.stdout.flush()
        sts = os.waitpid(p.pid,0)[1]
        if xopt_dict["verbose"] != 1:
            if sts == 0:
                print "\t\t%s...success%s"%(GREEN,NONE)
            else:
                print "\t\t%s...failed%s"%(RED,NONE)
                sys.exit(1)

        contigfile.close()
        
        contigfile = open("%s.contigs.fasta"%(outprefix),'r')
        contigdata = contigfile.read()
        ncontigs = contigdata.count(">")
        print "\t%s %s contigs%s"%(BLUE,ncontigs,NONE)
    if xopt_dict["all"] == 1 or len(opt_dict["printscaff"]) > 0:
        p = subprocess.Popen(AMOSDIR+"printScaff -e %s -s %s -l %s -f %s -merge -o %s"%(prefix+".scaff.evidence.xml",prefix+".scaff.out.xml",prefix+".scaff.library",outprefix+".contigs.fasta",outprefix+".scaffold"), shell=True, stdin=subprocess.PIPE, stdout=vtext, stderr=logfile)
        if xopt_dict["verbose"] == 1:
            print "8) running printScaff"
        else:
            print "8) running printScaff...",
        sys.stdout.flush()
        sts = os.waitpid(p.pid,0)[1]
        if xopt_dict["verbose"] != 1:
            if sts == 0:
                print "\t\t%s...success!%s"%(GREEN,NONE)
            else:
                print "\t\t%s...failed%s"%(RED,NONE)
                sys.exit(1)


