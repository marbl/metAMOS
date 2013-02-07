import os,string, sys
allfiles = []
MA_dir = ""

def myvisit(a, dir, files):
    global allfiles
    global MA_dir
    #print dir,": %d files"%len(files)                                                                                                                    
    for file in files:
        if not os.path.isdir(file):
            zf = dir+os.sep+file
            if "." not in zf.rsplit(os.sep,1)[-1]:
                continue
            if "~" in zf:
                continue
            if ".git" in zf:
                continue
            #print MA_dir
            #print zf
            zf = zf.replace(MA_dir,"")
            #print zf
            if zf[0] == "/":
                zf = zf[1:]
            if zf not in allfiles:
                allfiles.append(zf)

def validate_dir(MAdir,file_list):
    global allfiles
    global MA_dir
    MA_dir = MAdir
    file_dict = {}
    ff = open(file_list,'r')
    for line in ff.xreadlines():
        line = line.replace("\n","")
        file_dict[line] = 1

    os.path.walk(MA_dir,myvisit,None)
    cf = open('allfiles.txt','w')
    for file in allfiles:
        cf.write(file+"\n")
    cf.close()
    missing_files = []
    for key in file_dict.keys():
        #key = "/" + key
        key = key.replace("\n","")
        if key.strip() not in allfiles:
            print "ERROR: Required file %s not found! Please reinstall MetAMOS"%(MA_dir+os.sep+key)
            sys.exit(1)

    print "Installation in %s looks ok, have fun!"%(MA_dir)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print "usage: python check_install.py <MetAMOS install dir> <file list>"
        sys.exit(1)

    MA_dir = sys.argv[1]
    file_list = sys.argv[2]

    if not os.path.exists(MA_dir):
        print MA_dir, " does not exist!"
        sys.exit(1)

    if not os.path.exists(file_list):
        print file_list, "does not exist!"
        sys.exit(1)

    validate_dir(MA_dir,file_list)

