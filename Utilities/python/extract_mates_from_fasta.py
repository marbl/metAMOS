import string, sys

if __name__ == "__main__":

    infile = sys.argv[1]
    #outprefix = sys.argv[2]
    #library_min = sys.argv[3]
    #library_max = sys.argv[4]
    f1 = open(infile,'r')
    #outprefix = outprefix+infile
    f2 = open("%s.mates"%(infile),'w')
    #f2.write("library\t110110\t%s\t%s\n"%(library_min,library_max))
    #lines = f1.readlines()
    first = 1
    second = 0
    firstmate = ""
    linecnt = 0
    for line in f1.xreadlines():
        #if linecnt % 2 == 0:#">" not in line:
        if ">" in line:

            line = line.replace(">","")
            line = line.replace("\n","")
            data = line.split(" ")
            mate= data[0]
            mate = mate.strip()
        else:
            linecnt +=1
            continue
        if first:
            firstmate = mate
            first = 0        
            second = 1
        elif second:
            f2.write(firstmate+"\t"+mate+"\n")
            f2.flush()
            first = 1
            second = 0
        else:
            linecnt+=1
            continue
        linecnt +=1
    f1.close()
    f2.close()
