import sys, string, os

if __name__ == "__main__":
    infile = sys.argv[1]
    print infile
    f = open(infile,'r')
    fdata = f.readlines()
    fout = open(infile[:-3]+"fastq",'w')
    
    #fdata = fout.readlines()
    for line in fdata:
        if ">" in line:
            fout.write(line)
            continue
        else:
            seq = ""
            toks = line.split(" ")
            if len(toks) < 2:
                fout.write(line)
                continue
            #print toks
            for char in toks:
                if char == "\n":
                    seq += char
                    break
                elif char == "" or char == " ":
                    continue
                else:
                    seq += chr(int(char)+62)
            if "\n" not in seq:
                seq += "\n"
            fout.write(seq)
    f.close()
    fout.close()
                
