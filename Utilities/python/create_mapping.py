import os, sys, string

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print "usage: create_mapping.py <translation file> <phymm output> <output file>"
#    r = open("class_key.tab",'r')
    r = open(sys.argv[1],'r')
#    p = open("phymm.out",'r')
    p = open(sys.argv[2],'r')
#    out = open("s12.annots",'w')
    out = open(sys.argv[3],'w')
    class_dict = {}
    
    for line in r:
        line = line.replace("\n","")
        id,label = line.split("\t")
        class_dict[label] = id

    p.readline()
    out.write("contigID\tclassID\n")
    contig_dict = {}
    for line in p:
        id = ""
        data = line.split("\t")
        contig_id = data[0]
        phylum = data[11].split("-")[0]
        try:
            id = class_dict[phylum]
        except KeyError:
            print "phylum %s not found!"%phylum
            continue
        contig_dict[int(contig_id.split("_")[1])] = id
    keys = contig_dict.keys()
    keys.sort()
    for key in keys:
        out.write("contig_%d\t%s\n"%(key,contig_dict[key]))
    p.close()
    r.close()
    out.close()
    
       
