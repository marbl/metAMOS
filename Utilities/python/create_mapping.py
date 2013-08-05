import os, sys, string

#if __name__ == "__main__":
def create_mapping (first,second,third):
    #if len(sys.argv) < 3:
    #    print "usage: create_mapping.py <translation file> <metaphyler output> <output file>"

    r = open(first,'r')
    p = open(second,'r')
    out = open(third,'w')
    class_dict = {}
    
    for line in r:
        line = line.replace("\n","")
        id,label = line.split("\t")
        class_dict[string.upper(label)] = id
        #print label
    p.readline()
    out.write("contigID\tclassID\n")
    contig_dict = {}
    phylum_flag = 0
    for line in p:
        id = ""
        data = line.split("\t")
        contig_id = data[0]
        phylum = data[7].replace("\n","").replace(" ","")
        try:
            id = class_dict[string.upper(phylum)]
        except KeyError:
            print "phylum %s not found!"%phylum
            continue
        contig_dict[int(contig_id.split("_")[0])] = id
    keys = contig_dict.keys()
    keys.sort()
    for key in keys:
        out.write("contig_%d\t%s\n"%(key,contig_dict[key]))
    p.close()
    r.close()
    out.close()
    
       
