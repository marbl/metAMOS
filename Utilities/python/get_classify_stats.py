#!/usr/bin/python

import sys
import os

contigs_by_class = { }
origContigsByClass = { }
origClassifiedCount = 0
classifiedCount = 0
id_class = { }
id_class["0"] = "UNKNOWN"

orig_class_file = open(sys.argv[1])
class_file = open(sys.argv[2])
class_key = open(sys.argv[3])
#pass outdir as argument
out_dir = sys.argv[4]
orig_out = open(out_dir + os.sep + sys.argv[5], 'w')
out = open(out_dir + os.sep + sys.argv[6], 'w')

# parse in key file
for line in class_key:
    line = line.strip()
    fields = line.split("\t")
    # f1 is id, f2 is class name

    if len(fields) != 2:
        print "Error in file format\n"
    else:
        id_class[fields[0]] = fields[1]

# parse original file to identity ambiguous assignment (which is one more than max previous ID)
maxClassID = 0;
for line in orig_class_file:
    line = line.strip()
    fields = line.split()
    # f1 is contig, f2 is class

    if len(fields) != 2:
        print "Error in file format\n"

    elif maxClassID < int(fields[1]):
       maxClassID = int(fields[1])

    if origContigsByClass.has_key(fields[1]):
       origContigsByClass[fields[1]]+=1
    else:
       origContigsByClass[fields[1]] = 1
    origClassifiedCount += 1

id_class[str(maxClassID+1)] = "AMBIGUOUS"

# parse contig class file
for line in class_file:
    line = line.strip()
    fields = line.split()
    # f1 is contig, f2 is class

    if len(fields) != 2:
        print "Error in file format\n"
    elif contigs_by_class.has_key(fields[1]):
        contigs_by_class[fields[1]] += 1
    else:
        contigs_by_class[fields[1]] = 1

    if fields[1] > 0:
       classifiedCount += 1

# output stats
# todo: add info on ORFs and read counts
summary = []
summary.append("<p>Originally classified contigs:")
summary.append("<table>")
for key in origContigsByClass:
   class_name = id_class[key]
   summary.append("<tr>")
   summary.append("<td align=\"left\">%s</td><td align=\"right\">%d</td><td align=\"right\">%3.2f%%</td>"%(class_name, origContigsByClass[key], origContigsByClass[key]/float(origClassifiedCount)*100))
   summary.append("</tr>")
summary.append("</table>")
summary.append("Total classified: %d </p>"%(origClassifiedCount))
orig_out.write("var classifyHTML ='%s'"%("\\\n".join(summary)))
summary = []

summary.append("<p>Propagate classified contigs:")
summary.append("<table>")
for key in contigs_by_class:
    class_name = id_class[key]
    summary.append("<tr>")
    summary.append("<td align=\"left\">%s</td><td align=\"right\">%d</td><td align=\"right\">%3.2f%%</td>"%(class_name, contigs_by_class[key], contigs_by_class[key]/float(classifiedCount)*100))
    summary.append("</tr>")
summary.append("</table>")
additional = classifiedCount - origClassifiedCount 
if additional >= 0:
   summary.append("Total additional classified contigs: %d</p>"%(additional))
else:
   summary.append("Total contigs classified as unknown from known: %d</p>"%(abs(additional)))
out.write("var propagateHTML='%s'"%("\\\n".join(summary)))
