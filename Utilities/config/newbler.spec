#-m keep seq data in memory to speed up assembly
-m
#-a minimum contig length for 454AllContigs
-a 0
#-l minimum contig length for 454LargeContigs
-l 1000
#-ml minimum overlap length
-ml 60
#-mi minimum overlap % identity
-mi 90
#-ud treats each read separately, not grouping duplicates
-ud
