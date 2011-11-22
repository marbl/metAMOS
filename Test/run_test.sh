#/bin/sh
../initPipeline -f -m carsonella_pe_filt.fna -d test1  -i 150:450
../runPipeline -c metaphyler -p 15 -d test1 -k 45
