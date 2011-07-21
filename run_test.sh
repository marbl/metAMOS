#/bin/sh
python createProject.py -m carsonella_pe_filt.fna -d test1 -f -i 150,450
python runPipeline.py -c fcp -p 4 -d test1 -k 61 
