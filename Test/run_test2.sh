#/bin/sh
../createProject.py -m MA_test2.filt2.fna -d test2 -f -i 1300:1700
../runPipeline.py -c fcp -p 8 -d test2 -k 31 -f Scaffold,Propagate,Classify -v
