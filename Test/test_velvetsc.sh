#bzip2 -d test4.sff.bz2
../initPipeline -l flx -f -m carsonella_pe_filt.fna -d test_velvetsc -i 3000:4000 
../runPipeline -c fcp -p 8 -a velvet-sc -d test_velvetsc -k 55 -f Preprocess
