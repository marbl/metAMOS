#bzip2 -d test4.sff.bz2
../initPipeline -l flx -f -m carsonella_pe_filt.fna -d test_sparse -i 3000:4000 
../runPipeline -c phylosift -p 8 -a sparseassembler -d test_sparse -k 55 -f Preprocess
