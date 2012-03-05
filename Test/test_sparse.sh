#bzip2 -d test4.sff.bz2
../initPipeline -l flx -q -m carsonella_pe_filt.fq -d test_sparse -i 3000:4000 
../runPipeline -c phmmer -p 8 -a sparseassembler -d test_sparse -k 55 -t -f Preprocess
