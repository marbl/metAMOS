#bzip2 -d test4.sff.bz2
../initPipeline -l flx -f -m carsonella_pe_filt.fna -c test.asm -c carsonella_reference.fna -d test_isolate -i 3000:4000 
../runPipeline -u -r -v -I -c kraken -p 8 -a masurca,soap,velvet,ray -d test_isolate -k 51 -t -n FunctionalAnnotation -f MapReads,Postprocess -z phylum
