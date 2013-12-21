#bzip2 -d test4.sff.bz2
../initPipeline -l flx -f -m carsonella_pe_filt.fna -d test_kraken -i 3000:4000
../runPipeline -u -r -v -c kraken -t eautils -p 8 -a soap2,soap,velvet -d test_kraken -k 55 -n FunctionalAnnotation -f Postprocess -z phylum
