#bzip2 -d test4.sff.bz2
../initPipeline -q -1 lambda.sim.pac.fastq -d test_pbcr
../runPipeline -u -r -v -c kraken -t pbcr -p 16 -a soap2,soap,velvet,ca -d test_pbcr -n FunctionalAnnotation -f Postprocess -z phylum
