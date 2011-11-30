../initPipeline -f -m carsonella_pe_filt.fna -d test_amphora  -i 150:450
../runPipeline -c amphora -p 4 -d test_amphora -k 45 -t -f Preprocess
