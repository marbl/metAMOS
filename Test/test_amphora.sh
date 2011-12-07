../initPipeline -q -1 test1.fq -2 test2.fq  -d test_amphora  -i 150:450
../runPipeline -a soap -c amphora -p 4 -d test_amphora -k 45 -t -f Preprocess
