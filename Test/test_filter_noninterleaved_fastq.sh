../initPipeline -q -1 test1.fq -2 test2.fq -d test12 -i 150:450
../runPipeline -c fcp -p 4 -d test12 -k 61 -t -f Preprocess
