#/bin/sh
../initPipeline -q -1 SRR987657 -d test_sra -W iMetAMOS
../runPipeline -d test_sra -p 16
