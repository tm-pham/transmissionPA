#! /bin/bash
cd /home/julius_id/tpham/BedwiseExpMCMC
g++ -I/home/julius_id/tpham/library/boost -L/home/julius_id/tpham/library/boost readcsvfileto2Darrayofint.cpp relabel.cpp loglikelihood.cpp changeloglandI.cpp createusefuldata.cpp truepositiveandfalsenegative.cpp samplefromnormal.cpp pacquisition.cpp environment.cpp mainNonHamBed.cpp -o mainNonHamBed
./mainNonHamBed

