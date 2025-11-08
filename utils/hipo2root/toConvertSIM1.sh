#!/bin/csh -f
module use /scigroup/cvmfs/hallb/clas12/sw/modulefiles
module load clas12root
clas12root -q -b ana12GeVShortFCQA.C  --in=pass1_rga_inclusive_runs_valera_test.dat >> ../../logs/pass1_rga_inclusive_runs_valera_test.log
