#!/bin/bash
#
# Script to illustrate running batch jobs and passing in arguments.
#
# 
# This script assumes that the following has been run successfully:
# scons co=1 b=GccOpt projects/CryptCellSeparation/test/TestCryptSingleRun.hpp
#


LABELED_RATIO=0.33;


SEPARATION=1.0;
SEPARATION_STEP=0.1; 
MAX_SEPARAITION=4.0;

for SEPARATION in $(seq 0 0.1 4)
do
	echo "  SEPARATION " ${SEPARATION};
	# NB "nice -20" gives the jobs low priority (good if they are going to dominate the server and no slower if nothing else is going on)
	# ">" directs std::cout to the file.
	# "2>&1" directs std::cerr to the same place.
	# "&" on the end lets the script carry on and not wait until this has finished.
	nice -20 ../build/optimised/TestCryptSingleRunRunner -labeled_ratio ${LABELED_RATIO} -separation_multiplier ${SEPARATION} > output/Run_${i}_${LABELED_RATIO}_${SEPARATION}_Output.txt 2>&1 &
done

echo "Jobs submitted"


