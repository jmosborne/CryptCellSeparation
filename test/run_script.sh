#!/bin/bash
#
# Script to illustrate running batch jobs and passing in arguments.
#
# 
# This script assumes that the following has been run successfully:
# scons co=1 b=GccOpt projects/CryptCellSeparation/test/TestCryptSingleRun.hpp
#


num_sims=1;

RATIO[0]="0.0"
RATIO[1]="0.5"
RATIO[2]="1.0"

for (( i=0 ; i<${num_sims} ; i++))
do
    echo "Run " $i;
    for (( j=0 ; j<${#RATIO[*]} ; j++))
    do
    	echo "  Ratio " ${RATIO[$j]};
    	# NB "nice -20" gives the jobs low priority (good if they are going to dominate the server and no slower if nothing else is going on)
    	# ">" directs std::cout to the file.
    	# "2>&1" directs std::cerr to the same place.
    	# "&" on the end lets the script carry on and not wait until this has finished.
    	nice -20 ../build/optimised/TestCryptSingleRunRunner -labeled_ratio ${RATIO[$j]} > output/Run_${i}_${RATIO[$j]}_Output.txt 2>&1 &
    	#nice -20 ../build/optimised/TestCryptSingleRunRunner -sim_index $i -labeled_ratio ${RATIO[$j]} > output/Run_${i}_${RATIO[$j]}_Output.txt 2>&1 &
    done
done

echo "Jobs submitted"


