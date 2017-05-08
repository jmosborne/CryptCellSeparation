#!/bin/bash
#
# Script to illustrate running batch jobs and passing in arguments.
#
# 
# This script assumes that the following has been run successfully:
# scons co=1 b=GccOpt ts=projects/CryptProliferation2013/test/TestCryptTakeoverProbabilityLiteratePaper.hpp
#

num_sims=10;
num_sweeps=10;

for (( i=0 ; i<${num_sweeps} ; i++))
do
	start_sim=`expr $i \* $num_sims`;
	echo $start_sim

	# ">" directs std::cout to the file.
	# "2>&1" directs std::cerr to the same place.
	# "&" on the end lets the script carry on and not wait until this has finished.
	nice -20 ../build/optimised/TestClonalCryptRunRunner -num_runs $num_sims -run_index $start_sim > output/Run_${i}_Output.txt 2>&1 &
done

echo "Jobs submitted"
