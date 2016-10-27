#!/bin/bash
#
# Script to illustrate running batch jobs and passing in arguments.
#
# 
# This script assumes that the following has been run successfully:
# scons co=1 b=GccOpt projects/CryptCellSepaSEPARATIONn/test/TestCryptSingleRun.hpp
#


LABELED_RATIO=0.33;


SEPARATION=1.0;
SEPARATION_STEP=0.1; 
MAX_SEPARAITION=4.0;


# SEPARATION[0]="1"
# SEPARATION[1]="1.25"
# SEPARATION[2]="1.5"
# SEPARATION[3]="1.75"
# SEPARATION[4]="2"
# SEPARATION[5]="2.25"
# SEPARATION[6]="2.5"
# SEPARATION[7]="2.75"
# SEPARATION[8]="3"
# SEPARATION[9]="3.25"
# SEPARATION[10]="3.5"
# SEPARATION[11]="3.75"
# SEPARATION[12]="4"

# SEPARATION[0]="1"
# SEPARATION[1]="1.5"
# SEPARATION[2]="2"
# SEPARATION[3]="2.5"
# SEPARATION[4]="3"
# SEPARATION[5]="3.5"
# SEPARATION[6]="4"

# SEPARATION[0]="1.25"
# SEPARATION[1]="1.75"
# SEPARATION[2]="2.25"
# SEPARATION[3]="2.75"
# SEPARATION[4]="3.25"
# SEPARATION[5]="3.75"




for SEPARATION in $(seq 0 0.1 4)
do
	echo "  SEPARATION " ${SEPARATION};
	# NB "nice -20" gives the jobs low priority (good if they are going to dominate the server and no slower if nothing else is going on)
	# ">" directs std::cout to the file.
	# "2>&1" directs std::cerr to the same place.
	# "&" on the end lets the script carry on and not wait until this has finished.
	nice -20 ../build/optimised/TestCryptSingleRunRunner -labeled_ratio ${LABELED_RATIO} -separation_multiplier ${SEPARATION} > output/Run_${i}_${LABELED_RATIO}_${SEPARATION}_Output.txt 2>&1 &
	#nice -20 ../build/optimised/TestCryptSingleRunRunner -sim_index $i -labeled_SEPARATION $SEPARATION > output/Run_${i}_${SEPARATION[$j]}_Output.txt 2>&1 &
done

echo "Jobs submitted"


