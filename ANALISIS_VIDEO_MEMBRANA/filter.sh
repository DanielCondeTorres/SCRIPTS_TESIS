#!/bin/bash

# Script to help to filter xtc trajectory
# this code is provided "as is" without warranty of any kind.  
# SBCB, Oxford, March 2014 

# script to filter the data
# arg1: xtc file
# arg2: tpr file
# arg3: filter size, 5, 10, 20 etc ..

gmx filter -f  $1 -nf $3 -all -ol lowpass.xtc
echo 0 | gmx trjconv -f lowpass.xtc -pbc atom -o system.xtc -s $2
echo 0 | gmx trjconv -f system.xtc -pbc whole -o low_pass.xtc -s $2
filename="filter_$3.xtc"
mv low_pass.xtc $filename
rm system.xtc
