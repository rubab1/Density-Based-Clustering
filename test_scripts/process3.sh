#!/bin/sh
#
# PROCESS3.sh
#
# Run PROCESS3.m in Matlab batch mode.

# Shourov K. Chatterji
# shourov@ligo.caltech.edu
# 2006-Jul-24

# $Id:$

# change to analysis directory
cd `dirname $0`

# path to matlab executable
PATH=/lcdg/matlab_r13/bin:${PATH}

# run qtestcluster
echo "process3" | matlab -nosplash -nodesktop -nojvm

# echo newline
echo ""

# return
exit 0
