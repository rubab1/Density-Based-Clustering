#!/bin/sh
#
# PROCESS1.5.sh
#
# Concatanates singles and clusters files, to combine single triggers with clustered
# triggers. It also keeps a backup of the original singles and clusters truggers files.
#
# CAUTION: If this script is run multiples times in the same directory it will replace 
# backup with concatanated files, and put redundant singles and clusters entries in the
# processed files. 

# Rubab M. Khan
# rmk2109@columbia.edu
#
# Shourov K. Chatterji
# shourov@ligo.caltech.edu
# 2006-Aug-06

# $Id:$

# change to analysis directory
cd `dirname $0`

# path to post processing utilities
PATH=../../postprocess/bin:${PATH}

# determine list of waveforms to process
waveforms=`/bin/ls -1 *_singles.txt | sed -e 's|_singles.txt||'`

# begin loop over waveforms
for waveform in ${waveforms}; do

  # display status
  echo "copying ${waveform}..."

  cp ${waveform}_clusters.txt ${waveform}_clusters.orig.txt

  cat ${waveform}_singles.txt ${waveform}_clusters.orig.txt > ${waveform}_clusters.txt

# end loop over waveforms
done

# return
exit 0

