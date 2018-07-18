#!/bin/sh
#
# PROCESS2.5.sh
#
# Run this after running process2.

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

# desired minimum separation in time between triggers
minimumSeparation=0

# duration of detection window
detectionDuration=0.5

# write detection duration to file
echo ${detectionDuration} >detection.txt

# determine list of waveforms to process
waveforms=`/bin/ls -1 *_singles.txt | sed -e 's|_singles.txt||'`

# begin loop over waveforms
for waveform in ${waveforms}; do

  # display status
  echo "processing ${waveform}..."

  cp ${waveform}_clusters_sparse.txt ${waveform}_clusters_sparse1.txt

  cat ${waveform}_singles_sparse.txt ${waveform}_clusters_sparse1.txt > ${waveform}_clusters_sparse2.txt

  # thin cluster triggers
  sparse ${minimumSeparation} \
         ${waveform}_clusters_sparse2.txt \
         ${waveform}_clusters_sparse.txt

  # if waveform is an injection
  if [ -f ${waveform}_injections.txt ]; then

    # identify cluster detections of injections
    detect ${detectionDuration} \
           ${waveform}_injections.txt ${waveform}_clusters_sparse.txt \
           ${waveform}_clusters_detected.txt ${waveform}_clusters_measured.txt

  # otherwise continue
  fi

# end loop over waveforms
done

# return
exit 0
