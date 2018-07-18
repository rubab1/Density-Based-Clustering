#!/bin/sh
#
# PROCESS2.sh
#
# Post processes raw trigger files from QTESTCLUSTER to remove all but the
# most significant trigger on a 1 second time scale and test for detected
# injections on a 0.5 second time scale.

# Shourov K. Chatterji
# shourov@ligo.caltech.edu
# 2006-Jul-24

# $Id:$

# change to analysis directory
cd `dirname $0`

# path to post processing utilities
PATH=../../postprocess/bin:${PATH}

# desired minimum separation in time between triggers
minimumSeparation=1.0

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

  # thin single triggers
  sparse ${minimumSeparation} \
         ${waveform}_singles.txt \
         ${waveform}_singles_sparse.txt

  # thin cluster triggers
  sparse ${minimumSeparation} \
         ${waveform}_clusters.txt \
         ${waveform}_clusters_sparse.txt

  # if waveform is an injection
  if [ -f ${waveform}_injections.txt ]; then

    # identify single detections of injections
    detect ${detectionDuration} \
           ${waveform}_injections.txt ${waveform}_singles_sparse.txt \
           ${waveform}_singles_detected.txt ${waveform}_singles_measured.txt

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
