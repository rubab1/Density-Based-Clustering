#!/bin/sh
#
# PROCESS1.sh
#
# Run QTESTCLUSTER executable for the current analysis directory.

# Shourov K. Chatterji
# shourov@ligo.caltech.edu
# 2006-Jul-24

# $Id:$

# change to analysis directory
cd `dirname $0`

# path to qtestcluster executable
PATH=../../qtestcluster/bin:${PATH}

# run qtestcluster
qtestcluster.sh >log.txt 2>&1

# return
exit 0
