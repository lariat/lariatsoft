#!/usr/bin/env python
# Author: Jonathan Paley
#
# ./FindFhiclFile.py [file]
#
# Searches through FHICL_FILE_PATH to find 1st instance of the file.

import sys
import os

if len(sys.argv) < 2:
    print "Missing filename."
    exit(1)

fileName = sys.argv[1]

ffp = os.environ.get('FHICL_FILE_PATH')

dirList = ffp.split(':')

for dir in dirList:
    tf = dir + "/" + fileName
    if (os.path.isfile(tf)):
        print tf
        exit(0)

