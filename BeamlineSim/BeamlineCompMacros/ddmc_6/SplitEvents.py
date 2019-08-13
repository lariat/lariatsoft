#! /usr/bin/env python
import time, os, shutil, sys, gc
import pprint
import subprocess
import datetime, json
import argparse
import warnings
import signal


parser = argparse.ArgumentParser()
parser.add_argument("filename", help="this is the name of the file you want to generate metadata for")
args = parser.parse_args()
filename  = args.filename

print filename

outname_short = filename[:-4]

print outname_short

thousand = 0

print thousand

lineCounter = 0
with open(filename) as fp:
    for line in fp:
        lineCounter+=1
        outname = outname_short +"_" + str(thousand) + ".list"
        target = open(outname, 'aw+')
        target.write(line)
        target.close()
        if not lineCounter % 2000:
            thousand+=1
	    print thousand
