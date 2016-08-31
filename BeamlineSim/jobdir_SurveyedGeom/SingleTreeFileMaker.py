#!/usr/bin/env python
# Author: Jason St. John
# A very brief skeleton of a script to read in a file of TTrees and save certain ones to their own new files.
#
# Usage:
# python SingleTreeFileMaker.py <options> FileOfManySpillTrees.root 

import ROOT
import sys
import optparse
import commands
import os
import glob
from decimal import Decimal
import re
from array import *
import math
from ctypes import *
import random
import ctypes
import copy

parser = optparse.OptionParser("usage: %prog [options]<input file.ROOT> \n")
parser.add_option ('--maxspill', dest='maxspill', type='int',
                   default = -1,
                   help="Abbreviate processing to this many spill trees.")
parser.add_option ('--firstspill', dest='firstspill', type='int',
                   default = -1,
                   help="Spill number to process first.")
parser.add_option ('--lastspill', dest='lastspill', type='int',
                   default = -1,
                   help="Spill number to process last.")
parser.add_option ('-v', dest='debug', action="store_true", default=False,
                   help="Turn on verbose debugging.")

options, args = parser.parse_args()
maxspill = options.maxspill
firstspill = options.firstspill
lastspill = options.lastspill
debug = options.debug
infile = args[0]

####################################################
# Check files exist, get some strings, make a TFile
####################################################
if not os.path.isfile(infile) :
    print "Cannot find %s.  Exiting." % infile
    sys.exit()

infilename = infile
infiles = glob.glob(infile)
infilepath = ""  ## The path to the input file.
base = '' ## The input file name
for file in infiles:
    infilepath = os.path.splitext(file)[0]
    base = os.path.basename(file)
infile = ROOT.TFile(infile)

################################
# Set up ROOT
################################

#  Do ROOT initialization stuff here so that
#  my --help gets reported instead of ROOT --help.

ROOT.gROOT.SetBatch(0)  # Don't draw things. Slows us down.
#ROOT.gStyle.SetOptStat(0)
#ROOT.gStyle.SetPalette(1)
#ROOT.gStyle.SetNdivisions(405,"x")

#########################################
# Get into the file and see what's there
#########################################

InputSpillTrees = {}
infile.cd()
ROOT.gDirectory.cd()
abspath = base+":/"
if debug: print "|",abspath,"|"
ROOT.gDirectory.cd(abspath)
if debug: ROOT.gDirectory.ls()

dumtuple = ROOT.TTree() # Just need an instance of this class, it seems.

print "Looping over input file contents, getting spill trees."
# Read in all the single-detector TTree objects in the input file.
treenames = []
spillcount = 0
for key in ROOT.gDirectory.GetListOfKeys():
    if dumtuple.Class() == key.ReadObj().Class():
        name = key.GetName()
        if debug: print name,
        if debug: print key.ReadObj().ClassName()

        # Important to access only the most recent name cycle. (DAMNIT, ROOT!)
        if name in treenames: continue        
        treenames.append(name) # Remember this one for later

        # Extract the spill number from the spill tree's name
        spill = int(name.split("_Spill")[1])

        # Skip anything out of range (if range bounds are given)
        if firstspill > 0 and spill < firstspill : continue
        if lastspill > 0 and spill > lastspill : continue

        # OK, this is one to keep a pointer to.
        InputSpillTrees[spill] =  key.ReadObj() # Store pointer to this TTree in the dictionary
        spillcount += 1
        if maxspill > 0 and spillcount > maxspill: break

#### The Output: A file with a tree in it. ####
infilepath = os.path.abspath(infilename)
infilename = os.path.basename(infilename)

# Initialize a handy list
spillnums = []

# Loop over input TTree objects
for spill, intree in InputSpillTrees.iteritems():
    # Extract the spill number from the TTree's name. Comes after the _Spill
    print 'Spill:',spill,' ',intree.GetName()

    if firstspill > 0 and spill < firstspill: continue
    if debug: print '        Not before spill',firstspill
    if lastspill > 0 and spill > lastspill: continue
    if debug: print '        Not after spill',lastspill

    n_entries = intree.GetEntriesFast()
    if debug: print "Starting ",intree.GetName()," with ",n_entries," entries."
    spillnums.append(spill)
    outfilename = infilename.replace('.root','_OnlySpilltree'+str(spill)+'.root')
    outfile = ROOT.TFile(outfilename, 'RECREATE')

    outtree = intree.Clone()
    outfile.cd()
    outtree.Write()
    outfile.Close()

infile.Close()

