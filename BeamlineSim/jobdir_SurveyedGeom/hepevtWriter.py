#!/usr/bin/env python
# Author: Jason St. John
# A very brief skeleton of a script to read in a file of TTrees representing particle interactions in G4BL, all assumed to be part of
# the same simulated particle interactions such that EventID and TrackID combinations are unique, and save them in new trees, one for
# each <spillsize> range of EventID values.  Time values (t) will be offset by <spillinterval> * SpillID, plus a random offset chosen 
# from a distribution which mimics the Fermilab Test Beam Facility's 4.2 seconds of 53 MHz beam.
#
# Usage:
# python hepevtWriter.py <options> FileOfManySpillTrees.root 

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
parser = optparse.OptionParser("usage: %prog [options]<input file.ROOT> \n")
parser.add_option ('-o', dest='outfile', type='string',
                   default = 'SortedTrees.root',
                   help="Output filename (ends .root).  This option is ignored.")
parser.add_option ('--maxspill', dest='maxspill', type='int',
                   default = -1,
                   help="Abbreviate processing to this many spill trees.")
parser.add_option ('--spillinterval', dest='spillinterval', type='float',
                   default = 60.0,
                   help="The duration of a sub-run. (seconds)")
parser.add_option ('--driftinterval', dest='driftinterval', type='float',
                   default = 0.00039321600,
                   help="The duration of a TPC drift. (seconds)")
parser.add_option ('-v', dest='debug', action="store_true", default=False,
                   help="Turn on verbose debugging.")
parser.add_option ('--test', dest='test', action="store_true", default=False,
                   help="Abbreviate to first spill, first 0.2 seconds")


options, args = parser.parse_args()
outfile = options.outfile
maxspill = options.maxspill
spillinterval = options.spillinterval
driftinterval = options.driftinterval
debug = options.debug
test = options.test
infile = args[0]

## Constants and such.
OneDrift = 0.0003932160 # 128 ns * 3072 samples
neutrals = (22, 311, 130, 310, 2112)
coordshift = {} # mm.  Add to target coords to get to TPC coords.
coordshift['x'] = 1473.59
coordshift['y'] = 23.42
coordshift['z'] = -8559.57

# define dynamically a python class containing root Leaves objects
class PyListOfLeaves(dict) :
    pass

massesbyPDG = {} # In GeV/c/c
massesbyPDG[11]   = 0.000511
massesbyPDG[13]   = 0.105658
massesbyPDG[22]   = 0.
massesbyPDG[211]  = 0.139570
massesbyPDG[2212] = 0.938272
massesbyPDG[321]  = 0.493667

# gimmestr()
# 
# Return the hepevt line: (GeV, ns, cm)
# 1 [pdg] 0 0 0 0 px py pz E m x y z t
def gimmestr(pile, tzero):
    pdg = int(pile.PDGidStartLine.GetValue())                 #0
    if not abs(pdg) in massesbyPDG.keys():                   
        exit ('Unknown particle species ',pdg)               
    Px = pile.PxStartLine.GetValue()/1000.                    #1 MeV/c --> GeV/c
    Py = pile.PyStartLine.GetValue()/1000.                    #2 MeV/c --> GeV/c
    Pz = pile.PzStartLine.GetValue()/1000.                    #3 MeV/c --> GeV/c
    mass = massesbyPDG[abs(pile.PDGidStartLine.GetValue())] 
    energy = pow( pow(Px,2.0) + pow(Py,2.0) + pow(Pz,2.0) + pow(mass,2.0) , 0.5)
    E = energy                                #4
    m = mass                                  #5
    x = ( pile.xStartLine.GetValue()+coordshift['x'] )/10.0   #6
    y = ( pile.yStartLine.GetValue()+coordshift['y'] )/10.0   #7
    z = ( pile.zStartLine.GetValue()+coordshift['z'] )/10.0   #8
    t = ( pile.tStartLine.GetValue() - tzero )* 1.0e9        #9 NTS: Shift zero to trigger time. Pre-extract?
    hepstr = '1 {0:d} 0 0 0 0 {1} {2} {3} {4} {5} {6} {7} {8} {9}'.format(pdg, Px, Py, Pz, E, m, x, y, z, t)
    return hepstr

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
for key in ROOT.gDirectory.GetListOfKeys():
    if dumtuple.Class() == key.ReadObj().Class():
        name = key.GetName()
        if debug: print name,
        if debug: print key.ReadObj().ClassName()
        # Important to access only the most recent name cycle. (DAMNIT, ROOT!)
        if name in treenames: continue
        treenames.append(name) # Remember this one for later
        InputSpillTrees[name] =  key.ReadObj() # Store pointer to this TTree in the dictionary

# To speed up finding tracks by EventID and TrackID, 
# build an index with these as the major and minor indices.
for name, tuple in InputSpillTrees.iteritems():
    print "Building an index for ",name,"...",
    tuple.BuildIndex("tStartLine")
    print "Done."

#### The Output: A file with a tree in it. ####
infilepath = os.path.abspath(infilename)
infilename = os.path.basename(infilename)

outfilename = "SortedEvents"+infilename+".txt"
outfile = open(outfilename,'w')

# Initialize some counters
trackcount = 0
lastspill = 0
spillcount = 0
die = False
# Loop over input TTree objects
for intree in InputSpillTrees.values():
    n_entries = intree.GetEntriesFast()
    print "Starting ",intree.GetName()," with ",n_entries," entries."
    spillcount += 1
    # Extract the spill number from the TTree's name. Comes after the _Spill
    spill = int(intree.GetName().split("_Spill")[1])

    # Gonna need a handy little class hooked to this tree:
    # Define dynamically a python class containing root Leaves objects
    # get all leaves 
    leaves = intree.GetListOfLeaves()
    # create an instance
    pyl = PyListOfLeaves()
    # add leaves as attributes
    for i in range(0,leaves.GetEntries() ) :
        leaf = leaves.At(i)
        name = leaf.GetName()
        # add dynamically attribute to the baby class
        pyl.__setattr__(name,leaf)

    # Build a dictionary of index values, with their times as the keys
    allentriesbytime = {}
    # Also some useful other lists:
    entrytimes = []
    triggertimes = []
    triggerentrynums = []
    
    # First Loop over this tree: Get the entry numbers and tStartLine (if defined)
    for n in xrange(0, n_entries):
        intree.GetEntry(n) # Fill pyl with values from the entry at index n
        timeExists = pyl.TrackPresentStartLine.GetValue()
        time = float(pyl.tStartLine.GetValue())
        if not timeExists: continue
        if debug: print n,":",time
        entrytimes.append(time)

        # Make sure there's a list of entry numbers in the dictionary for this time
        if time not in allentriesbytime.keys(): allentriesbytime[time] = []
        allentriesbytime[time].append(n)

    # What did we get?  Any non-unique index values?
    timecount = len(entrytimes)
    uniqcount = len(set(entrytimes))
    if debug: print "Time count:",len(allentriesbytime.keys())," from a spill of size:", n_entries
    if debug: print "With entry.TrackPresentStartLine:",intree.GetEntries("TrackPresentStartLine")
    if debug: print "total time count: ",len(entrytimes)," (unique:",len(set(entrytimes)),", non-unique:",100.*(float(timecount)-float(uniqcount))/float(timecount),"%)"

    # Second Loop over the tree:
    # Visit entries in order by their times, and check for triggering particles.
    LastTriggerTime = float(-1)
    firsttime = 0.0
    for time in sorted(allentriesbytime.keys()):
        if die: break
        if LastTriggerTime == -1.: firsttime = time
        #Get the delta_t and update LastTriggerTime
        delta_t = time-LastTriggerTime

        # For this time value, visit one or more entries n:
        for n in allentriesbytime[time]:
            if die: break
            ret = intree.GetEntry(n)
            if debug: print "time: ",time,"   delta_t:",delta_t,"   n:",n            
            # Was this a triggering particle?
            pdg = pyl.PDGidStartLine.GetValue()
            if pdg in neutrals: continue # Photons? Don't care.
            if (LastTriggerTime == -1 or (LastTriggerTime > -1 and delta_t > 2.0*driftinterval) and
                # Wire Chambers
                pyl.TrackPresentDet1.GetValue() and 
                (pyl.TrackPresentDet2.GetValue() or pyl.TrackPresentDet3.GetValue()) and
                pyl.TrackPresentDet4.GetValue() and 
                # Time of Flight, too
                pyl.TrackPresentTOFus.GetValue() and pyl.TrackPresentTOFdsHorz.GetValue()):
                
                #..if so, add it to the list
                triggertimes.append(time) 
                triggerentrynums.append(n) 
                LastTriggerTime = time
                print "Triggering: ",n,":",time
                if test and time - firsttime > 0.2: 
                    die = True
                    break

    # Third and final loop: 
    #Collect the hepevt fields for all particles which might give signals in the triggered events.
    for time in triggertimes:
        delta_t = 0.
        offset = 0
        print "\n Collecting tracks around trigger at ",time
        # OK, where are we in the sorted list of particle times?
        triggerTimeIndex = sorted(allentriesbytime.keys()).index(time)
        checkindex = triggerTimeIndex - offset
        entry_t = sorted(allentriesbytime.keys())[checkindex]
        delta_t = abs(time - entry_t)
        if debug: print "offset: ",offset,"entry_t:",entry_t,"    delta_t:",delta_t

        # Print the hepevt lines: (GeV, ns, cm)
        # 1 [pdg] 0 0 0 0 px py pz E m x y z t
        # First at the time of the triggering particle:
        for n in allentriesbytime[time]:
            ret = intree.GetEntry(n)
            if ret == -1: exit ('No entry #'+str(n))
            txtstr = gimmestr(pyl, time)
            print txtstr

        # At which times were there particles BEFORE this trigger, close enough to put energy in the TPC?
        while (True):
            offset +=1
            checkindex = triggerTimeIndex - offset
            if checkindex < 0: break #Not too early!

            # Get the event time at this nearby index
            entry_t = sorted(allentriesbytime.keys())[checkindex]

            # How far away in time from our triggering particle? 
            delta_t = abs(time - entry_t)
            if delta_t > driftinterval: break # Don't go more than one driftinterval.

            if debug: print "offset: -",offset,"entry_t:",entry_t,"    delta_t:",delta_t

            # Loop over any entries which had this time value
            for n in allentriesbytime[entry_t]:
                if debug: print "Getting n:",n
                ret = intree.GetEntry(n) # Set pyl's pointers to this entry's leaf values
                if ret == -1: exit ('No entry #'+str(n))
                txtstr = gimmestr(pyl, time)
                print txtstr

        # At which times were there particles AFTER this trigger, close enough to put energy in the TPC?
        delta_t = 0.
        offset = 0
        while (True):
            offset +=1
            checkindex = triggerTimeIndex + offset
            if checkindex > len(allentriesbytime): break #Not too late!

            # Get the event time at this nearby index
            entry_t = sorted(allentriesbytime.keys())[checkindex]

            # How far away in time from our triggering particle? 
            delta_t = abs(time - entry_t)
            if delta_t > driftinterval: break # Don't go more than one driftinterval.

            if debug: print "offset: +",offset,"entry_t:",entry_t,"    delta_t:",delta_t

            # Loop over any entries which had this time value
            for n in allentriesbytime[entry_t]:
                if debug: print "Getting n:",n
                ret = intree.GetEntry(n) # Set pyl's pointers to this entry's leaf values
                if ret == -1: exit ('No entry #'+str(n))
                txtstr = gimmestr(pyl, time)
                print txtstr

outfile.Close()
infile.Close()

