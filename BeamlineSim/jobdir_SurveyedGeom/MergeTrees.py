#!/usr/bin/env python
# Author: Jason St. John
# A very brief skeleton of a script to read in a file of TTrees, all assumed to be part of
# the same simulated particle interactions, and save them in a single new tree, indexed by
# the unique combinations of their Track and Event numbers. 
#
# Usage:
# ./MergeTrees.py <options> results.root 

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

parser = optparse.OptionParser("usage: %prog [options]<input file.ROOT> \n")
parser.add_option ('--o', dest='outfile', type='string',
                   default = 'MergedTree_test.root',
                   help="Output filename (ends .root).  This option is ignored.")
parser.add_option ('-T', dest='starterTree', type='string',
                   default = 'BigDisk',
                   help="The one TTree whose tracks will be iterated over. Effectively requires tracks present here.")
parser.add_option ('--spillsize', dest='spillsize', type='int',
                   default = 300000,
                   help="The number of G4BL events (particles launched at the target) per spill.")
parser.add_option ('--spillinterval', dest='spillinterval', type='float',
                   default = 60.0,
                   help="The duration of a sub-run. (seconds)")
parser.add_option ('-l', dest='keepitlocal', action="store_true", default=False,
                   help="Keep the output file in the same directory as the input file.")


options, args = parser.parse_args()
outfile = options.outfile
debug = False
starterTree = options.starterTree
spillsize = options.spillsize
spillinterval = options.spillinterval
keepitlocal   = options.keepitlocal
infile = args[0]

################################
# Timing offset generator
################################
BatchesPerOrbit = 7
BucketsPerBatch = 84
BucketsPerOrbit = BatchesPerOrbit * BucketsPerBatch

bucketcenterspacing = float("18.8e-9")
bucketwidth = float("2.2e-9")
batchlength = bucketcenterspacing * float(BucketsPerBatch)
orbitlength = batchlength * float(BatchesPerOrbit)
spillduration = 4.2
OrbitsInSpill = spillduration / orbitlength
filledbatches = (1,2,3,4,5,6) # (out of BatchesPerOrbit)

# Function to return a time during the spill, weighted to get the 
# time structure of the Fermilab Test Beam Facility's beam
def RandomOffsetSeconds ():
    BucketInBatch = random.randint(1,BucketsPerBatch-1)
    BatchInOrbit = random.choice(filledbatches)
    OrbitInSpill = random.randint(0,int(OrbitsInSpill))
    
    offset = random.gauss(0,bucketwidth)
    offset += bucketcenterspacing * float(BucketInBatch)
    offset += batchlength * BatchInOrbit
    offset += orbitlength * OrbitInSpill
    return offset


##############################################
# Check files exist
##############################################
if not os.path.isfile(infile) :
    print "Cannot find %s.  Exiting." % infile
    sys.exit()

infiles = glob.glob(infile)
base = ""  ## The path to the file.
for file in infiles:
    base = os.path.splitext(file)[0]
    #print file, base


infilename = infile
infile = ROOT.TFile(infile)

################################
# Set up ROOT
################################

#  Do ROOT initialization stuff here so that
#  my --help gets reported instead of ROOT --help.

ROOT.gROOT.SetBatch(0)
#ROOT.gStyle.SetOptStat(0)
#ROOT.gStyle.SetPalette(1)
#ROOT.gStyle.SetNdivisions(405,"x")

#####################################
# Declare & fill dictionaries
#####################################
INtuples = {}

#########################################
# Get into the file and see what's there
#########################################
infile.cd()
ROOT.gDirectory.cd()
abspath = base+".root:/VirtualDetector/"
print "|",abspath,"|"
ROOT.gDirectory.cd(abspath)
ROOT.gDirectory.ls()

dumtuple = ROOT.TNtuple() # Just need an instance of this class, it seems.

if debug: print "Looping over input file contents, getting trees."

# Read in all the single-detector TTree objects in the input file.
for key in ROOT.gDirectory.GetListOfKeys():
    if dumtuple.Class() == key.ReadObj().Class():
        if debug: print key.GetName(),
        if debug: print key.ReadObj().ClassName()
        INtuples[key.GetName()] =  key.ReadObj()

# To speed up finding tracks by EventID and TrackID, 
# build an index with these as the major and minor indices.
for name, tuple in INtuples.iteritems():
    if debug: print "Building an index for ",name,"...",
    tuple.BuildIndex("EventID","TrackID")
    if debug: print "Done."

# Lists of variable names and TTree names to use in loops.
vars = ('x','y','z','t','Px','Py','Pz','PDGid','ParentID','EventID','TrackID')
StartLine = ('StartLine',)
WCs = ('Det1', 'Det2', 'Det3', 'Det4')
Scints = ('TOFus', 'Halo', 'HaloHole', 'TOFdsHorz', 'TiWindow','BigDisk')
Punch  =  ('PunchUL', 'PunchLL', 'PunchUR', 'PunchLR')

### One dictionary to rule them all. ##
## Unfortunately, ROOT won't process a single line defining a single struct for all these; too long.  We integrate the by parts.
detsysts = {} 
detsysts['StartLine'] = StartLine
detsysts['WCs'] = WCs
detsysts['Scints'] = Scints
detsysts['Punch'] = Punch

# Invent some types of struct for holding stuff from the tree,
# one for each "system", 
# and tell ROOT about them.
for systname,syst in detsysts.iteritems():
    linetoprocess = "struct "+systname+"stuff { Int_t EventID; Int_t TrackID; "
    ## Be sure to do the ints before the floats ##
    for var in ('EventID', 'TrackID'):
        for det in syst:
            linetoprocess += "Int_t " + var + det + "; "
    for var in vars:
        if var == "EventID" or var == "TrackID": continue  #Don't re-do these ints!
        for det in syst:
            linetoprocess += "Float_t " + var + det + "; "
    linetoprocess += "};"
    if debug: print linetoprocess,"\n\n"
    ROOT.gROOT.ProcessLine(linetoprocess)

#Make one of each struct to use
structs = {}
for systname in detsysts.keys():
    structs[systname] = eval("ROOT."+systname+"stuff()")

#Hook the variables to their branches, so they take on the right values when an entry is fetched.
for systname,syst in detsysts.iteritems():
    if debug: print "\n Syst:",systname
    for var in vars:
        if debug: print "    now on :",var
        for det in syst:
            if debug: print "       det:",det
            tuple = INtuples[det]
            tuple.SetBranchAddress(var, ROOT.AddressOf(structs[systname],var+det))

#### The Output: A file with a tree in it. ####
infilepath = os.path.abspath(infilename)
infilename = os.path.basename(infilename)

if keepitlocal:
    outfilename = infilepath+"MergedAt"+starterTree+infilename
else:
    outfilename = "MergedAt"+starterTree+infilename
outfile = ROOT.TFile(outfilename,"RECREATE")

newTree = ROOT.TTree("EventTree","EventTree")
## Two important branches to include:
SpillID_p = array( 'i', [ 0 ] )
EventID_p = array( 'i', [ 0 ] )
TrackID_p = array( 'i', [ 0 ] )
newTree.Branch("SpillID",SpillID_p,"SpillID/I")
newTree.Branch("EventID",EventID_p,"EventID/I")
newTree.Branch("TrackID",TrackID_p,"TrackID/I")
## ...and all the others we want:
pointers = {}
if debug: print "                 Filling Pointers:"
for syst in detsysts.values():
    for det in syst:
        name = 'TrackPresent'+det
        if debug: print "                              ",name
        pointers[name] = array( 'B', [ 0 ] )
        newTree.Branch(name,pointers[name],name+"/O")
for var in vars:
    if var == "EventID" or var == "TrackID" or var == "SpillID": continue  # Only one unique combo of these, already done above.
    # Loop over every system, every detector in that system, for this variable
    for syst in detsysts.values():
        for det in syst:
            name = var+det
            if debug: print "                              ",name
            pointers[name] = array( 'd', [ 0 ] )
            newTree.Branch(name,pointers[name],name+"/d")

outfile.cd()

## How to remove the trigger-like requirement?
## Maybe we loop over everything which encountered the Halo paddle (including the HaloHole volume),
## and also add whatever additional entries can be found in the first WireChamber.

trackcount = 0
lastspill = 0
for ds_track in INtuples[starterTree]:
    trackcount += 1
    ## Unique identifiers for this track:
    (event, track) = (int(ds_track.EventID), int(ds_track.TrackID))
    # Development purposes: (Abbreviated run)
    #if trackcount > 1000: break

    ## What I hate about this is that the newTree looks completely uninvolved.
    ## Of course, it's tied in above, at newTree.Branch(name,pointers[name],name+"/f")
    SpillID_p[0] = 1 + (event/spillsize) # Arbitrarily group tracks by EventID into spills.
    EventID_p[0] = event
    TrackID_p[0] = track

    if SpillID_p[0] != lastspill: 
        if debug: print lastspill,"==> ",SpillID_p[0],"  (EventID:", EventID_p[0],")"
        lastspill += 1

    #Momentum cut? if ds_track.Pz < 1000. : continue
    for tuplename,tuple in INtuples.iteritems():
        for syst in detsysts.values():
            if tuplename not in syst: continue
            EntryNum = tuple.GetEntryWithIndex(event,track)

            if EntryNum == -1: #Track Not Found
                pointers['TrackPresent'+tuplename][0] = False
            else:
                pointers['TrackPresent'+tuplename][0] = True
            for var in vars:
                if var == 'EventID' or var == 'TrackID': continue
                ## Pointers are named like 'PzWC4' etc..
                vardet = var + tuplename
                if EntryNum == -1: #Track Not Found
                    pointers[vardet][0] = -123456789  ##Dummy Value
                    if debug: print "(",event,",",track,"): ",vardet,"       =-1"
                else: # Track found.
                    for systname,syst in detsysts.iteritems():
                        if tuplename in syst:
                            if debug: print "Filling ",vardet," in ",systname
                            if var == 't':  # Convert times to seconds & add an offset mimicking spill time profile
                                random.seed(event) #Same random offset for each event.
                                spilltimeoffset = spillinterval * float (SpillID_p[0])
                                pointers[vardet][0] = getattr(structs[systname],vardet)*1e-9 + RandomOffsetSeconds() + spilltimeoffset
                            else:
                                pointers[vardet][0] = getattr(structs[systname],vardet)
                                break #tuplename will only be in one syst, not more. q
    newTree.Fill()
print trackcount, "  total tracks in ",starterTree
print lastspill, " total spills (sub-runs)"
outfile.cd()
newTree.Write()
outfile.Close()
#-------------- /Final Beam --------------------#

infile.Close()
if not keepitlocal: commands.getoutput("ln -s /lariat/data/users/MYUSERNAME/MergedAt"+starterTree+infilename+" MergedAt"+starterTree+infilename)
print "Here ya go: ",
commands.getoutput("ls -ltra")
