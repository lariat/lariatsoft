#!/usr/bin/env python
# Author: Jason St. John
# A very brief skeleton of a script to read in a file of TTrees representing particle interactions in G4BL, all assumed to be part of
# the same simulated particle interactions such that EventID and TrackID combinations are unique, and save them in new trees, one for
# each <spillsize> range of EventID values.  Time values (t) will be offset by <spillinterval> * SpillID, plus a random offset chosen 
# from a distribution which mimics the Fermilab Test Beam Facility's 4.2 seconds of 53 MHz beam.
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
import pickle

parser = optparse.OptionParser("usage: %prog [options]<input file.ROOT> \n")
parser.add_option ('--o', dest='outfile', type='string',
                   default = 'MergedTree_test.root',
                   help="Output filename (ends .root).  This option is ignored.")
parser.add_option ('-T', dest='starterTree', type='string',
                   default = 'StartLine',
                   help="The one TTree whose tracks will be iterated over. Effectively requires tracks present here.")
parser.add_option ('--maxspill', dest='maxspill', type='int',
                   default = -1,
                   help="Abbreviate processing to this many spills.")
parser.add_option ('--spillsize', dest='spillsize', type='int',
                   default = 300000,
                   help="The number of G4BL events (particles launched at the target) per spill.")
parser.add_option ('--spillinterval', dest='spillinterval', type='float',
                   default = 60.0,
                   help="The duration of a sub-run. (seconds)")
parser.add_option ('--subspillnumber', dest='subspillnumber', type='int',
                   default = 1,
                   help="If processing the spill in parallel tranches, which subspill number is this.")
parser.add_option ('--subspillcount', dest='subspillcount', type='int',
                   default = 1,
                   help="If processing the spill in parallel tranches, how many subspills to assume.")
parser.add_option ('--gammafloor', dest='gammacutoff', type='float',
                   default = 0.5,
                   help="In the starterTree gammas less than this energy will be ignored. default: 0.5 (MeV)")
parser.add_option ('-l', dest='keepitlocal', action="store_true", default=False,
                   help="Keep the output file in the same directory as the input file.")
parser.add_option ('-v', dest='debug', action="store_true", default=False,
                   help="Turn on verbose debugging.")


options, args  = parser.parse_args()
outfile        = options.outfile
debug          = options.debug
starterTree    = options.starterTree
maxspill       = options.maxspill
spillsize      = options.spillsize
spillinterval  = options.spillinterval
subspillnumber = options.subspillnumber
subspillcount  = options.subspillcount
keepitlocal    = options.keepitlocal
gammacutoff    = options.gammacutoff
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
    subspillduration = spillduration/subspillcount
    subspilltimewindow_early = subspillnumber * subspillduration
    subspilltimewindow_late = subspilltimewindow_early + subspillduration
    offset = -1
    while offset < subspilltimewindow_early or offset >= subspilltimewindow_late:
        BucketInBatch = random.randint(1,BucketsPerBatch-1)
        BatchInOrbit = random.choice(filledbatches)
        OrbitInSpill = random.randint(0,int(OrbitsInSpill))
        
        offset = random.gauss(0,bucketwidth)
        offset += bucketcenterspacing * float(BucketInBatch)
        offset += batchlength * BatchInOrbit
        offset += orbitlength * OrbitInSpill
    # exit loop when we finally get the right range
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

ROOT.gROOT.SetBatch(0)  # Don't draw things. Slows us down.
#ROOT.gStyle.SetOptStat(0)
#ROOT.gStyle.SetPalette(1)
#ROOT.gStyle.SetNdivisions(405,"x")


#########################################
# Get into the file and see what's there
#########################################

INtuples = {}
infile.cd()
ROOT.gDirectory.cd()
abspath = base+".root:/VirtualDetector/"
if debug: print "|",abspath,"|"
ROOT.gDirectory.cd(abspath)
if debug: ROOT.gDirectory.ls()

dumtuple = ROOT.TNtuple() # Just need an instance of this class, it seems.

print "Looping over input file contents, getting trees."
# Read in all the single-detector TTree objects in the input file.
for key in ROOT.gDirectory.GetListOfKeys():
    if key.GetName() in INtuples.keys(): continue # Do not pick up further, lower-cycle versions of trees already grabbed. 
    if dumtuple.Class() == key.ReadObj().Class():
        if debug: print key.GetName(),
        if debug: print key.ReadObj().ClassName()
        INtuples[key.GetName()] =  key.ReadObj()

# To speed up finding tracks by EventID and TrackID, 
# build an index with these as the major and minor indices.
for name, tuple in INtuples.iteritems():
    print "Building an index for ",name,"...",
    tuple.BuildIndex("EventID","TrackID")
    print "Done."

# Lists of variable names and TTree names to use in loops.
vars = ('x','y','z','t','Px','Py','Pz','PDGid','ParentID','EventID','TrackID')
StartLine = ('StartLine',)
WCs = ('Det1', 'Det2', 'Det3', 'Det4')
Scints = ('TOFus', 'TOFds') # Horz removed

## One dictionary to rule them all. ##
## Unfortunately, ROOT won't process a single line defining a single struct for all these; too long.  
## So here we integrate them by parts.
detsysts = {} 
detsysts['StartLine'] = StartLine
detsysts['WCs'] = WCs
detsysts['Scints'] = Scints

# Invent some types of struct for holding stuff from the tree,
# one for each "system", 
# and tell ROOT about them.
print "Making structs for ROOT...",
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
print "Done.",

#Hook the variables to their branches, so they take on the right values when an entry is fetched.
for systname,syst in detsysts.iteritems():
    if debug: print "\n Syst:",systname
    for var in vars:
        if debug: print "    now on :",var
        for det in syst:
            if debug: print "       det:",det
            tuple = INtuples[det]
            tuple.SetBranchAddress(var, ROOT.AddressOf(structs[systname],var+det))
print " Also hooked structs to input TTree branches. It's like a pointer jungle in here."

#### The Output: A file with a tree in it. ####
infilepath = os.path.abspath(infilename)
infilename = os.path.basename(infilename)

if keepitlocal:
    outfilename = infilepath+"MergedAt"+starterTree+infilename
else:
    outfilename = "MergedAt"+starterTree+infilename
outfile = ROOT.TFile(outfilename,"RECREATE")

timeindexf = outfilename.replace('root','pickle')
newTrees = {} # dictionary of new TTrees to save to outfile
pointers = {} # dictionary of pointers, to please ROOT, which loves pointers for its TTrees

#
# AddTree: A worker function to instantiate a new TTree for a Spill
# and supply it and its variable pointers to the  relevant dictionaries. 
# Also hooks the branches to their pointers.
#
def AddTree(treedict, spillnum, pointers, timeindex):
    # Add a dictionary for the time index of this spill
    timeindex[spillnum] = {}

    print "Adding tree for spill ",spillnum
    treename = "EventTree_Spill"+str(spillnum)
    # Make the new TTree ans store it in the TTree dictionary
    treedict[spillnum] = ROOT.TTree(treename,treename)
    # Make ROOT happy with some pointers
    pointers[spillnum, "SpillID"] = array( 'i', [ 0 ] )
    pointers[spillnum, "EventID"] = array( 'i', [ 0 ] )
    pointers[spillnum, "TrackID"] = array( 'i', [ 0 ] )
    #...and hook the pointers to the branches.
    treedict[spillnum].Branch("SpillID",pointers[spillnum, "SpillID"],"SpillID/I")
    treedict[spillnum].Branch("EventID",pointers[spillnum, "EventID"],"EventID/I")
    treedict[spillnum].Branch("TrackID",pointers[spillnum, "TrackID"],"TrackID/I")
    # Same thing with the booleans:
    for syst in detsysts.values():
        for det in syst:
            name = 'TrackPresent'+det
            if debug: print "                              ",name
            pointers[spillnum, name] = array( 'B', [ 0 ] ) ## Boolean (unsigned char actually), initialized to 0
            treedict[spillnum].Branch(name,pointers[spillnum, name],name+"/O")
            # and while we're at it, the physical variables, too
            for var in vars:
                if var == "EventID" or var == "TrackID" or var == "SpillID": continue  # Only one unique combo of these, already done above.
                name = var+det
                if debug: print "                              ",name
                pointers[spillnum, name] = array( 'd', [ 0 ] ) ## Double, initialized to 0
                treedict[spillnum].Branch(name,pointers[spillnum, name],name+"/d") ## Add this branch to the treedict
#End of AddTree definition

# Be sure we're working in the correct output file.
outfile.cd()

# Initialize some counters
trackcount = 0
lastspill = 0
spillcount = 0
entrytally = 0
timeindex = {}

if maxspill <= 0: print ("Making as many spills of %s events each as %s contains." % (spillsize,infilename))
else: print ("Processing %s into %s spills of %s events each." % (infilename, maxspill, spillsize ))

# Loop over input TTree
for ds_track in INtuples[starterTree]:
    trackcount += 1
    entrytally += 1
    ## Unique identifiers for this track:
    (event, track) = (int(ds_track.EventID), int(ds_track.TrackID))
    spill = 1 + (event/spillsize) # Arbitrarily group tracks by EventID into spills. 

    # Skip low-E photons
    if ds_track.PDGid == 22:
        E = pow( pow(ds_track.Px,2) + pow(ds_track.Py,2) + pow(ds_track.Pz,2), 0.5)
        if E < gammacutoff: 
            if debug: print "Gamma cutoff at ",gammacutoff
            continue 
    
    # New spill?  Needs a new TTree
    if spill not in newTrees.keys(): 
        AddTree(newTrees, spill, pointers, timeindex) # Also adds a dictionary of TTree entry times for this spill
        entrytally = 1
        spillcount += 1
        if maxspill > 0 and spillcount > maxspill: break
    # Development purposes: (Abbreviated run)
    #if trackcount > 1000: break

    ## What I hate about this is that the newTree looks completely uninvolved.
    ## Of course, it's tied in above, at newTree.Branch(name,pointers[spillnum, name],name+"/f")
    pointers[spill, "SpillID"][0] = spill 
    pointers[spill, "EventID"][0] = event
    pointers[spill, "TrackID"][0] = track

    # Loop over each input TTree, one for each detector in the simulation.
    for tuplename,tuple in INtuples.iteritems():
        for syst in detsysts.values():
            if tuplename not in syst: continue
            # Look for an entry from this event & track
            EntryNum = tuple.GetEntryWithIndex(event,track)
    
            if EntryNum == -1: #Track Not Found
                pointers[spill, 'TrackPresent'+tuplename][0] = False
            else:
                pointers[spill, 'TrackPresent'+tuplename][0] = True
            for var in vars:
                if var == 'EventID' or var == 'TrackID': continue
                ## Pointers are named like 'PzWC4' etc..
                vardet = var + tuplename
                if EntryNum == -1: #Track Not Found
                    pointers[spill, vardet][0] = -123456789  ##Dummy Value
                else: # Track found.
                    for systname,syst in detsysts.iteritems():
                        if tuplename in syst:
                            if debug: print "Filling ",vardet," in ",systname
                            if var == 't':  # Convert times to seconds & add an offset mimicking spill time profile
                                random.seed(event) #Same random offset for each event; progeny of the same pion on target have the same offset.
                                spilltimeoffset = spillinterval * float (spill)
                                newtime = getattr(structs[systname],vardet)*1e-9 + RandomOffsetSeconds() + spilltimeoffset
                                pointers[spill, vardet][0] = newtime
                                # Special need to track entrie by tTOFds
                                if tuplename == 'StartLine': 
                                    # Make a new list if this time has never been seen before
                                    if newtime not in timeindex[spill].keys():  timeindex[spill][newtime] = []
                                    # Add this entry number to the list of entries at this time. (For next stage of processing.)
                                    timeindex[spill][newtime].append(entrytally)
                            else:
                                pointers[spill, vardet][0] = getattr(structs[systname],vardet)
                            #-close else from var == 't'
                        #-close if tuplename in syst
                    #-close for systname,syst in detsysts.iteritems()
                #-close else from if EntryNum == -1
            #-close for var in vars
        #-Closes for syst in detsysts.values()
    #-Closes for tuplename,tuple in INtuples.iteritems()

    # Push the values to the TTree
    newTrees[spill].Fill()
    
#-Closes for ds_track in INtuples[starterTree]

print trackcount, "  total tracks in ",starterTree
print spillcount-1, " total spills (sub-runs)"
outfile.cd()

for newTree in newTrees.values(): newTree.Write()
outfile.Close()
#-------------- /Final Beam --------------------#

infile.Close()

print "Here ya go: ",
commands.getoutput("ls -ltra")

with open(timeindexf, 'wb') as handle:
  pickle.dump(timeindex, handle)
