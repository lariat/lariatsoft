#!/usr/bin/env python
# Author: Jason St. John
#
# Read in a file of one or more single-spill TTrees representing particle interactions in G4BL, all assumed to be part of
# the same simulated particle interactions such that EventID and TrackID combinations are unique, and that times have been set
# to reflect beam time structure. 
# 
# First Loop:  Fill a dictionary of lists of single-particle entries, keyed by their tStartLine.
# Second Loop: Find the triggering particles, simulating trigger system deadtime.
# Third Loop:  Print each triggering particle to the hepevt file, along with all particles in a time window around it. 
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
import copy
import pickle

parser = optparse.OptionParser("usage: %prog [options]<input file.ROOT> \n")
parser.add_option ('--maxspill', dest='maxspill', type='int',
                   default = -1,
                   help="Abbreviate processing to this many spill trees.")
parser.add_option ('--firstspill', dest='firstspill', type='int',
                   default = -1,
                   help="Spill number to start processing.")
parser.add_option ('--lastspill', dest='lastspill', type='int',
                   default = -1,
                   help="Spill number to start processing.")
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
parser.add_option ('--photoncutoff', dest='photoncutoff', type='float',
                   default =0.0005,
                   help="Minimum photon inclusion energy for textfile (GeV)")
parser.add_option ('--timeindex', dest='picklefilename', type='string',
                   default ='',
                   help="Name of pickle file containing time index.")
parser.add_option ('--gimmeTrigger', dest='gimmeTrigger',action="store_true", default=False,
                   help="Also return a text file of just the triggering particles.")
options, args = parser.parse_args()
maxspill       = options.maxspill
firstspill     = options.firstspill
lastspill      = options.lastspill
spillinterval  = options.spillinterval
driftinterval  = options.driftinterval
debug          = options.debug
test           = options.test
photoncutoff   = options.photoncutoff
picklefilename = options.picklefilename
gimmeTrigger   = options.gimmeTrigger
infile         = args[0]
os.system("date")
triggeronlytriggercount=0
## Attempt to make picklefilename automatically
if picklefilename == '': picklefilename=str(infile.split(".root")[0])+".pickle"
    
if os.path.isfile(picklefilename): timeindexfromfile = True
else: 
    timeindexfromfile = False
    print "Unable to find pickle file ",picklefilename,".  Will extract time index from ROOT file."

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
# Given a pile of leaves and the chosen t_0, 
# return the hepevt line: (GeV, ns, cm)
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
    t = ( pile.tStartLine.GetValue() - tzero )* 1.0e9        #9 NTS: Shift zero to trigger time.
    hepstr = '1 {0:d} 0 0 0 0 {1} {2} {3} {4} {5} {6} {7} {8} {9}\n'.format(pdg, Px, Py, Pz, E, m, x, y, z, t)
    if E<photoncutoff and pdg==22: return "no"
    return hepstr

# Given a pointer to an entry in the SpillTree, return 
# a boolean as to whether the trigger condition was met.
def triggercondition(pile):
    if debug:
        if pile.TrackPresentDet1.GetValue(): print 'Det1'
        if pile.TrackPresentDet2.GetValue(): print 'Det2'
        if pile.TrackPresentDet3.GetValue(): print 'Det3'
        if pile.TrackPresentDet4.GetValue(): print 'Det4'
        if pile.TrackPresentTOFus.GetValue(): print 'TOFus'
        if pile.TrackPresentTOFds.GetValue(): print 'TOFds' #Horz removed
        
    # Wire Chambers
    return (pile.TrackPresentDet1.GetValue() and 
            (pile.TrackPresentDet2.GetValue() or pile.TrackPresentDet3.GetValue()) and
            pile.TrackPresentDet4.GetValue() and 
            # Time of Flight, too
            pile.TrackPresentTOFus.GetValue() and pile.TrackPresentTOFds.GetValue() )

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
        spill = int(name.split("_Spill")[1])
        InputSpillTrees[spill] =  key.ReadObj() # Store pointer to this TTree in the dictionary
        spillcount += 1
        if maxspill > 0 and spillcount > maxspill: break

# To speed up finding tracks by EventID and TrackID, 
# build an index with these as the major and minor indices.
for name, tuple in InputSpillTrees.iteritems():
    print "Building an index for ",name,"...",
    tuple.BuildIndex("tStartLine")
    print "Done."

#### The Output: A file with a tree in it. ####
infilepath = os.path.abspath(infilename)
infilename = os.path.basename(infilename)

outfilename = "hepevt_"+infilename.replace('.root','.txt')
if gimmeTrigger:
  Triggeroutfile="TriggersOnly"+infilename.replace('.root','.txt')
  triggerfile=open(Triggeroutfile,'w')
outfile = open(outfilename,'w')

# Initialize a handy list of spill numbers visited
spillnums = []
die = False
# Loop over input TTree objects
for spill, intree in InputSpillTrees.iteritems():
    # Extract the spill number from the TTree's name. Comes after the _Spill
    if firstspill > 0 and spill < firstspill: continue
    if lastspill > 0 and spill > lastspill: continue

    n_entries = intree.GetEntriesFast()
    print "Starting ",intree.GetName()," with ",n_entries," entries."
    spillnums.append(spill)

    # Gonna need a handy little class hooked to this tree:
    # __Pile (pyl) of leaves representing a single particle__
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
    triggertimes = []
    triggerentrynums = []
    # Get the time index allentriesbytime from a pickle file if possible. 
    if timeindexfromfile:
        # Get a dictionary of index values, with their times as the keys
        allentriesbytime= pickle.load(open(picklefilename,"rb"))
        allentriesbytime=allentriesbytime[spill]
        entrytimes=allentriesbytime.keys()
        # Also some useful other lists:
        triggertimes = []
        triggerentrynums = []
    else: # Have to get it the old-fashioned way
        entrytimes=[]
	allentriesbytime={}
        # First Loop over this tree: Get the entry numbers and tStartLine (if defined)
        if debug: print '    Beginning 1st loop over', n_entries," entries."
        for n in xrange(0, n_entries):
            intree.GetEntry(n) # Fill pyl with values from the entry at index n
            # Has to occur in StartLine
            timeExists = pyl.TrackPresentStartLine.GetValue()
            if not timeExists: continue
            # Get the time value
            time = float(pyl.tStartLine.GetValue())
            if debug: print n,":",time
            entrytimes.append(time) # All tStartLine values. Values can be non-unique.
            # Make sure there's a list of entry numbers in the dictionary for this time
            if time not in allentriesbytime.keys():
	      allentriesbytime[time] = []
            # For each unique tStartLine, make a list of the entry numbers.
            allentriesbytime[time].append(n)

    # What did we get?  Any non-unique index values?
    timecount = len(entrytimes)
    uniqcount = len(set(entrytimes))
    if debug: print "Time count:",len(allentriesbytime.keys())," from a spill of size:", n_entries
    if debug: print "With entry.TrackPresentStartLine:",intree.GetEntries("TrackPresentStartLine")
    if debug: print "total time count: ",len(entrytimes)," (unique:",len(set(entrytimes)),", non-unique:",100.*(float(timecount)-float(uniqcount))/float(timecount),"%)"

    # Second Loop over the tree:
    if debug: print '    Beginning 2nd loop over all', len(allentriesbytime.keys()),"distinct particle times."
    # Visit entries in order by their times, and check for triggering particles.
    # Must be separate from above because triggers must be known to be >2 driftintervals after foregoing triggers.
    LastTriggerTime = float(-1)
    firsttime = sorted(allentriesbytime.keys())[0] # Keep for reference.

    # Loop over all the tStartLine values IN INCREASING ORDER.
    for time in sorted(allentriesbytime.keys()):

        #Get the delta_t and update LastTriggerTime
        delta_t = time-LastTriggerTime

        # For this time value, visit one or more entries n:
        for n in allentriesbytime[time]:
            ret = intree.GetEntry(n)
            if debug: print "time: ",time,"   delta_t:",delta_t,"   n:",n            
            if ret == -1: die('No entry number '+str(n))

            # Was this a triggering particle?
            pdg = pyl.PDGidStartLine.GetValue()
            if pdg in neutrals: continue # Photons? Don't care.

            # For the first trigger time, or of it's a later trigger and it has been long enough to trigger again:
            timeok = (LastTriggerTime == -1 or (LastTriggerTime > -1 and delta_t > 2.0*driftinterval))
            if timeok and triggercondition(pyl): # Check the physical trigger condition

                #..if so, add it to the list
                triggertimes.append(time) 
                triggerentrynums.append(n) 
                LastTriggerTime = time
                print "Triggering: ",n,":",time,"PDG:  ",pdg
                if gimmeTrigger:
		  triggeronlytriggertime=pyl.tTOFds.GetValue()
		  triggeronlytxtstr=gimmestr(pyl,triggeronlytriggertime)
		  triggerfile.write(str(triggeronlytriggercount)+" 1\n")
		  triggeronlytriggercount+=1
		  triggerfile.write(triggeronlytxtstr)
        # Testing mode. Had enough? Break the loop.
        if test and time - firsttime > 0.1: 
            die = True
            break

    # Third and final loop: 
    if debug: print "    Beginning 3rd loop over", len(triggertimes), "triggers."
    #Collect the hepevt fields for all particles which might give signals in the triggered events.
    eventnum = 0 # Unique within the spill (only).
    for time in triggertimes:
        particlelines = {} # A dictionary where we can keep the lines for this hepevt
        particlecount = 0
        delta_t = 0.       # 
        offset = 0
        print "\n Collecting tracks around trigger at ",time
        # OK, where are we in the sorted list of particle times?
        triggerTimeIndex = sorted(allentriesbytime.keys()).index(time)
        checkindex = triggerTimeIndex - offset
        entry_t = sorted(allentriesbytime.keys())[checkindex]
        delta_t = abs(time - entry_t)
        if debug: print "time:",time,"  offset: ",offset,"entry_t:",entry_t,"    delta_t:",delta_t

        # To mimic trigger decision latency, grab the tTOFds of the triggering particle
        # This will be subtracted off the tStartLine of all particles in the event window.
        tTriggers = []
        if debug: print 'Grabbing TOFds values among',len(allentriesbytime[time]),'values:',
        for n in allentriesbytime[time]: # In case of multiple trigger-time particles
	    if debug: print "Entry",n,
            ret = intree.GetEntry(n)
            if ret == -1: exit ('No entry #'+str(n))
            # Only interested in particles which might have triggered
            # (and not merely coincident with the triggering particle)
            if triggercondition(pyl): 
                # Trigger t_zero chosen at DS TOF paddle:
                tTriggers.append(pyl.tTOFds.GetValue()) 
            print "\n"

        # Take the earliest tTOFds of any trigger particle candidates at the trigger time.
        # (Nearly always there is one and only one candidate triggering particle. Paranoia.)
	if not len(tTriggers)>0: 
            print "tTriggers zero length" 
            continue # Maybe the next trigger will be okay
        tTrigger = sorted(tTriggers)[0]
        if debug: print "    Using tTrigger:",tTrigger

        # Print the hepevt lines: (GeV, ns, cm)
        # 1 [pdg] 0 0 0 0 px py pz E m x y z t

        # At which times were there particles BEFORE this trigger, close enough to put energy in the TPC?
        while (True):
            offset +=1
            checkindex = triggerTimeIndex - offset
            if checkindex < 0: break #Not too early!

            # Get the event time at this nearby index
            entry_t = sorted(allentriesbytime.keys())[checkindex]

            # How far away in time from our triggering particle? 
            delta_t = abs(time - entry_t)
            if delta_t > driftinterval: break # Don't go more than one driftinterval before the triggering particle.

            if debug: print "offset: -",offset,"entry_t:",entry_t,"    delta_t:",delta_t

            # Loop over any entries which had this time value
            for n in allentriesbytime[entry_t]:
                if debug: print "Getting n:",n
                ret = intree.GetEntry(n) # Set pyl's pointers to this entry's leaf values
                if ret == -1: exit ('No entry #'+str(n))
                txtstr = gimmestr(pyl, tTrigger)
		if txtstr=="no": continue
                # Copy the values, so that txtstr can change later without affecting this entry in particlelines
                particlelines[-1 * particlecount] = copy.copy(txtstr) 
                particlecount += 1
                print txtstr,

        # Precisely at the time of the triggering particle:
        for n in allentriesbytime[time]:
            ret = intree.GetEntry(n)
            if ret == -1: exit ('No entry #'+str(n))
            txtstr = gimmestr(pyl, tTrigger)
	    if txtstr=="no": continue
            # Copy the values, so that txtstr can change later without affecting this entry in particlelines
            particlelines[particlecount] = copy.copy(txtstr)
            particlecount += 1
            print txtstr,

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
            if delta_t > driftinterval: break # Don't go more than one driftinterval after the triggering particle.

            if debug: print "offset: +",offset,"entry_t:",entry_t,"    delta_t:",delta_t

            # Loop over any entries which had this time value
            for n in allentriesbytime[entry_t]:
                if debug: print "Getting n:",n
                ret = intree.GetEntry(n) # Set pyl's pointers to this entry's leaf values
                if ret == -1: exit ('No entry #'+str(n))
                txtstr = gimmestr(pyl, tTrigger)
		if txtstr=="no": continue
                # Copy the values, so that txtstr can change later without affecting this entry in particlelines
                particlelines[particlecount] = copy.copy(txtstr)
                particlecount += 1
                print txtstr,

        if debug: print "\nPrint to file:"
        # Print event header, which must list the particle count.
        line = str(eventnum)+' '+str(particlecount)+'\n'
        outfile.write(line) 
        eventnum += 1
        if debug: print line
        # For each line we want to print in this event,
        for num in sorted(particlelines.keys()):
            line = particlelines[num]
            # print the line for this particle
            outfile.write(line)
            if debug: print line

outfile.close()
if gimmeTrigger:
  triggerfile.close()
lo_spill = sorted(spillnums)[0]
hi_spill = sorted(spillnums)[len(spillnums)-1]
spillrange_str = 'Spill_'+str(lo_spill)+'thru'+str(lo_spill)
if not (hi_spill-lo_spill == 0): 
    os.rename(outfilename,outfilename.replace('.txt',spillrange_str+'.txt')) # List both spill numbers if the difference is not 0.
infile.Close()
os.system("date")
