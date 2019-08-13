# Copyright (C) 2017-2017 Elena Gramellini <elena.gramellini@yale.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You can read copy of the GNU General Public License here:
# <http://www.gnu.org/licenses/>.

"""@package docstring
Scope of this python script: convert ROOT TTree to HEPEvt format for LArIAT MC
Author: Elena Gramellini
Creation Date: 2016-12-03 
Version 0 
-----------------------------------------------------------------------
Input: name of the root file, name of the ttree, pdg code for particle to generate
Output: a text file with tabulated HEPEvt
"""   


import ROOT
from ROOT import *
import sys,os
import math
import argparse
import random
#+100A
piparams=[-4.16E-7, 9.08E-4, 0.395]
muparams=[-1.11E-7, 1.5E-4,  0.037]
eparams= [2.83E-7, -6.77E-4, .422]  

#+60A
#piparams=[-2.06E-6,   3.07E-3,   -0.312]
#muparams=[-1.01E-7,   1.55E-4,    0.037]
#eparams= [ 2.16E-6,  -3.22E-3,     1.27]
# This function assigns the mass given the pdg code
def pdg2Mass(pdg) :
    mass = 0;
    pdg = math.fabs(float(pdg))
    if (pdg == 11):
        mass = 0.000511
    if (pdg == 13):
        mass = 0.105658374
    if (pdg == 211):
        mass = 0.139570
    if (pdg == 321):
        mass = 0.493667
    if (pdg == 2212):
        mass = 0.93828

    if (mass):
        return mass
    sys.exit("no valid pdg") 

# This function assigns the energy given the mass and momentum
# Relativistic formula is used 
def EnergyCalc(mass, momentum):
    E = math.sqrt(momentum*momentum + mass*mass);
    return  E  


def ProbofSpecies(pz):
  pdg=0
  piprob=piparams[0]*pow(pz,2)+piparams[1]*pz+piparams[2]
  muprob=muparams[0]*pow(pz,2)+muparams[1]*pz+muparams[2]
  eprob=eparams[0]*pow(pz,2)+eparams[1]*pz+eparams[2]
  
  totalp=piprob+muprob+eprob
  #print "pi: "+str(piprob)+" mu: "+str(muprob)+" e: "+str(eprob)+" Total: "+str(totalp)
  if (piprob<0):
    piprob=0
  if(muprob<0):
    muprob=0
  if(eprob<0):
    eprob=0
  scaledpiprob=piprob/totalp
  scaledeprob=eprob/totalp
  scaledmuprob=muprob/totalp
  
  r=random.random()
  if (r<scaledpiprob):
    pdg =211
  elif (r>scaledpiprob and r<(scaledpiprob+scaledmuprob)):
    pdg=-13
  elif (r>(scaledpiprob+scaledmuprob) and r <(scaledpiprob+scaledmuprob+scaledeprob)):
    pdg=-11
  else:
    print "Total Probability is greater than 1. Uh oh."
    sys.exit("Particle pdg choice by momentum failed")
  return pdg 
    




# This code takes as an argument the file 
# we need to generate metadata for
parser = argparse.ArgumentParser()
parser.add_argument("pdg", help="insert pdg you want to keep")
parser.add_argument("fname"   , nargs='?', default = 'referenceTree.root', type = str, help="insert fileName")
parser.add_argument("treeName", nargs='?', default = 'momentum'          , type = str, help="insert treeName")
parser.add_argument("--NoProb", dest="bProb", action="store_false", help="Give a boolean to use momentum dependent probability for a given species")
parser.add_argument("--Prob", dest="bProb", action="store_true", help="Give a boolean to use momentum dependent probability for a given species")
#fname    = "referenceTree.root"
#treeName = "momentum"

args = parser.parse_args()
pdgwanted       = args.pdg
fname     = args.fname   
treeName  = args.treeName
bProb     = args.bProb
print fname, treeName

f = ROOT.TFile(fname)
t = f.Get(treeName)

hwcPzGen    = TH1F("hwcPzGen","hwcPzGen",200,0,2000)
filename = "LArIATHepEvt_pdg_"+ str(pdgwanted)  +".txt"
outrootname="LArIATHepEvt_pdg_"+ str(pdgwanted)  +".root"
outfile=TFile(outrootname,"RECREATE")
print "Opening the file..."
target = open(filename, 'w')

i = 0
for event in t:
    if i%10000==0:
      print i
    pdgtemp=0
    Px      = float(event.momentumX)/1000.
    Py      = float(event.momentumY)/1000.
    Pz      = float(event.momentumZ)/1000.
    if (bProb):
      if(abs(int(pdgwanted))!=211 and abs(int(pdgwanted))!=11 and abs(int(pdgwanted))!=13):
        print("You've tried to use --Prob but aren't generating pions, muons or electrons. Use --NoProb instead.")
	sys.exit("Particle pdg choice by momentum failed")     
      pdgtemp     = ProbofSpecies(event.momentumZ)
    if not bProb:
      pdgtemp=pdgwanted 
    if (abs(int(pdgtemp))==abs(int(pdgwanted))):       
      Mass    = pdg2Mass(pdgtemp)
      Energy  = EnergyCalc(Mass,float(event.momentumTot)/1000. )
      X       = event.WC4X
      Y       = event.WC4Y
      hwcPzGen.Fill(event.momentumZ)
      line1 = str(i) + " 1"
      target.write(line1)
      target.write("\n")
      line2 = "1 " + str(pdgtemp) + " 0 0 0 0 "+ str(Px) + " "+ str(Py) + " "+ str(Pz) + " "+ str(Energy) + " "+ str(Mass) + " "+ str(X) + " "+ str(Y) + " -100.0 0.0"
      target.write(line2)
      target.write("\n")
      i += 1

outfile.Add(hwcPzGen)
outfile.Write()
    
