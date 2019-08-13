##############################################
###      This doesn't work as I thouth     ###
### cause I hate probability distributions ###
##############################################

import ROOT
from ROOT import *
import sys,os
import math
import random
from array import array
import time
from time import gmtime, strftime
print strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())
import argparse


# This code takes as an argument the file
# we need to generate metadata for
parser = argparse.ArgumentParser()
parser.add_argument("fname"   , nargs='?', default = '../MomentumManipulation/anaTree_postTOFReco.root', type = str, help="insert fileName")
parser.add_argument("mommin", nargs='?', default=0, type=int, help="Minimum momentum to draw from when selecting event kinematics")
parser.add_argument("mommax", nargs='?', default=2000, type=int, help="Maximum momentum to draw from when selecting event kinematics")
args = parser.parse_args()
fname     = args.fname
mommin    = args.mommin
mommax    = args.mommax
f = ROOT.TFile(fname)

# Take data histograms
hwcPTot     = f.Get("hwcPTot"    )
hwcPxTot    = f.Get("hwcPxTot"   )
hwcPyTot    = f.Get("hwcPyTot"   )
hwcPzTot    = f.Get("hwcPzTot"   )
hwcThetaTot = f.Get("hwcThetaTot")
hwcPhiTot   = f.Get("hwcPhiTotNotFUpDecentBinning"  )
hwcPhiVsThetaVsP = f.Get("hwcPhiVsThetaVsP")

# Position Generation
hwcXVsYVsPTot = f.Get("hwcXVsYVsP")
hWC4XTot      = f.Get("hWC4X")
hWC4YTot      = f.Get("hWC4Y")


# Normalize them
hwcPTot.Scale(1./hwcPTot.Integral())
hwcPxTot.Scale(1./hwcPxTot.Integral())
hwcPyTot.Scale(1./hwcPyTot.Integral())
hwcPzTot.Scale(1./hwcPzTot.Integral())
hwcThetaTot.Scale(1./hwcThetaTot.Integral())
hwcPhiTot.Scale(1./hwcPhiTot.Integral())
hwcPhiVsThetaVsP.Scale(1./hwcPhiVsThetaVsP.Integral())
hwcXVsYVsPTot.Scale(1./hwcXVsYVsPTot.Integral())
hWC4XTot.Scale(1./hWC4XTot.Integral())
hWC4YTot.Scale(1./hWC4YTot.Integral())

# Prettyfy them
hwcPTot.SetLineColor(4)
hwcPxTot.SetLineColor(4)
hwcPyTot.SetLineColor(4)
hwcPzTot.SetLineColor(4)
hwcThetaTot.SetLineColor(4)
hwcPhiTot.SetLineColor(4)
hWC4XTot.SetLineColor(4)
hWC4YTot.SetLineColor(4)

# Declare Generated histos
hwcPGen     = TH1F("hwcPGen","hwcPGen",200,0,2000)
hwcPxGen    = TH1F("hwcPxGen","hwcPzGen",200,-200,0)
hwcPyGen    = TH1F("hwcPyGen","hwcPyGen",200,-100,100)
hwcPzGen    = TH1F("hwcPzGen","hwcPzGen",200,0,2000)
hwcThetaGen = TH1F("hwcThetaGen","hwcThetaGen",200,0,0.2)
hwcPhiGen   = TH1F("hwcPhiGen","hwcPhiGen",200,2.,4.)
hwcPhiVsThetaGen  = TH3F("hwcPhiVsThetaGen","hwcPhiVsTheta",200,0,0.2,200,2.,4.,200,0,2000)
hWC4XGen    = TH1F("hWC4XGen","hWC4XGen"  , 200 , 20., 40.)
hWC4YGen    = TH1F("hWC4YGen","hWC4YGen"  , 140 , -7.,  7.)



# Declare outFile and arrays
outFile = TFile("GeneratedEvents"+fname,"RECREATE")
t = TTree( 'momentum', 'tree' )

momentumTot = array( 'f', [ 0 ] )
momentumX   = array( 'f', [ 0 ] )
momentumY   = array( 'f', [ 0 ] )
momentumZ   = array( 'f', [ 0 ] )
thetaA      = array( 'f', [ 0 ] )
phiA        = array( 'f', [ 0 ] )
WC4X        = array( 'f', [ 0 ] )
WC4Y        = array( 'f', [ 0 ] )


t.Branch( 'momentumTot', momentumTot, 'momentumTot/F' )
t.Branch( 'momentumX'  , momentumX  , 'momentumX/F' )
t.Branch( 'momentumY'  , momentumY  , 'momentumY/F' )
t.Branch( 'momentumZ'  , momentumZ  , 'momentumZ/F' )
t.Branch( 'theta'  , thetaA  , 'theta/F' )
t.Branch( 'phi'    , phiA    , 'phi/F' )
t.Branch( 'WC4X'   , WC4X    , 'WC4X/F' )
t.Branch( 'WC4Y'   , WC4Y    , 'WC4Y/F' )




P = 0.
n_tries = 300000000
print mommin,mommax
for i in xrange(n_tries):
    if(i%100000)==0:
      print str(i)+": Percentage: "+str(float(i)/float(n_tries)*100)+" Time: "+strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())
    randP          = random.uniform(mommin,mommax)
    randPhi        = random.uniform(200., 400.)
    randTheta      = random.uniform(0., 106.)
    randPhiThetaP   = random.uniform(0, 0.00013)
    import math
    binPhi   = math.ceil(randPhi-200)
    binTheta = math.ceil(randTheta)
    binP     = math.ceil(randP/10.)
#    print randPhi, binPhi," ::: ",  randTheta,  binTheta

    if (randPhiThetaP > hwcPhiVsThetaVsP.GetBinContent(int(binTheta),int(binPhi),int(binP))):
        continue
    

    P = randP
    phi   = float(randPhi)/100.
    theta = float(randTheta)/1000.

    Pz = P*ROOT.TMath.Cos(theta)
    Py = P*ROOT.TMath.Sin(theta)*ROOT.TMath.Sin(phi)
    Px = P*ROOT.TMath.Sin(theta)*ROOT.TMath.Cos(phi)


    #  X and Y generation
    randX    = random.uniform(21.8, 34.5) 
    randY    = random.uniform(-6.1,  6.7)
    binX     = math.ceil((randX - 20.)*10,)   #####
    binY     = math.ceil((randY +  7.)*10.)   #####
    randXYP  = random.uniform(0, 2e-5)       # NEEDS CHECKING
    attempts=0     
    while (randXYP > hwcXVsYVsPTot.GetBinContent( int(binX), int(binY), int(binP) ) ):
        attempts+=1
        randXYP  = random.uniform(0, 0.001)    # NEEDS CHECKING
        randX    = random.uniform(21.8, 34.5) 
        randY    = random.uniform(-6.1,  6.7)
        binX     = math.ceil((randX - 20.)*10,)   #####
        binY     = math.ceil((randY +  7.)*10.)   #####


    X = randX
    Y = randY


    hwcPGen    .Fill(randP)
    hwcPxGen   .Fill(Px)
    hwcPyGen   .Fill(Py)
    hwcPzGen   .Fill(Pz)
    hwcThetaGen.Fill(theta)
    hwcPhiGen  .Fill(phi)
    hwcPhiVsThetaGen.Fill(theta,phi,randP)
    hWC4XGen.Fill(X)
    hWC4YGen.Fill(Y)


    momentumTot[0] = randP
    momentumX[0]   = Px 
    momentumY[0]   = Py 
    momentumZ[0]   = Pz 
    thetaA[0]      = theta
    phiA[0]        = phi
    WC4X[0]        = X
    WC4Y[0]        = Y
    t.Fill()
    


#outFile.cd()
outFile.Add(hwcPGen    )
outFile.Add(hwcPxGen   )
outFile.Add(hwcPyGen   )
outFile.Add(hwcPzGen   )
outFile.Add(hwcThetaGen)
outFile.Add(hwcPhiGen  )
outFile.Add(hwcPhiVsThetaGen)
outFile.Add(hWC4XGen)
outFile.Add(hWC4YGen)



outFile.Add(hwcPTot    )
outFile.Add(hwcPxTot   )
outFile.Add(hwcPyTot   )
outFile.Add(hwcPzTot   )
outFile.Add(hwcThetaTot)
outFile.Add(hwcPhiTot  )
outFile.Add(hwcPhiVsThetaVsP)
outFile.Add(hWC4XTot)
outFile.Add(hWC4YTot)

outFile.Write()
print strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())
#outFile.Close()

#raw_input()  

 

 


