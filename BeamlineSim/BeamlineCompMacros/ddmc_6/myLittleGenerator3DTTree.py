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

fname = "simpleGen.root"
f = ROOT.TFile(fname)

# Take data histograms
hwcPTot     = f.Get("hwcPTot"    )
hwcPxTot    = f.Get("hwcPxTot"   )
hwcPyTot    = f.Get("hwcPyTot"   )
hwcPzTot    = f.Get("hwcPzTot"   )
hwcThetaTot = f.Get("hwcThetaTot")
hwcPhiTot   = f.Get("hwcPhiTotNotFUpDecentBinning"  )
hwcPhiVsThetaVsP = f.Get("hwcPhiVsThetaVsP")

# Normalize them
hwcPTot.Scale(1./hwcPTot.Integral())
hwcPxTot.Scale(1./hwcPxTot.Integral())
hwcPyTot.Scale(1./hwcPyTot.Integral())
hwcPzTot.Scale(1./hwcPzTot.Integral())
hwcThetaTot.Scale(1./hwcThetaTot.Integral())
hwcPhiTot.Scale(1./hwcPhiTot.Integral())
hwcPhiVsThetaVsP.Scale(1./hwcPhiVsThetaVsP.Integral())

# Prettyfy them
hwcPTot.SetLineColor(4)
hwcPxTot.SetLineColor(4)
hwcPyTot.SetLineColor(4)
hwcPzTot.SetLineColor(4)
hwcThetaTot.SetLineColor(4)
hwcPhiTot.SetLineColor(4)

# Declare Generated histos
hwcPGen     = TH1F("hwcPGen","hwcPGen",200,0,2000)
hwcPxGen    = TH1F("hwcPxGen","hwcPzGen",200,-200,0)
hwcPyGen    = TH1F("hwcPyGen","hwcPyGen",200,-100,100)
hwcPzGen    = TH1F("hwcPzGen","hwcPzGen",200,0,2000)
hwcThetaGen = TH1F("hwcThetaGen","hwcThetaGen",200,0,0.2)
hwcPhiGen   = TH1F("hwcPhiGen","hwcPhiGen",200,2.,4.)
hwcPhiVsThetaGen  = TH3F("hwcPhiVsThetaGen","hwcPhiVsTheta",200,0,0.2,200,2.,4.,200,0,2000)


# Declare outFile and arrays
outFile = TFile("results.root","RECREATE")
t = TTree( 'momentum', 'tree' )

momentumTot = array( 'f', [ 0 ] )
momentumX   = array( 'f', [ 0 ] )
momentumY   = array( 'f', [ 0 ] )
momentumZ   = array( 'f', [ 0 ] )
thetaA       = array( 'f', [ 0 ] )
phiA         = array( 'f', [ 0 ] )

t.Branch( 'momentumTot', momentumTot, 'momentumTot/F' )
t.Branch( 'momentumX'  , momentumX  , 'momentumX/F' )
t.Branch( 'momentumY'  , momentumY  , 'momentumY/F' )
t.Branch( 'momentumZ'  , momentumZ  , 'momentumZ/F' )
t.Branch( 'theta'  , thetaA  , 'theta/F' )
t.Branch( 'phi'    , phiA    , 'phi/F' )




P = 0.
for i in xrange(200000000):
#for i in xrange(200000):
    randP          = random.uniform(0., 2000.)
    randPhi        = random.uniform(200., 400.)
    randTheta      = random.uniform(0., 106.)
    randPhiThetaP   = random.uniform(0, 0.001)
    import math
    binPhi   = math.ceil(randPhi-200)
    binTheta = math.ceil(randTheta)
    binP     = math.ceil(randP/10.)
#    print randPhi, binPhi," ::: ",  randTheta,  binTheta

    if (randPhiThetaP > hwcPhiVsThetaVsP.GetBinContent(int(binTheta),int(binPhi),int(binP))):
        continue
    
    #hwcPhiVsTheta.Fill(event.wctrk_theta[0],event.wctrk_phi[0]+3.14159265*2)
    P = randP
    phi   = float(randPhi)/100.
    theta = float(randTheta)/1000.

    Pz = P*ROOT.TMath.Cos(theta)
    Py = P*ROOT.TMath.Sin(theta)*ROOT.TMath.Sin(phi)
    Px = P*ROOT.TMath.Sin(theta)*ROOT.TMath.Cos(phi)



    hwcPGen    .Fill(randP)
    hwcPxGen   .Fill(Px)
    hwcPyGen   .Fill(Py)
    hwcPzGen   .Fill(Pz)
    hwcThetaGen.Fill(theta)
    hwcPhiGen  .Fill(phi)
    hwcPhiVsThetaGen.Fill(theta,phi,randP)

    momentumTot[0] = randP
    momentumX[0]   = Px 
    momentumY[0]   = Py 
    momentumZ[0]   = Pz 
    thetaA[0]       = theta
    phiA[0]         = phi
    t.Fill()
    


#outFile.cd()
outFile.Add(hwcPGen    )
outFile.Add(hwcPxGen   )
outFile.Add(hwcPyGen   )
outFile.Add(hwcPzGen   )
outFile.Add(hwcThetaGen)
outFile.Add(hwcPhiGen  )
outFile.Add(hwcPhiVsThetaGen)

outFile.Add(hwcPTot    )
outFile.Add(hwcPxTot   )
outFile.Add(hwcPyTot   )
outFile.Add(hwcPzTot   )
outFile.Add(hwcThetaTot)
outFile.Add(hwcPhiTot  )
outFile.Add(hwcPhiVsThetaVsP)
outFile.Write()
#outFile.Close()

#raw_input()  

 

 


