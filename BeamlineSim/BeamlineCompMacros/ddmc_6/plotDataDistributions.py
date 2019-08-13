import ROOT
from ROOT import *
import sys,os
import argparse


# This code takes as an argument the file 
# we need to generate metadata for
parser = argparse.ArgumentParser()
parser.add_argument("fname"   , nargs='?', default = '../MomentumManipulation/anaTree_postTOFReco.root', type = str, help="insert fileName")
parser.add_argument("treeName", nargs='?', default = 'anatree/anatree'          , type = str, help="insert treeName")

args = parser.parse_args()
fname     = args.fname   
treeName  = args.treeName


f = ROOT.TFile(fname)
t = f.Get(treeName)


hwcPTot = TH1F("hwcPTot","hwcPTot",200,0,2000)
hwcPxTot= TH1F("hwcPxTot","hwcPzTot",200,-200,0)
hwcPyTot= TH1F("hwcPyTot","hwcPyTot",200,-100,100)
hwcPzTot= TH1F("hwcPzTot","hwcPzTot",200,0,2000)
hwcThetaTot= TH1F("hwcThetaTot","hwcThetaTot",200,0,0.2)
hwcPhiTot  = TH1F("hwcPhiTot","hwcPhiTot",800,-4.,4.)
hwcPhiTotNotFUp  = TH1F("hwcPhiTotNotFUp","hwcPhiTotNotFUp",6400,0.,6.4)
hwcPhiTotNotFUpDecentBinning  = TH1F("hwcPhiTotNotFUpDecentBinning","hwcPhiTotNotFUpDecentBinning",200,2.,4.)

hwcPhiVsTheta  = TH2F("hwcPhiVsTheta","hwcPhiVsTheta",200,0,0.2,200,2.,4.)

hwcThetaVsP  = TH2F("hwcThetaVsP","hwcThetaVsP"   ,200,0,0.2,200,0,2000)
hwcPhiVsP    = TH2F("hwcPhiVsP","hwcPhiVsP",200,2.,4.,200,0,2000)
hwcPhiVsThetaVsP  = TH3F("hwcPhiVsThetaVsP","hwcPhiVsThetaVsP",200,0,0.2,200,2.,4.,200,0,2000)

hWC4XvsY    = TH2F("hWC4XvsY","hWC4XvsY"  , 200 , 20., 40., 140, -7., 7.)
hPvsWC4Y    = TH2F("hPvsWC4Y","hPvsWC4Y"  , 140 , -7.,  7., 200,  0., 2000)
hPvsWC4X    = TH2F("hPvsWC4X","hPvsWC4X"  , 200 , 20., 40., 200,  0., 2000)
hwcXVsYVsP  = TH3F("hwcXVsYVsP","hwcXVsYVsP", 200 , 20., 40., 140, -7., 7., 200 , 0. , 2000.)


hWC4X    = TH1F("hWC4X","hWC4X"  , 200 , 20., 40.)
hWC4Y    = TH1F("hWC4Y","hWC4Y"  , 140 , -7.,  7.)

print "I'm looping on your tree, this might take some minutes"
stupidCounter = 0 
for event in t:
    stupidCounter += 1
    if not stupidCounter % 10000.:
        print "Event: ", stupidCounter 
    if event.wcP[0] > 10:
        hwcPTot.Fill(event.wcP[0])
        hwcPxTot.Fill(event.wcPx[0])
        hwcPyTot.Fill(event.wcPy[0])
        hwcPzTot.Fill(event.wcPz[0])
        hwcThetaTot.Fill(event.wctrk_theta[0])
        hwcPhiTot  .Fill(event.wctrk_phi[0])

        WC4x = event.WC4xPos[0]
        WC4y = event.WC4yPos[0]

        hWC4XvsY.Fill(WC4x,WC4y)
        hPvsWC4Y.Fill(WC4y,event.wcP[0])
        hPvsWC4X.Fill(WC4x,event.wcP[0])
        hwcXVsYVsP.Fill(WC4x,WC4y,event.wcP[0])
        hWC4X.Fill(WC4x)
        hWC4Y.Fill(WC4y)

        if event.wctrk_phi[0]>0:
            hwcPhiTotNotFUp  .Fill(event.wctrk_phi[0]) 
            hwcPhiTotNotFUpDecentBinning  .Fill(event.wctrk_phi[0]) 
            hwcPhiVsTheta.Fill(event.wctrk_theta[0],event.wctrk_phi[0])
            hwcPhiVsP.Fill(event.wctrk_phi[0],event.wcP[0])
            hwcThetaVsP.Fill(event.wctrk_theta[0],event.wcP[0])
            hwcPhiVsThetaVsP.Fill(event.wctrk_theta[0],event.wctrk_phi[0],event.wcP[0])
        else:
            hwcPhiTotNotFUp  .Fill(event.wctrk_phi[0]+3.14159265*2) 
            hwcPhiTotNotFUpDecentBinning  .Fill(event.wctrk_phi[0]+3.14159265*2) 
            hwcPhiVsTheta.Fill(event.wctrk_theta[0],event.wctrk_phi[0]+3.14159265*2)
            hwcPhiVsP.Fill(event.wctrk_phi[0]+3.14159265*2,event.wcP[0])
            hwcThetaVsP.Fill(event.wctrk_theta[0],event.wcP[0])
            hwcPhiVsThetaVsP.Fill(event.wctrk_theta[0],event.wctrk_phi[0]+3.14159265*2,event.wcP[0])

outFile = TFile("Gen"+fname,"RECREATE")
outFile.Add(hwcPTot)
outFile.Add(hwcPxTot)
outFile.Add(hwcPyTot)
outFile.Add(hwcPzTot)
outFile.Add(hwcThetaTot)
outFile.Add(hwcPhiTot)
outFile.Add(hwcPhiTotNotFUp)
outFile.Add(hwcPhiTotNotFUpDecentBinning)
outFile.Add(hwcPhiVsTheta)
outFile.Add(hwcPhiVsThetaVsP)
outFile.Add(hwcPhiVsP)
outFile.Add(hwcThetaVsP)
outFile.Add(hWC4XvsY)  
outFile.Add(hPvsWC4Y) 
outFile.Add(hPvsWC4X)
outFile.Add(hwcXVsYVsP)

outFile.Add(hWC4X)
outFile.Add(hWC4Y)
outFile.Write()




    
