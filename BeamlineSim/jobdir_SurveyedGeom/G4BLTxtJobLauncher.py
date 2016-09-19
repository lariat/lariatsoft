#######################################################################
##                                                                    #
##  This script reads in all G4BL root files from a given directory   #
##  and runs a jobsub_submit script to create a text file             #
##  combining particles together into events to run in LArG4          #
##                                                                    #
##  usage: python G4BLTxtJobLauncher.py <path/to/root/file/directory/ #
##  Greg Pulliam gkpullia@syr.edu                                     #
#######################################################################

import os
import optparse
#Get options from command line

parser = optparse.OptionParser("usage: %prog [options] /directoryinput/ \n")  
options, args = parser.parse_args()
inpath = args[0]

# Use arguments to find files for parsing
print inpath
filestoparse=[]
filenames=os.listdir(inpath)
for file in filenames:
# Is the file one we want to run over? Don't run over an already created text file. 
  if file.count('OnlySpilltree')>0 and file.count('hepevt')==0:
    filestoparse.append(file)
    
# Get spill numbers    
spillnumbers=[]
for file in filestoparse:
  #Split the file name to just get the spill number. Remove everything before (and including) "OnlySpilltree", and ".root"
  spillnum=file.split("OnlySpilltree",1)[1].split(".root",1)[0]
  spillnumbers.append(spillnum)
    
#Create a template for the name of the file before the spill number, and then append with .root after    
templatebefore=filestoparse[0].split(spillnumbers[0],1)[0]#.split(".root",1)[0]
templateafter=".root"
for num in spillnumbers:
  # using the template, loop over all the spill numbers you will see and make strings for file names and file paths.
  filetorun=templatebefore+num+templateafter
  commandstring="jobsub_submit.py -G lariat --OS=SL5,SL6 --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC file://$PWD/Treescript.sh"
  print "Going to run jobsubmit to process: "+inpath+filetorun
  fullinpath=inpath
  fullinpath=fullinpath.replace("/","\/")
  textfilename="hepevt_"+templatebefore+num+"Spill_"+num+"thru"+num+".txt"

##Make Treescript.sh edits to have the root file and output text file formatted. 
  os.system("sed -i 's/INPUTFILEPATH/"+fullinpath+"/g' Treescript.sh") 
  os.system("sed -i 's/INPUTFILE/"+filetorun+"/g' Treescript.sh")
  os.system("sed -i 's/OUTFILE/"+textfilename+"/g' Treescript.sh")
  os.system(commandstring)
  os.system("sed -i 's/"+textfilename+"/OUTFILE/g' Treescript.sh")
  os.system("sed -i 's/"+filetorun+"/INPUTFILE/g' Treescript.sh")
  os.system("sed -i 's/"+fullinpath+"/INPUTFILEPATH/g' Treescript.sh")
  
  #os.system(commandstring)
  #os.system(commandstring)
