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
from datetime import tzinfo, timedelta, datetime
#Get options from command line

parser = optparse.OptionParser("usage: %prog [options] /directoryinput/ \n")  
parser.add_option('--heppath', dest='heppath', type='str', default=os.getcwd(), help="The absolute path to where hepevtWriter.py can be found")
options, args = parser.parse_args()
heppath=options.heppath.rstrip("\n").rstrip("/")
#heppath=options.heppath.rstrip("/")
inpath = args[0].rstrip("/")


## Use arguments to find files for parsing
print inpath

filestoparse=[]
filenames=os.listdir(inpath)
for file in sorted(filenames):
# Is the file one we want to run over? Don't run over an already created text file. 
  if file.count('Spill')>0 and file.count('hepevt')==0 and file.count(".root")==1:
    filestoparse.append(file)
# Get spill numbers    
spillnumbers=[]
for file in sorted(filestoparse):
  #Split the file name to just get the spill number. Remove everything before (and including) "OnlySpilltree", and ".root"
  spillnum=file.split("Spill",1)[1].split(".root",1)[0]
  spillnumbers.append(spillnum)
    
#Create a template for the name of the file before the spill number, and then append with .root after    
templatebefore=filestoparse[0].split(".root",1)[0]
templateafter=".root"
datenow=str(datetime.now().strftime('%Y-%m-%dT%H.%M.%S'))
print datenow
clusternumbersfile="jobs_"+str(datenow)
#os.system("touch "+clusternumbersfile)
print clusternumbersfile
for file in filestoparse:
  filetorun=str(file)
  pickletorun=str(file).split(".root",1)[0]+".pickle"
  # using the template, loop over all the spill numbers you will see and make strings for file names and file paths1.
  commandstring="jobsub_submit.py -G lariat --OS=SL5,SL6 --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC --expected-lifetime=10h file://$PWD/Treescript.sh | grep -e 'cluster'| awk '{print $0}'  >> "+clusternumbersfile
  #print "Going to run jobsubmit to process: "+inpath+filetorun

  print "ROOTPATH: "+inpath
  print "ROOTFILE: "+filetorun
  print "PICKLEFILE: "+pickletorun
  print "HEPPATH: "+heppath
  command1="sed -i 's:ROOTPATH:"+inpath+":g' Treescript.sh"
  command2="sed -i 's:ROOTFILE:"+filetorun+":g' Treescript.sh"
  command3="sed -i 's:PICKLEFILE:"+pickletorun+":g' Treescript.sh"
  command4="sed -i 's:HEPPATH:"+heppath+":g' Treescript.sh"
  command5="sed -i 's:"+heppath+":HEPPATH2:g' Treescript.sh"
  command6="sed -i 's:"+filetorun+":ROOTFILE2:g' Treescript.sh"
  command7="sed -i 's:"+pickletorun+":PICKLEFILE2:g' Treescript.sh"
  command8="sed -i 's:"+inpath+":ROOTPATH2:g' Treescript.sh"
##Make Treescript.sh edits to have the root file and output text file formatted. 
  os.system(command1)
  os.system(command2)
  os.system(command3)
  os.system(command4)
  os.system(commandstring)
  os.system(command5)
  os.system(command6)
  os.system(command7)
  os.system(command8) 

  #os.system(commandstring)
  #os.system(commandstring)
