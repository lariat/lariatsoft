import sys
import optparse
import commands
import os
import pickle

parser = optparse.OptionParser("usage: %prog /inputdirectory/ \n")
parser.add_option ('-v', dest='debug', action="store_true", default=False,
                   help="Turn on verbose debugging.")
parser.add_option ('-f', dest='force', action="store_true", default=False,
                   help="Force overwrite of pickle and root file.")
parser.add_option ('-o', dest='outputfilename', default="", type='string',
                   help="The path and output file name for the .root and .pickle.")
parser.add_option ('-n', dest='jobsperspill', default=-1, type='int',
                   help="Number of jobs needed to make a spill")		   	  
		   
options, args = parser.parse_args()
inpath= args[0]
debug=options.debug
force=options.force
outputfilename=options.outputfilename
jobsperspill=options.jobsperspill

if outputfilename=="":
  exit("No output file name specified")

if jobsperspill==-1:
  exit("Didn't specify number of jobs per spill. Use -n")
  
    
outputfilebase=outputfilename.rsplit(".",1)[0]
outputrootfilename=outputfilebase+".root"
outputpicklefilename=outputfilebase+".pickle"
# Don't overwrite pickle/root file if not supposed to and exit.
if not force and os.path.isfile(outputrootfilename):
  exit("Will not overwrite root file {}. Do you mean -f? ".format(outputrootfilename))
  
if not force and os.path.isfile(outputpicklefilename):
  exit("Will not overwrite pickle file {} Do you mean -f?".format(outputpicklefilename))

SpilltoProcessDict={}


MergeDict = {} #The final merged dictionary. Can contain multiple spills.
#output=open(outputfilename,"rw")
# Because we will add the last "/" later, make sure it's removed if in the input directory
if len(inpath.split("/")[-1])==0:
  inpath=inpath.rsplit("/",1)[0]
print inpath
filelist=os.listdir(inpath)
filebasename=""
#Given Jobsperspill, group process numbers together so the root/pickle files get appended into the correct spill file. 
#Mathematically: For spill, S, process, P, and jobsperspill,J, S=1+((P-1)-((P-1)%J))/J. Requires you to know J for your MC sample, else spills will be merged incorrectly  
#This works if you want P and S to start at 1. To start both at 0, remove all 1s from the above equation.
for file in sorted(filelist):
  if file.count("Amps")!=1: continue
  filebasename=file.split("Amps")[0]+"Amps"
  Process=file.split(".")[0].split("Amps")[1]
  Spill=((int(Process)-1)-(int(Process)-1)%jobsperspill)/jobsperspill+1
  print Process, Spill
  if not Spill in SpilltoProcessDict.keys():
    SpilltoProcessDict[Spill]=[]
  if not Process in SpilltoProcessDict[Spill]:
    SpilltoProcessDict[Spill].append(Process)
print SpilltoProcessDict
print filebasename

#
##Merge pickle dictionaries into one file.
#for file in sorted(filelist):
#  if file.count(".pickle")!=1: continue
#  Jobdict=pickle.load(open(inpath+"/"+file,"rb")) #open the individual pickle file for the job
#  for spillnum in Jobdict.keys(): #looping over job dictionary
#    if not spillnum in MergeDict.keys():
#      MergeDict[spillnum]={} #add dictionary for spill if it doesn't already exist
#    for time, entries in Jobdict[spillnum].iteritems():  #loop over the time and entries in the spill dict
#      if not time in MergeDict[spillnum].keys():
#        MergeDict[spillnum][time]=[]  # make a new list if the time is new
#      MergeDict[spillnum][time].append(entries) #add tree entry numbers to the list for this time
#      #print spillnum, len(MergeDict[spillnum])
#    outputpicklefilename=outputfilebase+"Spill"+str(spillnum)+".pickle"
#    print outputpicklefilename
#    with open(outputpicklefilename,"wb") as f:
#      pickle.dump(MergeDict,f)
#    pickle.close()
#
for spill in SpilltoProcessDict.keys():
  print "Starting to merge pickle for Spill: "+str(spill)
  outputpicklefilename=outputfilebase+"Spill"+str(spill)+".pickle"
  MergeDict = {} #The final merged dictionary. Can contain multiple spills.
  if spill not in MergeDict.keys():
    MergeDict[spill]={}
  for process in SpilltoProcessDict[spill]:
    processfilename=inpath+"/"+filebasename+process+".pickle"
    if not os.path.isabs(processfilename):
      exit("Trying to open a pickle file, {}, that does not exist.".format(processfilename))
    Processdict=pickle.load(open(processfilename,"rb")) #open the individual pickle file for the job  
    for time, entries in Processdict[spill].iteritems(): # Process spills start at 1, everything else starts at 0
      if not time in MergeDict[spill].keys():
        MergeDict[spill][time]=[]
	for iEntry in entries:
          MergeDict[spill][time].append(iEntry)
  with open(outputpicklefilename,"wb") as f:
    pickle.dump(MergeDict,f)
  f.close()
  
for spill in SpilltoProcessDict.keys():
  outputrootfilename=outputfilebase+"Spill"+str(spill)+".root"
  if force: 
    commandstr="hadd -f "+outputrootfilename+" "
  else: 
    commandstr="hadd "+outputrootfilename+" "
  for process in SpilltoProcessDict[spill]:
    processfilename=inpath+"/"+filebasename+process+".root "
    if not os.path.isabs(processfilename):
      exit("Trying to hadd a root file, {}, that does not exist.".format(processfilename))
    commandstr=commandstr+processfilename
  print commandstr
  #os.system(commandstr)   
      #print entries
    
  
#    basename=file.split("Amps")[0]+"Amps"
#    print basename
#    numbertomerge=file.split("jobs")[0].split("SurveyedGeom_")[1]
#    print numbertomerge
  
#MergeDir = {}
#for i in range(0,int(numbertomerge)):
#  print inpath+"/"+basename+str(i)+".pickle"
#  Singledir=pickle.load(open(inpath+"/"+basename+str(i)+".pickle","rb"))
#  print Singledir
#  numberofkeys=len(Singledir[i+1].keys())
#  MergeDir.append(Singledir[i+1])
#  print len(MergeDir.keys())
#MergeDir=pickle.load(open(picklefilename,"rb"))
#
#
#MergedAtStartLinesim_LAriaT_13degProdxn_10degAna_SurveyedGeom_10000jobsof15k_64GeV_pos40Amps1000.pickle
#output.close()
#
#
#while 1
#  try:
#    MergeDir = []
#    MergeDir=pickle.load(open(picklefilename,"rb"))
#    
#    MergeDir.append(pickle.load())
#
#    pickle.dump(MergeDir, output)
#    output.close()
