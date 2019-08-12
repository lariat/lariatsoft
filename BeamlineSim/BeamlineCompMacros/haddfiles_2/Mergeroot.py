import os

#This section loops through the grid output files to find any that are missing. If there is a file missing, perhaps due to the job failing, the spill that the job would have contributed to is declared "bad" and will not be merged.

#This code assumes your job outputs exist in the directory above where this script is. Change it if necessary. Also, make sure this is right name base. Your files should be base<SOMENUMBER>.root 
base='/pnfs/lariat/scratch/users/gpulliam/MCdata/Neg60A_2121_1143/sim_LAriaT_13degProdxn_10degAna_SurveyedGeomMagColl_10000jobsof35k_64GeV_neg60Amps'

#This is where the output will be stored. This directory needs to exist already.
outdir='/lariat/data/users/gpulliam/BeamSim_May2018/TestingProductionScripts/neg60A/'

jobsperspill=10 #The number of jobs to merge into one spill. Default: 10

spills=1000 # The number of spills to merge. Default: 1000. Set lower for testing.

filesizecut=10 #Number of MB you expect one job output to be. 35k POT jobs have around 10.5-11MB. Flag files that have less memory than this as broken files.
badspill=[]
weirdsize=[]
for i in xrange(1,jobsperspill*spills): #Assumes you had 10k jobs. If you had more, increase this. If you had less, then fine, it'll just loop some unneccesarily.
  rootfile=base+str(i)+".root" 
  if not os.path.isfile(rootfile):
    if(i%jobsperspill==0): 
      brokenspill=(i-i%jobsperspill)/jobsperspill 
      print str(i)+" "+str(brokenspill)
      if brokenspill not in badspill:
        badspill.append(brokenspill)
    if(i%jobsperspill!=0): 
      brokenspill=((i-i%jobsperspill)/jobsperspill+1) 
      print str(i)+" "+str(brokenspill)	
      if brokenspill not in badspill:
        badspill.append(brokenspill)      
print "Here are the "+str(len(badspill))+" that are broken due to missing job outputs"

print badspill

#If you think a particular subspill output exists, but is not viable, you can add it by hand in this list. Perhaps the job finished, but the memory of the output G4BL file is way smaller than the other good files. 
#If you like the memory footprint of all your files, this list should be empty.
spillnum=[]

for i in xrange(1,jobsperspill*spills): #Assumes you had 10k jobs. If you had more, increase this. If you had less, then fine, it'll just loop some unneccesarily.
  rootfile=base+str(i)+".root" 
  if os.path.isfile(rootfile) and os.path.getsize(rootfile)<filesizecut*1E6:
    print str(i)+": "+str(os.path.getsize(rootfile)/1E6)+" MB"
    if(i%jobsperspill==0): 
      brokenspill=(i-i%jobsperspill)/jobsperspill 
      if brokenspill not in spillnum:
        spillnum.append(brokenspill)
    if(i%jobsperspill!=0): 
      brokenspill=((i-i%jobsperspill)/jobsperspill+1) 
      if brokenspill not in spillnum:
        spillnum.append(brokenspill)  
#Using this hand-added list, find which spills should be also declared bad.
print "Spills that seem broken due to size"
for i in spillnum:
  if i not in badspill:
    weirdsize.append(int(i))
weirdsizesort=sorted(weirdsize)
print weirdsizesort

#Now, create the hadded spill files if they aren't in these two bad lists

filelist=os.listdir(os.getcwd())
for i in xrange(1,spills+1):
  if i not in badspill:
    if i not in spillnum:
      filetocheck=outdir+"Spill"+str(i)+".root" 
      if not os.path.isfile(filetocheck):    
        command="hadd -f "+outdir+"Spill"+str(i)+".root " #Hadded files will be stored in the same directory as this script. Change if necessary.
        for j in xrange(1,jobsperspill+1): #This should be from 1 to jobsperspill+1. So this is also a 10 check
          subspill=jobsperspill*(i-1)+j 
          filestr=base+str(subspill)+".root "
          command+=filestr
        os.system(command)
	#print command
	  
  
#Also, remove those existing, yet untrusted files. If you have to re-hadd, the respective spill will get picked up in "badspill" instead of you having to add by hand.
for i in spillnum:
  rmcom="rm "+outdir+"Spill"+str(i)+".root"
  os.system(rmcom)
  #print command
  
  
  
