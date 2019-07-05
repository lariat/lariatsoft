import os
import ROOT
minspill=1
maxspill=1000

# The next two lines define the input and output paths. I've commented them for now to require you to uncomment and change them to your lariat/data area.
#indirbase="/lariat/data/users/gpulliam/BeamSim_May2018/TestingProductionScripts/pos60A/triggers/TriggerTree" #name of input files, up to Spill<number>.root. The script expects Spill<number>.root to complete the file name.
#outdirbase="/lariat/data/users/gpulliam/BeamSim_May2018/TestingProductionScripts/pos60A/tracks/" #output directory for track files
for i in xrange(minspill,maxspill+1):
  print "seds for spill "+str(i)
  sedstring1="sed -i 's;SpillSPILLNUM;Spill"+str(i)+";g' BeamCompHY60.C"
  sedstring2="sed -i 's;SpillSPILLNUM;Spill"+str(i)+";g' BeamCompHY60.h"
  sedstring3="sed -i 's;BeamCompHY60SPILLNUM;BeamCompHY60"+str(i)+";g' BeamCompHY60.C"
  sedstring4="sed -i 's;BeamCompHY60SPILLNUM;BeamCompHY60"+str(i)+";g' BeamCompHY60.h"
  sedstring5="sed -i 's;indirbaseSpillSPILLNUM;"+indirbase+"Spill"+str(i)+".root;g' BeamCompHY60.h"
  sedstring6="sed -i 's;outdirbase;"+outdirbase+";g' BeamCompHY60.C"
  
  #os.system(sedstring1)
  #os.system(sedstring2)
  os.system(sedstring3)
  os.system(sedstring4)
  os.system(sedstring5)
  os.system(sedstring6)
  os.system(sedstring1)
  os.system(sedstring2)
  cpstr1="cp BeamCompHY60.C BeamCompHY60Spill"+str(i)+".C"
  cpstr2="cp BeamCompHY60.h BeamCompHY60Spill"+str(i)+".h"
  os.system(cpstr1)
  os.system(cpstr2)
  sedstring1="sed -i 's;Spill"+str(i)+";SpillSPILLNUM;g' BeamCompHY60.C"
  sedstring2="sed -i 's;Spill"+str(i)+";SpillSPILLNUM;g' BeamCompHY60.h"
  sedstring3="sed -i 's;BeamCompHY60"+str(i)+";BeamCompHY60SPILLNUM;g' BeamCompHY60.C"
  sedstring4="sed -i 's;BeamCompHY60"+str(i)+";BeamCompHY60SPILLNUM;g' BeamCompHY60.h"
  sedstring5="sed -i 's;"+indirbase+"Spill"+str(i)+".root;indirbaseSpillSPILLNUM;g' BeamCompHY60.h"
  sedstring6="sed -i 's;"+outdirbase+";outdirbase;g' BeamCompHY60.C" 

  #os.system(sedstring1)
  #os.system(sedstring2)
  os.system(sedstring3)
  os.system(sedstring4)
  os.system(sedstring5)
  os.system(sedstring6)
  os.system(sedstring2)
  os.system(sedstring1)  
with open("Run.sh","w") as outfile:
  for i in xrange(minspill,maxspill+1):
    filenametorun=indirbase+"Spill"+str(i)+".root"
    if os.path.isfile(filenametorun):
      print "root for spill "+str(i)
      command="root -b -q BeamCompHY60Spill"+str(i)+".C \n"
      outfile.write(command)
  outfile.close()  
  
