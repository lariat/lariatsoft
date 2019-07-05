import os
import ROOT
minspill=1
maxspill=1000
# The next two lines define the input and output paths. I've commented them for now to require you to uncomment and change them to your lariat/data area.
#indirbase="/lariat/data/users/gpulliam/BeamSim_May2018/TestingProductionScripts/neg100A/triggers/TriggerTree" #name of input files, up to Spill<number>.root. The script expects Spill<number>.root to complete the file name.
#outdirbase="/lariat/data/users/gpulliam/BeamSim_May2018/TestingProductionScripts/neg100A/tracks/" #output directory for track files
for i in xrange(minspill,maxspill+1):
  print "seds for spill "+str(i)
  sedstring1="sed -i 's;SpillSPILLNUM;Spill"+str(i)+";g' BeamCompHY100.C"
  sedstring2="sed -i 's;SpillSPILLNUM;Spill"+str(i)+";g' BeamCompHY100.h"
  sedstring3="sed -i 's;BeamCompHY100SPILLNUM;BeamCompHY100"+str(i)+";g' BeamCompHY100.C"
  sedstring4="sed -i 's;BeamCompHY100SPILLNUM;BeamCompHY100"+str(i)+";g' BeamCompHY100.h"
  sedstring5="sed -i 's;indirbaseSpillSPILLNUM;"+indirbase+"Spill"+str(i)+".root;g' BeamCompHY100.h"
  sedstring6="sed -i 's;outdirbase;"+outdirbase+";g' BeamCompHY100.C"
  
  #os.system(sedstring1)
  #os.system(sedstring2)
  os.system(sedstring3)
  os.system(sedstring4)
  os.system(sedstring5)
  os.system(sedstring6)
  os.system(sedstring1)
  os.system(sedstring2)
  cpstr1="cp BeamCompHY100.C BeamCompHY100Spill"+str(i)+".C"
  cpstr2="cp BeamCompHY100.h BeamCompHY100Spill"+str(i)+".h"
  os.system(cpstr1)
  os.system(cpstr2)
  sedstring1="sed -i 's;Spill"+str(i)+";SpillSPILLNUM;g' BeamCompHY100.C"
  sedstring2="sed -i 's;Spill"+str(i)+";SpillSPILLNUM;g' BeamCompHY100.h"
  sedstring3="sed -i 's;BeamCompHY100"+str(i)+";BeamCompHY100SPILLNUM;g' BeamCompHY100.C"
  sedstring4="sed -i 's;BeamCompHY100"+str(i)+";BeamCompHY100SPILLNUM;g' BeamCompHY100.h"
  sedstring5="sed -i 's;"+indirbase+"Spill"+str(i)+".root;indirbaseSpillSPILLNUM;g' BeamCompHY100.h"
  sedstring6="sed -i 's;"+outdirbase+";outdirbase;g' BeamCompHY100.C" 

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
      command="root -b -q BeamCompHY100Spill"+str(i)+".C \n"
      outfile.write(command)
  outfile.close()  
  
