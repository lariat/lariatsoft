import os
import ROOT
minspill=1
maxspill=1000

# The next two lines define the input and output paths. I've commented them for now to require you to uncomment and change them to your lariat/data area.
#indirbase="/lariat/data/users/gpulliam/BeamSim_May2018/TestingProductionScripts/pos60A/Spill" 
#outdirbase="/lariat/data/users/gpulliam/BeamSim_May2018/TestingProductionScripts/pos60A/triggers/"
for i in xrange(minspill,maxspill+1):
  print "seds for spill "+str(i)
  sedstring1="sed -i 's;SpillSPILLNUM;Spill"+str(i)+";g' lariat_spill.C"
  sedstring2="sed -i 's;(ShiftSPILLNUM-1);("+str(i)+"-1);g' lariat_spill.C"  
  sedstring3="sed -i 's;lariat_spillSPILLNUM;lariat_spill"+str(i)+";g' lariat_spill.C"
  sedstring4="sed -i 's;spill=RawSPILLNUM;spill="+str(i)+";g' lariat_spill.C"
  sedstring5="sed -i 's;indirbase;"+indirbase+str(i)+".root;g' lariat_spill.C"
  sedstring6="sed -i 's;outdirbase;"+outdirbase+";g' lariat_spill.C"
  os.system(sedstring1)
  os.system(sedstring2)
  os.system(sedstring3)
  os.system(sedstring4)
  os.system(sedstring5)
  os.system(sedstring6)
  cpstr1="cp lariat_spill.C lariat_spill"+str(i)+".C"
  os.system(cpstr1)
  sedstring1="sed -i 's;Spill"+str(i)+";SpillSPILLNUM;g' lariat_spill.C"
  sedstring2="sed -i 's;("+str(i)+"-1);(ShiftSPILLNUM-1);g' lariat_spill.C"
  sedstring3="sed -i 's;lariat_spill"+str(i)+";lariat_spillSPILLNUM;g' lariat_spill.C"
  sedstring4="sed -i 's;spill="+str(i)+";spill=RawSPILLNUM;g' lariat_spill.C"
  sedstring5="sed -i 's;"+indirbase+str(i)+".root;indirbase;g' lariat_spill.C"
  sedstring6="sed -i 's;"+outdirbase+";outdirbase;g' lariat_spill.C"
  #os.system(sedstring1)
  os.system(sedstring2)
  os.system(sedstring3)
  os.system(sedstring4)
  os.system(sedstring5)
  os.system(sedstring6)
  os.system(sedstring1)
with open("Run.sh","w") as outfile:
  for i in xrange(minspill,maxspill+1):
    filenametoskip=outdirbase+"TriggerTreeSpill"+str(i)+".root"
    if not os.path.isfile(filenametoskip):
      print "root for spill "+str(i)
      command="root -l -b -q lariat_spill"+str(i)+".C++ \n"
      outfile.write(command)
  outfile.close()  
  
