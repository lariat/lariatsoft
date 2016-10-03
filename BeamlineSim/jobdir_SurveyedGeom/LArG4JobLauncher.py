###########################################################################################################################################################################################################
## This is a python script that takes as input a directory containing a set of                                                                                                                           
## text files to use as input to prodtext_lariat.fcl.                                               
##
## For each text file in that directory this script does the following:
##
## 1. Edit LArG4_Example.xml, changing the number of events to match the number of events in the txt file
##
## 2. Edit LArG4_Example.xml, changing the work and output directories to be different so the grid will work
##
## 3. Edit LArG4_Example.xml, changing the inputfile to be the text file to process
##
## 4. Edit prodtext_lariat.fcl, changing the SubRun number for every job to avoid a degeneracy run/subrub/event numbers in the output from multiple jobs
##
##
## 5. With the edited .xml and .fcl, submit a job to the grid to process the text file
##
## 6. ????????????
##
## 7. PROFIT!
##
##
## NOTE:  IT IS UP TO THE USER TO DO SOME PRE-EDITS TO THE XML/FCL FILE.
##
## 1. You must edit the "FclDir" in the XML "entity" section to point to your "JobConfigurations" directory in lariatsoft.  IF NOT: You are going to edit MY version of prodtext_lariat.fcl in my JobConfigurations
##
## 2. Change the base paths for the output and work directories for the job. This script only makes new directories at the end.  IF NOT: You are going to send files to my pnfs area.      
##
## 3  Edit prodtext_lariat.fcl, changing the Magnetic Field Descriptions at the bottom to have a B field consistent with the G4BL run. 
##    Somewhere in the text file name should be the current setting. For now, 100A=.35T, and positive polarity corresponds to a negative B_y (check yourself, with F=q(vxB))
##
## usage: python LArG4JobLauncher.py /folder/where/you/have/saved/text/files/
##
## Author: Greg Pulliam gkpullia@syr.edu

import optparse
import commands
import os


parser=optparse.OptionParser("usage: %prof [options] /directoryinput/ \n")
options, args=parser.parse_args()
inpath=args[0]

# Some dictionaries to fill
filestorun=[]
namesoffiles=[]
numberofeventsinfile=[]
numberoflinesoutfile=[]
filenames=os.listdir(inpath)


# The txt files should have hepevt and be a .txt file.  Also, as not to pick up a concatenated text file (output from hepevtConcatenator.py), required "Merged" as well
for file in filenames:
  if file.count('hepevt')>0 and file.count('.txt')>0 and file.count('MergedAtStartLine')>0:
    filestorun.append(file)

#Find how many events are in each file so we know what number to pass to the grid in the .xml file  
for file in filestorun:
  with open(inpath+file, "r") as infile:
    counter=0
    for line in infile:
      if line.count(".")==0:  #particle lines will have at least 1 "." in it when printing various kinematic variables. We dont want those, leaving only the "header" lines for each event. How many lines that pass will be the number of events in the file.
        counter+=1
  numberofeventsinfile.append(counter)


for iter in range(0, len(filestorun)):
  print "Getting ready to launch job to process "+inpath+filestorun[iter]
  with open("LArG4_Example.xml","r") as infile:
    with open("LArG4.xml","w") as outfile:
      for line in infile:
        if line.count("ENTITY")==1 and line.count("jobnumber")==1:
          newline='<!ENTITY jobnumber "'+str(iter)+'">\n'
          outfile.write(newline)
	  #print newline
        elif line.count("numevents")>0:
	  newline='  <numevents>'+str(numberofeventsinfile[iter])+'</numevents>' 
          outfile.write(newline)
	  #print newline
	elif line.count("inputfile")>0:
	  newline="    <inputfile>"+inpath+filestorun[iter]+"</inputfile>\n"
	  outfile.write(newline)
	# While we're writing the new xml, find the fcldir where your prodtext_lariat.fcl file is so we can edit that. 
	elif line.count("ENTITY fclDir")>0:
	  newline=line.split('"',1)[1].split('">')[0]
	  fcldir=newline
	  outfile.write(line)
	else:
	  #print line
	  outfile.write(line)
    outfile.close()
    
# We now have created a new .xml by editing the example. This new file will have the correct number of events to process and will have changed changed the work and output directories to be different (Else Grid get angry, smash things)
# It would still be up to the user to edit other parts of the example file to point to where you want the different work/out directories to be.  
# Now we get the .fcl file and edit it to change the SubRun 

  with open(fcldir+"prodtext_lariat.fcl","r") as fclfile:
    with open(fcldir+"prodtext_lariat_edit.fcl","w") as outfcl:
      print "In fcl file"
      for line in fclfile:
        if line.count("firstRun")>0:
          #inputline=line.split(":",1)[0]
	  #replacedline=line.replace(line,inputline+': "'+inpath+filestorun[iter]+'"\n')
	  #outfcl.write(replacedline)
	  replacedline=line.split(":")[0]+": "+str(iter+1)+"\n"
	  outfcl.write(replacedline)
        else:
          outfcl.write(line)

# Now we have everything edited and ready to go, submit a grid job for each text file
  os.system("project.py --xml LArG4.xml --stage LArG4 --clean")
  os.system("project.py --xml LArG4.xml --stage LArG4 --submit")

  


 
#listcounter=0
#for file in filestorun:
#  with open(inpath+file, "r") as infile:
#    linecounter=0
#    with open(inpath+file.split(".txt",1)[0]+"stripped.txt","w") as outfile:
#      namesoffiles.append(inpath+file.split(".txt",1)[0]+"stripped.txt")
#      for line in infile:
#        linecounter+=1
#	if line.count(".")==0:
#	  numberofevts=line.split(" ",1)[0]	  
#        newline=line.rstrip()
#	if numberoflinesinfile[listcounter]==linecounter:
#	  outfile.write(newline)
#	if numberoflinesinfile[listcounter]!=linecounter:
#	  outfile.write(newline+"\n")
#    listcounter+=1
#    
#    print file+ ": "+numberofevts	
#  infile.close()
#  outfile.close()
#  
#for file in filestorun:  
#  with open(inpath+file.split(".txt",1)[0]+"stripped.txt","r") as infile:
#    outlinecounter=0
#    for line in infile:
#      outlinecounter+=1
#        #print line
#  numberoflinesoutfile.append(outlinecounter)
##In case the file is already there, delete it and re-write it:
#os.remove(inpath+"LArG4FileList.txt")  
#with open(inpath+"LArG4FileList.txt","w") as output:
#  for file in namesoffiles:
#    output.write(file + '\n')
#output.close()


#Ok, we're ready to submit to the grid. If you use the example .xml file in this directory, this works out of the box. 
#Otherwise, you're going to have to edit this line to have your .xml and stage name. Also, you may need --clean if you already have
#stuff in your outdir 

#os.system("project.py --xml LArG4_Example.xml --stage LArG4 --submit")
