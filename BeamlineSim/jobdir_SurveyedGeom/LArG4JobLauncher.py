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
  if file.count('hepevt')>0 and file.count('.txt')>0 and file.count('Merged')>0 and file.count("stripped")<1:  
  #"Stripped" Files have already been processed to have whitespace removed. We don't want those.
    filestorun.append(file)
    
  if file.count('hepevt')>0 and file.count('.txt')>0 and file.count('Merged')>0 and file.count("stripped")==1: 
  # So if we find a stripped file, delete it. The file will be made again anyway. This bypasses read-write privilege problems.
    os.remove(inpath+file)

#Find how many events are in each file so we know what number to pass to the grid in the .xml file  
print "Here are the text files you're about to process :\n" 
for file in filestorun:
  with open(inpath+file, "r") as infile:
    counter=0
    for line in infile:
      if line.count(" ")==1:
        counter+=1
  numberofeventsinfile.append(counter)

#print the number of events to see if it makes sense
for num in numberofeventsinfile:
  print num
 
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
  
for file in filestorun:  
  with open(inpath+file.split(".txt",1)[0]+"stripped.txt","r") as infile:
    outlinecounter=0
    for line in infile:
      outlinecounter+=1
        #print line
  numberoflinesoutfile.append(outlinecounter)
#In case the file is already there, delete it and re-write it:
os.remove(inpath+"LArG4FileList.txt")  
with open(inpath+"LArG4FileList.txt","w") as output:
  for file in namesoffiles:
    output.write(file + '\n')
output.close()

#Ok, we're ready to submit to the grid. If you use the example .xml file in this directory, this works out of the box. 
#Otherwise, you're going to have to edit this line to have your .xml and stage name. Also, you may need --clean if you already have
#stuff in your outdir 

#os.system("project.py --xml LArG4_Example.xml --stage LArG4 --submit")
