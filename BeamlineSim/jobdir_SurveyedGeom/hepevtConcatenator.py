

import optparse
import commands
import os

parser=optparse.OptionParser("usage: %prof [options] /directoryinput/ \n")
options, args=parser.parse_args()
inpath=args[0]

filestocat=[]
filenames=os.listdir(inpath)
# The txt files should have hepevt and be a .txt file.
for file in filenames:
  if file.count('hepevt')>0 and file.count('.txt')>0 and file.count('Merged')>0:
    filestocat.append(file)
#####################################################################################
# Figure out what the current setting was for the G4BL that generated the text file #
# This should be "posXXXA" or "negXXXA" somewhere in the text file name.            #
#####################################################################################

# Split the file name into strings separated by "_". One of these will contain the current settings
stringsinfilename=filestocat[0].split("_")

#loop over these strings and find the one that has the current in it.
for string in stringsinfilename:
  if string.count('neg')>0 or string.count('pos')>0:
    currentstring=string

# With the current setting found, create "hepevt_<NSpills>_<currentsetting>.txt and combine the text files from N spills
# output concatenated textfile will be in the same place as the input text files.
with open(inpath+'hepevt_'+str(len(filestocat))+"Spills_"+currentstring+'.txt', 'w') as outfile:
  print outfile
  for file in filestocat:
    with open(inpath+file) as infile:
      for line in infile:
        outfile.write(line)
	

