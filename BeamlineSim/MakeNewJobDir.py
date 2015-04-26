#!/usr/bin/env python
# Author: Jason St. John
#
# ./MakeNewJobDir.py <<olddir>> [<<stringtoreplace>>] <<namemod>>
#
# Makes a copy of olddir with a single, uniform change to all the names.
#      Appends to the existing file names, or if stringtoreplace is present,
#      replaces that substring instead of appending.
#
# - Copy over the g4bl input script, renaming it to match.
# - Copy over the condor script which runs the g4bl script, rename it,
#   and modify it to point to the new g4bl script.
# - Copy over the RUN script, make it executable, and modify it to
#   point to the new condor script.
# Whew!

import commands
import optparse
import re

## Get the command line options ##
parser = optparse.OptionParser("usage: %prog [options] olddir [stringtoreplace] newstring\n")
parser.add_option("-v", action="store_true", dest="verbose", default=False,
                  help="Print debugging output.")
parser.add_option ('-B', dest='signedMagField', default="", type="string",
                   help="{+/-}0.00 format Bfield in Tesla.")
parser.add_option ('-E', dest='beamNrg', default="", type="string",
                   help="2-digit beam energy in GeV")
parser.add_option ('--jobcount', dest='jobcount', default=10, type="int",
                   help="Count of jobs to run")
parser.add_option ('--jobsize', dest='jobsize', default=30, type="int",
                   help="Beam events per job")

options, args  = parser.parse_args()
verbose        = options.verbose
signedMagField = options.signedMagField
beamNrg        = options.beamNrg
jobcount     = options.jobcount
jobsize      = options.jobsize
olddir = args[0]

##__Use the number of arguments to decide if a replacement string has been given.__##
replacing = False
namemod = ""
if len(args) < 3:
    namemod = ""
else:
    replaceme = args[1]
    namemod = args[2]
    replacing = True

###################################################################
## Shorten the job size with SI prefixes                         ##
## ...and add it to the dir/ and file names.                     ##
SIprefixes = {}
SIprefixes[0] = ""
SIprefixes[1] = "k"
SIprefixes[2] = "M"
SIprefixes[3] = "G"
thoupow = str(jobsize).count("000")
shortnum = str(jobsize).replace("000","")+SIprefixes[thoupow]
namemod = namemod+"_"+str(jobcount)+"jobsof"+shortnum

###################################################################
## Add the beam and megnet parameters to the name modification   ##
if beamNrg != "":
    str_nrg = beamNrg+"GeV"
    namemod = namemod +"_"+ str_nrg
magsign = ""
if signedMagField != "":
    if float(signedMagField) >= 0: magsign = "pos"
    else: magsign = "neg"
    unsignedMagField = signedMagField.replace("+","").replace("-","")
    str_mag = magsign+unsignedMagField+"Tesla"
    namemod = namemod +"_"+ str_mag

print "                  namemod will be ",namemod

##############################
##__Make the new directory__##
newdirname = olddir+namemod
if replacing:
    newdirname = olddir.replace(replaceme,namemod)
print "newdirname:",newdirname,":"
commands.getoutput('mkdir '+newdirname)
if verbose: print "mkdir "+newdirname
commands.getoutput('cp '+olddir+'/MergeTrees.py '+newdirname+'/.')

################################################################
## Make the new g4bl script. Have to modify this on your own. ##
lsd = commands.getoutput('ls '+olddir+'/LAriaT_*.in')
oldlar = lsd.split("/")[1]
if replacing:
    newlar = oldlar.replace(replaceme,namemod)
else:
    newlar = oldlar.replace(".in",namemod+".in")
sedlar = "sed 's/\.root/"+namemod+"\.root/g"

if beamNrg != "":
    momMeV = str(float(beamNrg)*1000.0)
    sedlar = sedlar + ";s/meanMomentum=8000.0/meanMomentum="+momMeV+"/g"
if signedMagField != "":
    Bfrac = str(float(signedMagField)/0.35)
    sedlar = sedlar + ";s/Bscale=1.0/Bscale="+Bfrac+"/g;"

sedlar +="' "+olddir+"/"+oldlar+" > "+newdirname+"/"+newlar
if verbose: print "sedlar:\n",sedlar
commands.getoutput(sedlar)


#############################################################################
## Make the new condor script, and sed it to point to the new g4bl script. ##
c_script = commands.getoutput('ls '+olddir+'/condor_LAriaT_*.sh')
oldscr = c_script.split("/")[1]
if replacing:
    newscr = oldscr.replace(replaceme,namemod)
else:
    newscr = oldscr.replace(".sh",namemod+".sh")
sedscr = "sed 's/"+oldlar+"/"+newlar+"/g' "+olddir+"/"+oldscr+" > "+newdirname+"/"+newscr
if verbose: print "sedscr:\n",sedscr
sedout = commands.getoutput(sedscr)

###################################################################
## Make the new RUN file, and point it to the new condor script. ##
## Also replace the .root with <<namemod>>.root                  ##
##
## njobs=10
## partperjob=30
sedrun = "sed 's/"+oldscr+"/"+newscr+"/g;s/\.root/"+namemod+"\.root/g"
sedrun = sedrun + ";s/njobs=10/njobs="+str(jobcount)+"/g"
sedrun = sedrun + ";s/partperjob=30/partperjob="+str(jobsize)+"/g"
sedrun = sedrun + ";s/sim30k/sim"+shortnum+"/g"
sedrun = sedrun + ";s/condor_TheJobName.sh/"+newscr+"/g"
sedrun = sedrun + ";s/"+oldlar+"/"+newlar+"/g"
sedrun = sedrun +"' "+olddir+"/RUN > "+newdirname+"/RUN"
if verbose: print "sedrun:\n",sedrun
commands.getoutput(sedrun)
## make it executable.
hey = commands.getoutput("chmod u+x "+newdirname+"/RUN")
if verbose: print hey

########################
## Show what we made. ##
gotls = commands.getoutput('ls -l '+newdirname)
print gotls


