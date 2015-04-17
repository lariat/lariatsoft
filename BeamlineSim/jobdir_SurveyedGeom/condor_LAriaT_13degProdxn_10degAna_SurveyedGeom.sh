#!/bin/bash
#
# A script to run  framework (ART) jobs on the local cluster/grid.
#
# It takes 4 arguments, the latter two of which are specified. First two are inherited.:
#   cluster    - the condor job cluster number
#   process    - the condor job process number, within the cluster
#   user       - the username of the person who submitted the job
#   submitdir  - the directory from which the job was submitted ( not used at this time).
#
# Outputs:
#  - All output files are created in the grid scratch space.  At the end of the job
#    all files in this directory will be copied to:
#      /lariat/data/users/MYUSERNAME/outstage/
#     USED TO BE /grid/data/lariat/outstage/$user/${cluster}_${process}
#    This includes a copy of the input files.
#
# Notes:
#
# 1) For documentation on using the grid, see
#      http://mu2e.fnal.gov/atwork/computing/fermigrid.shtm
#
umask 0002
verbose=T

# Copy arguments into meaningful names.
cluster=${CLUSTER}
process=${PROCESS}
user=$1
jobsize=$2

echo "Input arguments:"
echo "Cluster:    " $cluster
echo "Process:    " $process
echo "User:       " $user
echo "jobsize:    " $jobsize

ORIGDIR=`pwd`
echo "ORIGDIR will be "$ORIGDIR
TMP=`mktemp -d ${_CONDOR_SCRATCH_DIR:-/var/tmp}/working_dir.XXXXXXXXXX`
TMP=${TMP:-${_CONDOR_SCRATCH_DIR:-/var/tmp}/working_dir.$$}

{ [[ -n "$TMP" ]] && mkdir -p "$TMP"; } || \
  { echo "ERROR: unable to create temporary directory!" 1>&2; exit 1; }
trap "[[ -n \"$TMP\" ]] && { cd ; rm -rf \"$TMP\"; }" 0

cd $TMP
echo "Beginning work in "`pwd`
# Establish environment and run the job.
export GROUP=lariat
export EXPERIMENT=lariat
export EXTRA_PATH=lariat
export HOME=`cd ../;pwd`
echo "Sourcing setup"
source /cvmfs/oasis.opensciencegrid.org/fermilab/products/common/etc/setup
echo "Setting up ifdh"
setup ifdhc

source /grid/fermiapp/products/lariat/setup
setup G4beamline v2_16 -q e6:prof:nu
setup root v5_34_21b -q e6:nu:prof

echo "Bringing over from " ${ORIGDIR} " all this:"
ls -ltra $ORIGDIR/*
ifdh cp -r $ORIGDIR/* .

# Directory in which to put the output files.
outstage=/lariat/data/users/MYUSERNAME/outstage/

# run jobs
cd $TMP

echo ${HOSTNAME} " is the host and we're in " `pwd`
first=$(( ${process} * ${jobsize} ))
last=$(( ${first} + $jobsize - 1 ))

echo "First and last events: " ${first} " " ${last}
echo "Might find the cardfile in "${CONDOR_DIR_INPUT}
ls -ltra $CONDOR_DIR_INPUT
echo "Gonna try to run."
g4bl ${CONDOR_DIR_INPUT}/LAriaT_13degProdxn_10degAna_SurveyedGeom.in first=$first last=$last
chmod u+x ${CONDOR_DIR_INPUT}/MergeTrees.py
${CONDOR_DIR_INPUT}/MergeTrees.py -l *.root
echo "     FINISHED! Now we have:"
ls -ltra

# Make sure the user's output staging area exists.
test -e $outstage || ifdh mkdir $outstage
if [ ! -d $outstage ];then
   echo "${outstage} is not a directory."
   outstage=/lariat/data/users/MYUSERNAME/condor-tmp/
   echo "Changing outstage directory to: " $outstage 
   exit
fi

# Make a directory in the outstage area to hold all files from this job.
ifdh mkdir ${outstage}/${cluster}_${process}

# Copy all files from the working directory to the output staging area.
echo "Copying *.root and *.out to " ${outstage}/${cluster}_${process}

ifdh cp -D *.root ${outstage}/${cluster}_${process}/.

echo "Maybe we should find the output in "${CONDOR_DIR_INPUT}
ls -ltra $CONDOR_DIR_INPUT


# Make sure EXPERIMENTana (under which name the jobs run on the grid)
# writes these as group rw, so you can rm 'em, etc.
chmod -R g+rw $outstage/${cluster}_${process}

exit 0
