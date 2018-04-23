#!/bin/bash

source  /cvmfs/mu2e.opensciencegrid.org/artexternals/setup
source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setup
source /cvmfs/lariat.opensciencegrid.org/setup_lariat.sh
export GROUP=lariat
export JOBSUB_GROUP=lariat
setup jobsub_client v1_2_6_2
setup G4beamline v2_16a -q e6:prof:nu
setup ifdhc v2_3_0 
setup root v5_34_23 -q e6:prof
jobsize=1000
SUBspillcount=1
first=$((${PROCESS}*${jobsize}))
last=$(( ${first} + $jobsize - 1 ))
echo "PROCESS is: $PROCESS"
echo "jobsize is: $jobsize" 
echo "first = $first"
echo "last = $last"
SUBSPILL=$((${PROCESS}+1 ))
Particlesperspill=$((${jobsize} * ${SUBspillcount}))
ifdh cp path/input input
ifdh cp path/MergeTrees.py MergeTrees.py
ls -lrth
g4bl input first=$first last=$last
ls -lrth

chmod 777 sim_input.root
chmod 777 MergeTrees.py
./MergeTrees.py sim_input.root --subspillnumber $SUBSPILL --subspillcount $SUBspillcount --spillsize $Particlesperspill
chmod 777 MergedAtStartLinesim_input.root
chmod 777 MergedAtStartLinesim_input.pickle
ls -lrth

REALUSER=`basename ${X509_USER_PROXY} .proxy | grep -o -P '(?<=_).*(?=_)'`
echo '$USER: ' $USER
echo '$REALUSER: ' $REALUSER

ifdh cp MergedAtStartLinesim_input.root /pnfs/lariat/scratch/users/$REALUSER/MCdatatest/MergedAtStartLinesim_input$SUBSPILL.root
ifdh cp MergedAtStartLinesim_input.pickle /pnfs/lariat/scratch/users/$REALUSER/MCdatatest/MergedAtStartLinesim_input$SUBSPILL.pickle
ls -lrth
echo $CONDOR_DIR_INPUT

#touch MyFile
#ifdh cp MyFile 
