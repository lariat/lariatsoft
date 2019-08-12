source  /cvmfs/mu2e.opensciencegrid.org/artexternals/setup
source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setup
source /cvmfs/lariat.opensciencegrid.org/setup_lariat.sh
export GROUP=lariat
export JOBSUB_GROUP=lariat
setup jobsub_client v1_2_6_2
setup G4beamline v2_16a -q e6:prof:nu
setup ifdhc v2_3_0 
setup root v5_34_23 -q e6:prof
jobsize=35000
SUBspillcount=10
SUBSPILLNUM=$((${PROCESS}%10)) #Processes ending in 0 is the first subspill for its spill. 
Spill=$((((${PROCESS}-${SUBSPILLNUM})/10)+1)) 
#Spill number starts at 1, so Process 0-9 go to spill 1, and so on
first=$((${SUBSPILLNUM}*${jobsize} + ${Spill}))
last=$(( ${first} + $jobsize - 1  ))
echo "PROCESS is: $PROCESS"
echo "SUBSPILLNUM is: $SUBSPILLNUM"
echo "jobsize is: $jobsize" 
echo "first = $first"
echo "last = $last"
echo "Spill = $Spill"
SUBSPILL=$((${PROCESS}+1)) #Later parts of the code starts counting at 1 (Spill1 is subspills 1-10). Purely a naming convention for files to be hadded later.
Particlesperspill=$((${jobsize} * ${SUBspillcount}))
ifdh cp path/input input
ls -lrth
command1="sed -i 's/randomseed Set PROCESSJOBNUMBER/randomseed Set "$PROCESS"/g' input"
echo $command1
eval $command1
g4bl input first=$first last=$last
ls -lrth


#chmod 777 sim_input.root


REALUSER=`basename ${X509_USER_PROXY} .proxy | grep -o -P '(?<=_).*(?=_)'`
echo '$USER: ' $USER
echo '$REALUSER: ' $REALUSER

ifdh cp sim_input.root path/sim_input$SUBSPILL.root
ls -lrth
echo $CONDOR_DIR_INPUT

#touch MyFile
#ifdh cp MyFile 
