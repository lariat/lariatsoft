source /grid/fermiapp/products/common/etc/setups.sh
source /grid/fermiapp/products/larsoft/setup
export GROUP=lariat
export JOBSUB_GROUP=lariat
export PRODUCTS=/grid/fermiapp/products/lariat/:${PRODUCTS}
setup jobsub_client
setup G4beamline v2_16 -q e6:prof:nu
setup git

setup ifdhc

jobsize=Size
SUBspillcount=subspillcountn
first=$((${PROCESS}*${jobsize}))
last=$(( ${first} + $jobsize - 1 ))
echo "PROCESS is: $PROCESS"
echo "jobsize is: $jobsize" 
echo "first = $first"
echo "last = $last"
SUBSPILL=$((${PROCESS}+1 ))
Particlesperspill=$jobsize * $SUBspillcount
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

ifdh cp MergedAtStartLinesim_input.root /pnfs/lariat/scratch/users/$REALUSER/MCdata/MergedAtStartLinesim_input$SUBSPILL.root
ifdh cp MergedAtStartLinesim_input.pickle /pnfs/lariat/scratch/users/$REALUSER/MCdata/MergedAtStartLinesim_input$SUBSPILL.pickle
ls -lrth
echo $CONDOR_DIR_INPUT

#touch MyFile
#ifdh cp MyFile 
