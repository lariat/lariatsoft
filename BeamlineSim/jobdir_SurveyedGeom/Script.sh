source /grid/fermiapp/products/common/etc/setups.sh
#source /cvmfs/oasis.opensciencegrid.org/fermilab/products/common/etc/setup
source /grid/fermiapp/products/larsoft/setup
export GROUP=lariat
export JOBSUB_GROUP=lariat
export PRODUCTS=/grid/fermiapp/products/lariat/:${PRODUCTS}
setup jobsub_client
setup G4beamline v2_16 -q e6:prof:nu
setup git

setup ifdhc

jobsize=100000
first=${PROCESS}*${jobsize}
last=($PROCESS+1)*${jobsize}-1

echo "PROCESS is: $PROCESS"
echo "jobsize is: $jobsize" 
echo "first = $first"
echo "last = $last"

ifdh cp /lariat/app/users/$USER/newtemp1/lariatsoft/BeamlineSim/input input
ls -lrth
g4bl input first=$first last=$last
ls -lrth
chmod 777 sim_input.root
outstage=/pnfs/lariat/scratch/users/$USER

ifdh cp sim_input.root /pnfs/lariat/scratch/users/$USER/MCdata/sim_input$PROCESS.root
ls -lrth
echo $CONDOR_DIR_INPUT

#echo "Hi" > tmp.txt
#ifdh cp tmp.txt /lariat/data/users/soubasis/0slabs/tmp.txt
