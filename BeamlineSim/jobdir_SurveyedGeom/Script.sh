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
first=${PROCESS}*${jobsize}
last=($PROCESS+1)*${jobsize}-1

echo "PROCESS is: $PROCESS"
echo "jobsize is: $jobsize" 
echo "first = $first"
echo "last = $last"

ifdh cp path/input input
ls -lrth
g4bl input first=$first last=$last
ls -lrth
chmod 777 sim_input.root


REALUSER=`basename ${X509_USER_PROXY} .proxy | grep -o -P '(?<=_).*(?=_)'`
echo '$USER: ' $USER
echo '$REALUSER: ' $REALUSER

ifdh cp sim_input.root /pnfs/lariat/scratch/users/$REALUSER/MCdata/sim_input$PROCESS.root
ls -lrth
echo $CONDOR_DIR_INPUT

#touch MyFile
#ifdh cp MyFile 
