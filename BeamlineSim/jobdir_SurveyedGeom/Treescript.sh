source /grid/fermiapp/lariat/setup_lariat.sh
export GROUP=lariat
export JOBSUB_GROUP=lariat
export PRODUCTS=/grid/fermiapp/products/lariat/:${PRODUCTS}
setup jobsub_client
setup git
setup ifdhc
setup lariatsoft v05_15_00 -q e9:debug

#jobsize=100000
#first=${PROCESS}*${jobsize}
#last=($PROCESS+1)*${jobsize}-1

#echo "PROCESS is: $PROCESS"
#echo "jobsize is: $jobsize" 
#echo "first = $first"
#echo "last = $last"


ifdh cp INPUTFILEPATH/hepevtWriter.py hepfile.py
ifdh cp INPUTFILEPATHINPUTFILE INPUTFILE  #the / between path and file is included in the inputfilepath string
ls -lrth
python hepfile.py INPUTFILE
ls -lrth
chmod 777 OUTFILE


REALUSER=`basename ${X509_USER_PROXY} .proxy | grep -o -P '(?<=_).*(?=_)'`
echo '$USER: ' $USER
echo '$REALUSER: ' $REALUSER

ifdh cp OUTFILE /pnfs/lariat/scratch/users/$REALUSER/LArG4Files/OUTFILE
ls -lrth
echo $CONDOR_DIR_INPUT

#touch MyFile
#ifdh cp MyFile 
