source /grid/fermiapp/lariat/setup_lariat.sh
export GROUP=lariat
export JOBSUB_GROUP=lariat
export PRODUCTS=/grid/fermiapp/products/lariat/:${PRODUCTS}
setup jobsub_client
setup git
setup ifdhc
setup lariatsoft v06_07_00 -q e10:prof

#jobsize=100000
#first=${PROCESS}*${jobsize}
#last=($PROCESS+1)*${jobsize}-1

#echo "PROCESS is: $PROCESS"
#echo "jobsize is: $jobsize" 
#echo "first = $first"
#echo "last = $last"


ifdh cp HEPPATH/hepevtWriter.py hepfile.py
ifdh cp ROOTPATH/ROOTFILE ROOTFILE  #the / between path and file is included in the inputfilepath string
echo ifdh is done: here are the contents of the directory:
ls -lrth
python hepfile.py ROOTFILE
echo python done: here are the new contents of the directory:
ls -lrth
chmod 777 hepevt*.txt


REALUSER=`basename ${X509_USER_PROXY} .proxy | grep -o -P '(?<=_).*(?=_)'`
echo '$USER: ' $USER
echo '$REALUSER: ' $REALUSER

ifdh cp -D hepevt*.txt /pnfs/lariat/scratch/users/$REALUSER/LArG4Files/redo64GeVPos100A
ls -lrth
echo $CONDOR_DIR_INPUT

#touch MyFile
#ifdh cp MyFile 
