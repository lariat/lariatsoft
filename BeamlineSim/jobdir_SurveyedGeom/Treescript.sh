source /grid/fermiapp/lariat/setup_lariat.sh
export GROUP=lariat
export JOBSUB_GROUP=lariat
export PRODUCTS=/grid/fermiapp/products/lariat/:${PRODUCTS}
setup jobsub_client
setup git
setup ifdhc
setup lariatsoft v06_07_00 -q e10:prof


ifdh cp HEPPATH/hepevtWriter.py hepfile.py
ifdh cp ROOTPATH/ROOTFILE ROOTFILE #the / between path and file is included in the inputfilepath string
ifdh cp ROOTPATH/PICKLEFILE PICKLEFILE  #the / between path and file is included in the inputfilepath string
echo ifdh is done: here are the contents of the directory:
ls -lrth
python hepfile.py ROOTFILE
echo python done: here are the new contents of the directory:
ls -lrth
chmod 777 hepevt*.txt


REALUSER=`basename ${X509_USER_PROXY} .proxy | grep -o -P '(?<=_).*(?=_)'`
echo '$USER: ' $USER
echo '$REALUSER: ' $REALUSER

ifdh cp -D hepevt*.txt /pnfs/lariat/scratch/users/$REALUSER/MCdata/350kSpills/LArG4/txtwithtriggers/
ifdh cp -D TriggersOnly*.txt /pnfs/lariat/scratch/users/$REALUSER/MCdata/350kSpills/LArG4/txtwithtriggers/
ls -lrth
echo $CONDOR_DIR_INPUT

