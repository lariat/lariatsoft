sed -i "s/$1/$2/g" Script.sh
sed -i "s;path;$3;g" Script.sh
sed -i 's/Amps.in.root/Amps.root/g' Script.sh
sed -i 's/Amps.in.pickle/Amps.pickle/g' Script.sh
sed -i 's/Amps.in$SUBSPILL/Amps$SUBSPILL/g' Script.sh
sed -i 's/subspillcountn/jobsperspill/g' Script.sh
jobsub_submit.py -G lariat --memory=4GB --expected-lifetime=25h -N X --OS=SL5,SL6 --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC file://$PWD/Script.sh
sed -i 's/Amps.root/Amps.in.root/g' Script.sh
sed -i 's/Amps.pickle/Amps.in.pickle/g' Script.sh
sed -i 's/Amps$SUBSPILL/Amps.in$SUBSPILL/g' Script.sh
sed -i "s;$3;path;g" Script.sh
sed -i "s/$2/$1/g" Script.sh
