sed -i "s/$1/$2/g" Script.sh
sed -i "s;path;$3;g" Script.sh
sed -i 's/a.in.root/a.root/g' Script.sh
sed -i 's/a.in$PROCESS/a$PROCESS/g' Script.sh
jobsub_submit.py -G lariat -N X --OS=SL5,SL6 --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC file://$PWD/Script.sh
sed -i 's/a.root/a.in.root/g' Script.sh
sed -i 's/a$PROCESS/a.in$PROCESS/g' Script.sh
sed -i "s;$3;path;g" Script.sh
sed -i "s/$2/$1/g" Script.sh
