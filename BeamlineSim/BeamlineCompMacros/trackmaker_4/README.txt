Author: Greg Pulliam (gkpullia@syr.edu)

At this stage, you should have a series of root files containing all the triggers from a given spill. These files should be stored on /lariat/data.
The next step is to preform a ROOT macro version of TOF and WCReco over these triggers, as well calculate some WCQuality variables that you can cut on in the next step.

In this directory are two python scripts: MakeManyScripts60.py and 100.py. Each is to be used with the corresponding current setting for your simulation. Magnet polarity does not matter at this
stage, just the magnetude of the current.

Each python script works similarly to the previous stage. For each spill between minspill and maxspill, edits are made to the template scripts, BeamCompHY60/100.C/.h to point to the correct input
file and output directory. A copy of the .C/.h is made, one for each spill to be processed. Finally, a bash script, Run.sh, is made to execute the analysis.

All you need to specify is:
min/max spill, 
indirbase, which is the absolute path to the files, not including Spill<NUMBER>.root. The Script will assume this is how the file name ends and will append accordingly
outdirbase, which is the directory to output files. This directory needs to exist beforehand.

usage:
python MakeManyScripts.py
bash Run.sh

Script Details:

For each spill, there are some variables that get set at the preamble, and can be tuned to your liking:

B and L should already be known.

Magnet_midZ: the z coorindate of a plane to perform the midplane match. Usually the average of the z position of the center of the magnets.

middiff sigmascale and shift: After running with the default parameters (sigma=1, shift=0), you can look at the MidDiffX/Y between sim and data to decide how to smear/shift the sim distributions,
if you want to try to make them agree.
 
TOFdeadtime: How long after a TOF hit to wait before allowing another hit.

WCSmear: How much to smear the position of each WChit. Mimics noise effects.

Magnet/Col boundaries: Used to check whether tracks project to the steel of the magnet/collimator. Have to be hard coded, and is informed from the G4BL geometry you've used.


For each event:
Smear each WChit by WCSmear, and perform WCReco based on how many WCs were hit, using the same methods as in data. Once a track is selected, calculate the momentum, as well as determine whether the
track is made from the same particle (pure track) or not (frankentrack). Also note whether the track is picky/HY and 3/4 pt. Check whether track projection impinges on the magnets as well as
calculate the midplane match variables. Also perform a simplified TOF reco that finds a unique TOF between 10-100ns, subject to the TOFdeadtime.

