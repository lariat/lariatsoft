Author: Greg Pulliam gkpullia@syr.edu

At this stage, you should have a series of root files containing all particles created in a given spill, separated into trees, one for each detector: TOFUS, WC1, WC2, WC3, WC4, and DSTOF.
The goal of this stage is to process each spill to find triggers, similar to how it is done in data. For each trigger, the hits in all detectors in the readout window surrounding the
trigger time are stored. Each trigger is the equivalent to a "sliced event" in data.

Each spill must be processed individually. There is a template script, lariat_spill.C, which gets duplicated for each spill, with some edits made for file path names. MakeManyScripts.py
does this script duplication for you, making a lariat_spill.C for each spill to be processed. You'll have many (~1000) lariat_spill.C files when this is done.

The only edits to MakeManyScripts.py you must make is to alter the max and min spill to process (minspill, maxspill), as well as indirbase and outdirbase. indirbase is the absolute path name to where your spill files are, up to and including "Spill". The
script will automatically append the spill number and ".root" to the file name when it processes.
Outdirbase is the output directory for your processed trigger files. The file names in that directory will be TriggerTreeSpill<NUMBER>.root. Again the scripts will append with <NUMBER> for you.

Once executed, MakeManyScripts.py will make a bash script that you can execute to process over all the necessary spill files.

Usage: python MakeManyScript.py, then bash Run.sh.



Note: This stage can take a few hours to process linearly. However, you can parallelize by executing the python script multiple times, each with a different range of (minspill,maxspill), and a different name for Run.sh (Run1.sh, Run2.sh, etc). Then each
Run.sh can be run separately on multiple nodes to parallelize and save time.

lariat_spill.C:

This script does a lot.
It loads in each of the 6 trees, and applies some time shifts.
First, a random time is drawn between 0 and 4.2 s, one for each mother id pion number, 350k RNG offsets in total. The accelerator timing structure informs this, such that the allowed
offsets match the timing structure from the accelerator.
Across all detectors, all daughters of that mother pion are offset by that time. This is the equivalent of starting the mother pion at t=tRNG.
Next, a cable delay is added, based on detector. USTOF is offset the most (18ns, by default), with smaller offsets for later detectors. This matches in data, where we added cable delays
such that the trigger card will see the hits across the detectors at the same time. If not used, USTOF hits would be ~20-40ns earlier than DSTOF hits, making matching difficult.

With these offsets done, the script loops over the combination of hits in the DSTOF, and for each opens a 100ns gate (detSigWidth[5]). Then for each hit in WC4, a 100ns gate is opened
(detSigWidth[4]). A coincidence is found if these gates overlap by 10ns (coincWidth). Once this match is made, hits in WC3 are checked in the same fashion with the DSTOF-WC4 matched pair.
This continues until a match is found in all 6 detectors. Similar code exists to create the WC2-missed and WC3-missed triggers.

It should be noted that we restrict possibly matched particle hits to ignore photons and neutrinos, as they are uncharged and leave no signal in the detectors. This greatly speeds up
processing time, as most of the particles found in the detector trees are photons.

If two or more triggers are too close in time, the readout windows will overlap. Particles in that overlapping time window will be stored in multiple triggers. In data, we protected
against this by implementing a dead time equal to the readout window. Similarly, this is done in sim, with the readout window =1024 samples at 1.117ns per sample. If two triggers are within a
readout window, the second trigger is ignored.

For each remaining trigger, information on all charged particles found within the corresponding readout window are saved in a tree for WC/TOF Reco.







