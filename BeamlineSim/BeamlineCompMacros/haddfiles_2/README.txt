Author: Greg Pulliam gkpullia@syr.edu

At this stage, you should have the raw G4BL output files from your jobs that have processed on the grid, stored (hopefully) on the pnfs scratch space.

However, each job output only processes a fraction of an actual spill, so you need to combine these files, usually in multiples of 10, to create a full spill for analysis.

This script loops over the directory containing the G4BL files, and merges them in groups of 10, and stores them for the next stage.
There are plenty of comments in the code itself, but in short:

Using a base file name, find all files with that base name, appended by the subspill number.
If a subspill file cannot be found, the spill associated to this subspill will be declared bad and won't be created. This ensures you only have full spills, and not partial
spills.

You may also declare subspill files bad as well, to have the spill file skipped. This is useful if the subspill file was created, but is broken in some way, usually indicated by
having a lower memory footprint than good subpsill files.

These lists of bad spills get printed for you.

For spills where all the subspills files exist and haven't been vetoed by you, the spill is created by hadding the appropriate subspill files.
If you already have spill files stored in the output directory, this script will not re-hadd them. This saves time by not hadding files that already exist.
This script will return the hadded files in the same directory as the script. Therefore, either cp Mergeroot.py to some place capable of storing ~100GB of simulation or edit the
script to send output there.
