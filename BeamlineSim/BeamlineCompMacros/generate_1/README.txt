Author: Greg Pulliam (gkpullia@syr.edu)
NOTE: READ COMPLETELY BEFORE EXECUTING.

This stage generates the raw Geant4BeamLine (g4bl) files that will used for the beamline compositon studies. 
While some machinery may be in place by other collaborators to automate more of the process, this file explains the process as it stands Aug 12th, 2019. 

In this folder are 4 .in files. Each .in file is the configuration file for each of the four current settings used for analyses: pos/neg 60/100 Amps.
For whichever current setting you are running, the process is similar.

As we generate a lot of simulation (350M pions on target per current setting), we parallelize the process on the grid.
For a given spill (4.2s of beam per supercycle), we will simulate 350k pions on target. However, that would take too long on a single grid node. Instead, each grid node processes 1/10 of
a spill (35k pions on target). Therefore, to create 350M pions on target at 35k pions on target per job, 10k jobs are needed.

The .in file is written to perform simulation of 35k pions on target. Therefore, this .in will be copied to each grid node.
However, as this isn't a larbatch submission (like with project.py), you need to supply a script to the grid node so that it sets up the necessary environment, executes the code and
returns the output to you. This is the purpose of Script.sh.

Script.sh sets up the environment to run, copies over the .in file, then uses the process number (0-9999) to figure out where to start the simulation in two ways.
1: G4BL uses a first and last to iterate how many pions to process. As a precaution, we require each pion to have a different ID number (0-350M). Using the process number, we uniquely
set first and last to make sense. For example: process 0 creates the first 35k pions (first=0, last=34999), which is the first tenth of the first spill. Process 1 (first=35k, last=74999) is the second tenth of
the first spill and so on to process 9999, which is the last tenth of spill 1000 (first=350M-35k to 350M-1).

2: using the process number, we set the RNG seed for the simulation. This is the purpose of "command1". It ensures the simulation across jobs will not be identical, which was a problem before, editing the .in
script to use the process number as the seed.

Finally, g4bl is exectuted, the 35k pions are simulated and returned to the user.

You (hopefully) notice that Script.sh needs some edits to work. This is where JobSubmit.sh is used.
JobSubmit.sh performs a couple edits to Script.sh (discussed below), and submits the jobs to the grid to be processed. It's a short script, and if you know sed, you can translate.
If not, each command does an edit to Script.sh

1. Replace the second argument with the first (explained below)
2. Replace "path" with the third argument
3-5. Fix the file suffix issue that is created by the first two commands.
6: Actually redunant. But I'll leave it so this readme makes sense.
7: Submit to the grid, using lariat as the group, 500MB of memory, a job-lifetime of 18h, 10k jobs, SL5/6 as the OS, using certain resources, and uses Script.sh as the script to process.
8-12: Undo all those previous seds.

Usage:
./JobSubmit.sh input <whatever>.in /some/directory/in/pnfs/scratch/space

Therefore, in Script.sh, "input" will be replaced with the input file, "path" will be replaced with the output directory you specify. You must create the output directory, using mkdir, before you submit.

Note: line 25 of script.sh copies the .in from scratch. Therefore, you must copy the .in to your output directory. As the script copies over the .in as part of its execution, the .in needs to be where the grid can see
it, and your output directory is as good a place as any, so you know what the input was when you get all the output.


Other things to be cognizant of, and why more improvements to automate all this is needed:

1. If you want to change the numbers for processing (particles/job, jobs/spill, total jobs), you need to know where all those numbers come into play with these scripts to change them, which probably will change file
names, leading to bullet 2.
2. When G4BL runs, the .in specifies the name of the file to create (line 14 of the .in) The script expects the following convention: outputfile =sim_lariat_<inputfile.in>.root, removing the .in. Be careful when playing
with file names, else the seds (sed 3-5) will fail.

As a suggestion, if you want to play with file names, stash the template script so you dont lose it if you mess up, comment out the seds after the job submission, run JobSubmit.sh and see if the Script.sh makes sense.
