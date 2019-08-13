Code Author: Elena Gramellini (elenag@fnal.gov)

README Author and code contributor: Greg Pulliam (gkpullia@syr.edu)

NOTE: Most of this work was done by Elena, and the scripts in this folder are from her github (https://github.com/ElenaGramellini/LArIATSimpleGen), with a couple upgrades from me.


This stage uses distributions from data, coupled with the beamline composition analyze stage to generate MC events to simulate into the TPC. 
Instructions on Elena's GitHub explains the basis of this code, and should be the first reference. I'll review usage syntax for each stage and note changes I've made.

Step 1: python plotDataDistributions.py <filename> <Treename>
Makes histograms from anatree. Just a plotting script

Step 2: python XYMomentumGenTTree.py <filename> (optional)<mommin> (optional)<mommax>
Randomly draws XY,P,theta, phi, then checks whether that combination matches data, and creates a fills a tree when that occurs.

I've made an upgrade here. The initial momentum drawn is between 0-2000 MeV. We don't observe all those momenta in data, so this can be constrained by deciding a mommin and mommax to bracket the actual range of momentum in
the data sample being used as the probability distribution. This increases the probability a correct momentum is drawn, and increases the yield of simulated events for the next stage, without affecting the goal of the script.
Also, I've found occurances where a momentum in the tail of the data distribution is drawn, and the script gets stuck in a near infinite loop later on in the script when trying to pick a valid XY to go with the momentum. If
you notice that the script hangs in the while loop at line 142, it might be because of this. Therefore, you can change mommin and mommax to exclude the last few bins on the momentum spectra where there are only a couple fills
in data. You only need to exclude bins with fewer than 5-10 entries and the script will be quick to find a proper combination and avoid getting stuck. The more restrictive you are, the more momemtum bins you aren't properly
sampling, so only use this trick to remove the extreme, super rare bins in the tail.

Step 3: python generateHEPEvt.py [-h] [--NoProb] [--Prob] pdg [fname] [treeName]

Using the .root file created in the previous step, creates a large textfile to simulate particles of given pdg, in the HEPEvt format. 
Here is another upgrade and where the beamsim composition comes in. There is a flag you can set, --Prob or --NoProb, which will use the beamline composition function as a probability weighting. If you use --NoProb, you will
generate text files for as many entries as were in the tree created in the previous step, and the momentum spectrum of the generated particles will match as well. This is the "old" way of doing the simulation, and will
generate the same number of pi/mu/e/p

However, if you've run the beamline composition macro in stage 5, and have the composition fit functions, you can generate the events to match the momentum spectrum from the beamline composition study as well as represent the
true global composition. At the top of the script are 3 arrays to store the fit parameters from the composition fits. For example, the default parameters from when I performed the beamline compsosition in Aug 2019 had the
following for +60AHY:

piparams=[-2.06E-6,   3.07E-3,   -0.312]
muparams=[-1.01E-7,   1.55E-4,    0.037]
eparams= [ 2.16E-6,  -3.22E-3,     1.27]

This means that for pions, the composition function was -2.06E-6 p^2 + 3.07E-3 p -0.312, and similar functions for muons and electrons. 

If you generate other particles, don't use this, and instead keep --NoProb. The script will return an error if you aren't generating pi/mu/e and use --Prob.

If you reperform the composition you can change these numbers. Then by using the --Prob flag, it will use these parameters and weigh the particles properly. I've done the +60/+100A HY versions and have the default parameters. It
is left to the user to obtain the +60/+100 picky track parameters and -60/-100 picky/HY parameters. When complete, you will have a .root file with the histograms of all the generated particles, as well a text file of all the
particles generated.

Step 4: python SplitEvents.py <filename>
This is an add-on script. When you run the last step, you're going to have a text file with hundreds of thousands, if not millions of particles. However, as you will submit these to the grid to process, and the grid requires an
individual text file per job, the file from the previous step needs to be split into a bunch of individual text files, one for each grid job. 

As a default, we make each job on the grid generate 1k events. This script will take the input text file, and create a text file for each batch of 1k particles in the file. If there are N particles in the input text file, there will be
ceiling(N/1000) files created. Example, if you have 1001 particles in the input file, you will have 2 text files output, one with particles 0-999, and a second file with only 1 particle. 

NOTE: 
Delete that last text file if it doesn't have exactly 1k particles, which is likely (999/1000 likely). When you send these to the grid, you will specify in your xml the number of particles per textfile, which must be held
constant. The job will expect to see EXACTLY that many particles in the file. Therefore, if you submit the last text file to the grid, which will not have 1k particles, the job will fail. 


And that's it. You now have a bunch of text files, each representing 1k particles to generate. Stash them all in an empty /pnfs/scratch/, create a global list of all the textfiles using ls on that directory, and use that as a
file list to tell the grid what to process. That's all LArSoft stuff, which has tutorials and wikis for how to move forward, so I end the README here.

