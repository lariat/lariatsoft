This directory contains all the FCL files used in the low-energy electron calorimetry analysis
using Michel electrons. Brief descriptions shown below. The "#" shown in file names signifies which
run period the script is for, either 1 or 2b in this case.
- Will Foreman (wforeman @ uchicago.edu)

///////////////////////////////
Data Scripts

The general flow is:
[digit-level data] --> Michel_Reco#.fcl --> Michel_Ana#.fcl

Michel_Reco#.fcl    : Reconstruction stages, including wire-by-wire calibrations, gaushit, 
                      a hit number filter (to exclude events with absurdly huge numbers of hits,
                      which are unlikely to be Michel electron events and which significantly slow
                      down later clustering stages), clustering, tracking, and calorimetry.

Michel_Ana#.fcl     : Runs the main analysis module, MichelAna_module. ALL settings related to this
                      module are stored as separate run-period-specific parameter sets in the file
                      michelana.fcl. The output is a ROOT file that contains the Michel analysis
                      TTree which is used in later plot-making macros stored in directory
                      LArIATAnaModule/MichelAnaMacros.

///////////////////////////////
MC Scripts

The general flow is:
[gen data] --> Michel_GenMC#.fcl --> Michel_RecoMC#.fcl --> Michel_AnaMC#.fcl

or if you want to run everything in a single go:
[gen data] --> MichelGenRecoAna#.fcl

Michel_GenMC#.fcl   : Generates through-going muons simulated with a cosmic-muon-like angular 
                      distribution. It is run in series with a filter module, MichelMCFilter_module,
                      which filters out only muons that stop and decay to electrons with a decay 
                      time of 300ns < T < 7300us to mimic the hardware trigger on double PMT flashes,
                      which uses a seris of gates to impose a delayed onset of the coincidence
                      window (see detector paper).

                      Generator-level MC muons (simb::MCTruth) have already been generated with the
                      appropriate angular distribution, and with an energy range of p=50-300 MeV/c to
                      maximize the probability of stopping (so NOT a true cosmic flux spectrum, 
                      which is OK for these purposes, as higher-energy muons pass right through the 
                      TPC without stopping.) These files are stored in the /pnfs/persistent area and
                      a file list must be fed as input to largeant:
                      /pnfs/lariat/persistent/users/wforeman/mcgen/cosmicmu/files.list

                      In case you're curious, the scripts and instructions for creating these 
                      generator-level files is found in
                      /pnfs/lariat/persistent/users/wforeman/mcgen/prodlist_scripts

Michel_RecoMC#.fcl  : Same reco procedure for data (sourced from Michel_Reco.fcl) except uses the
                      lariat_simulation_services, no wire-by-wire corrections, slightly different
                      hit-finding settings due to smaller hit widths in sim, and with MC-specific
                      calorimetry constants.

Michel_AnaMC#.fcl   : Same analysis module as for data (sourced from Michel_Ana.fcl) except uses
                      MC-specific settings. Use the "_nosmearing" version if you want to turn off
                      the artificial PMT smearing which was tuned to match data using separate
                      analysis macros stored in LArIATAnaModule/MichelAnaMacros.

Michel_GenRecoAna#.fcl: Does it all in one swoop!

michelana.fcl   :   Contains all the 


///////////////////////////////
Calibration Scripts

The general flow is:
[digit-level data] --> Michel_Ana#_calib.fcl (or opdetser_#.fcl)

Michel_Ana#_calib.fcl : This runs the main Michel analysis module, MichelAna_module, but turns
                        OFF features that would require reconstructed quantities, and filters
                        to look only at cosmic data. Averaged PMT waveforms are saved in the 
                        TTree.

opdetser_#.fcl        : Runs the photoelectron calibration to make SPE histograms.

