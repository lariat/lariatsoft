Author: Greg Pulliam. gkpullia@syr.edu.

At this stage, G4BL events should have a reconstructed WCTrack and TOF.
The goal is to take these WCTracks and TOFs, apply WCQuality filters, then make composition plots by WCTrack condition: picky track or high yield

There are scripts for each of the four settings being used for pi analyses: Neg/Pos 60/100A. Each had a corresponding .C and .h file.

Some edits are required by user to make this work, primarily in the names of files the script will read in:

The .h file needs to be edited (~line 114) to point to the file where your reconstructed events are. 
The .C file needs to be edited (~line 146-154) to point to data files to compare to your simulation. However you've organized your data files, you need to edit these lines to match. 

Data information needed: WCPicky Tracks with WC-Filter. WCHY Tracks with WC-Filter. 3pt Tracks with WC-Filter. All TOF, no filtering. All mass, no filtering.

There are comments in the code for what each section does, but the short version is:

1. Loop over the Sim events. Ignore events that either dont pass the WCQuality cuts, isn't straight enough in the WCTrack straightness requirement, or doesn't have a TOF reco'd.
WCQuality/WCTrack cut values can be set at the top of the code. Set to values from +100A analysis (as of May 2019)

2. Smear momentum and TOF by some percentage defined at the top of the code. Set to 0 as default. You can smear the TOF by different amounts depending on TOF. A "Low" and "High" range for TOF smearing are allowed for now.

3. Calculate the mass from the momentum and TOF.

4. Fill momentum, TOF, mass plots for each particle species. In sim, an event is associated to a given particle species based on which particle made the hit in WC4 for the track.

5. With those in hand, make a bunch of plots comparing simulation to data. 


One plot that gets created is the fractional composition of pi/mu/e, as a function of momentum. This is useful for stage 6. For each of the three species, apply a second order polynomial fit over the momentum range, and note/save the fit parameters. They
can be used to create a more realistic DDMC. This plot uses the HY tracks, so if you want the picky track version, you need to change the script to use the picky track individual particle plots instead of the HY. 
