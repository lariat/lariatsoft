Author: Greg Pulliam (gkpullia@syr.edu)

This README gives an overview of the beamline simulation process. In each directory, there is another README that gives more details about that particular stage.

The goal of these stages is generate G4 Beamline similation (G4BL), merge files together, process for triggers, reconstruct a WCTrack and TOF for each trigger, then analyze the WCTracks to create a beamline composition at WC4.

Stage 1: Generation

This stage creates the raw G4BL files, and stores them in /pnfs/scratch/. If you use the default parameters, this will create 10k root TTree files.

Stage 2: Hadding

This stage takes those 10k root files, and combines them to make a "spill", which is the equivalent to 4.2s of beam per supercycle, the equivalent of a data "subrun". Using default parameters, this will merge in groups of 10, creating
1000 spill files, which can be stored in /lariat/data temporarily, as the memory footprintis large.

Stage 3: Trigger Finding

This stage loops over the spill, performs some timing offsets to get beam timing correct, then sorts through the particles in the spill to find coincidences across detectors, simulating the trigger process in data. TTree files are stored in
/lariat/data, one file per spill, and the file size is small enough to leave there long term.

Stage 4: Track Building

For all the triggers in a spill, reconstruct the WCTrack and TOF using the sim information, storing in a TTree, one file per spill. File size is small enough that they can be stored in /lariat/data long term.

Stage 5: Analyze

With all the tracks reconstructed, hadd them together, and run the appropriate script to return the beamline composition for that sample.
