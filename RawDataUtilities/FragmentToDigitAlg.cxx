////////////////////////////////////////////////////////////////
//                                                            //
// This is a class implementation for fragment-to-digit       //
// conversion.                                                //
//                                                            //
// Author: Brian Rebel, brebel@fnal.gov                       //
//                                                            //
//                                                            //
////////////////////////////////////////////////////////////////

// C++ includes
#include <bitset>

// LArIAT
#include "LArIATFragments/LariatFragment.h"
#include "LArIATFragments/WUTFragment.h"
#include "LArIATFragments/CAENFragment.h"
#include "LArIATFragments/V1495Fragment.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "Utilities/DatabaseUtilityT1034.h"

// LArSoft
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/AuxDetDigit.h"
#include "lardataobj/RawData/OpDetPulse.h"
#include "lardataobj/RawData/TriggerData.h"

#include "RawDataUtilities/FragmentToDigitAlg.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

enum {
  V1740_N_CHANNELS = 64,
  V1751_N_CHANNELS = 8,
  WUT_N_TDC_CHANNELS = 16,
  WUT_MAX_HITS = 128,
};


typedef std::map< int, std::map< int, std::vector< std::map< unsigned int, std::vector<unsigned int> > > > > match_maps;
typedef std::map< int, std::map< int, std::vector< std::pair<double, double> > > > fit_params_maps;
typedef std::vector<TDCFragment::TdcEventData> TDCDataBlock; 



//=====================================================================
//Constructor
FragmentToDigitAlg::FragmentToDigitAlg( fhicl::ParameterSet const& pset )
{
  this->reconfigure(pset);

}

//=====================================================================
FragmentToDigitAlg::~FragmentToDigitAlg()
{
}

//=====================================================================
void FragmentToDigitAlg::reconfigure( fhicl::ParameterSet const& pset )
{
  fRawFragmentLabel       = pset.get< std::string >("RawFragmentLabel",       "daq"  );
  fRawFragmentInstance    = pset.get< std::string >("RawFragmentInstance",    "SPILL");
  fTriggerDecisionTick    = pset.get< unsigned int>("TriggerDecisionTick",    135    ); 
  fTrigger1740Pedestal    = pset.get< float       >("Trigger1740Pedestal",    2000.  );
  fTrigger1740Threshold   = pset.get< float       >("Trigger1740Threshold",   0.     );
}

//=====================================================================
void FragmentToDigitAlg::makeTheDigits(std::vector<CAENFragment> caenFrags,
                                       std::vector<TDCDataBlock> tdcDataBlocks,
                                       std::vector<raw::AuxDetDigit> & auxDigits,
                                       std::vector<raw::RawDigit> & rawDigits,
                                       std::vector<raw::OpDetPulse> & opPulses )

{
  //Resetting the digits
  auxDigits.clear();
  rawDigits.clear();
  opPulses .clear();
  
  this->makeTPCRawDigits(caenFrags, rawDigits);
  this->makeOpDetPulses(caenFrags, opPulses);
  this->makeMuonRangeDigits(caenFrags, auxDigits);
  this->makeTOFDigits      (caenFrags, auxDigits);
  this->makeAeroGelDigits  (caenFrags, auxDigits);
  this->makeHaloDigits     (caenFrags, auxDigits);
  this->makeTriggerDigits  (caenFrags, auxDigits);
  for( size_t iDB = 0; iDB < tdcDataBlocks.size() ; ++iDB )
    this->makeMWPCDigits( tdcDataBlocks.at(iDB), auxDigits );

  return;  
}

//===============================================================----------
uint32_t FragmentToDigitAlg::triggerBits(std::vector<CAENFragment> const& caenFrags)
{

  // the trigger bits are piped into the V1740 board in slot 7, inputs 48 to 63
  // after run 6154 the bits were piped into a V1740 in slot 24, inputs 48 to 63
  // these are example connections as of May 08, 2015
  // 0   WC1      | OR of 2 X view TDCs ANDed with OR of 2 Y
  // 1   WC2      | "                                      " 
  // 2   WC3      | "                                      " 
  // 3   WC4      | "                                      " 
  // 4   BEAMON   | Spill gate : STARTs on $21, STOPs on $36 (cable says $26 but Bill says $36)
  // 5   USTOF    | OR of 4 PMTs
  // 6   DSTOF    | OR of 2 PMTs
  // 7   PUNCH    | OR of 2 X view paddles ANDed with OR of 2 Y
  // 8   HALO     | OR of 2 PMTs
  // 9   PULSER   |
  // 10  COSMICON | Cosmic gate : STARTs on $36, STOPs on $00 (not optimal, would like to stop before $00)
  // 11  COSMIC   | the trigger signal from the cosmic rack
  // 12  PILEUP   | Coincidence of any later LARSCINT with a delayed gate initiated by itself. Higher discrimination thresh. 
  // 13  MICHEL   | Coincidence of two light flashes in TPC (LARSCINT) occurring within a 5us time window
  // 14  LARSCINT | Coincidence of Hamamatsu and ETL PMTs (discriminated)
  // 15  MuRS     | Any coincidence of two planes.  Each plane is the OR of the discriminated pulses of 4 paddles. 

  // Each waveform corresponds to a single trigger channel.  If the (pedestal subtracted?) value of any ADC
  // in a waveform is less than 0, then the trigger for that channel fired

  // Need database eventually to set this correctly for different data-taking periods.
  // would set the fTriggerDecisionTick, fTrigger1740Pedestal, fTrigger1740Threshold values
  // in the beginRun method
  
  //std::bitset<16> triggerBits;
  std::bitset<32> triggerBits;

  size_t minChan  = 48;
  size_t maxChan  = 64;
  bool isBeamon = false;

  for(auto const& frag : caenFrags){

    // \todo Need to get this board ID information and run number from the database
    if     (frag.header.boardId != 7  && fRunNumber <  6155) continue;
    else if(frag.header.boardId != 24 && fRunNumber >= 6155) continue;

    for(size_t chan = minChan; chan < maxChan; ++chan){ 
      
      if(chan > frag.waveForms.size() )
      throw cet::exception("FragmentToDigitAlg") << "attempting to access channel "
						<< chan << " from 1740 fragment with only "
						<< frag.waveForms.size() << " channels";
       
      // require that waveform is large enough such that we don't encounter error
      // when trying to access time ticks around the fTriggerDecisionTick
      if(frag.waveForms[chan].data.size() > fTriggerDecisionTick+15 - 1) {
     
        // The trigger decision tick seems to jump around a bit over time, ranging between 125-145.  So,
        // let's check +/-15 samples around set fTriggerDecisionTick to be sure we don't miss it. (Considering 
        // the full drift is ~4000 ticks, any triggers occurring this close in time might as well be associated 
        // with the event anyway.)
        for(int i=-15; i<15; ++i){
          if( fTrigger1740Pedestal - frag.waveForms[chan].data[fTriggerDecisionTick+i] > fTrigger1740Threshold ){
            
            triggerBits.set(chan - minChan);
            LOG_VERBATIM("FragmentToDigitAlg")<<"*** FOUND TRIGGER BIT: "<<chan-minChan;
            
            // For Run 2 (Run > 8013), BEAMON was not being fed into any V1740 boards.  But it's input #16 on the 
            // V1495 so we'll just hack it by checking if any WC or TOF trigger bits were fired.
            if( (fRunNumber >= 8013) && ( (chan-minChan)==3 || (chan-minChan)==5 || (chan-minChan)==6 ) && !isBeamon) {
              triggerBits.set(16);
              LOG_VERBATIM("FragmentToDigitAlg")<<"*** FOUND TRIGGER BIT: 16";
              isBeamon = true;
            }
            
            break;   
          } // endIf over threshold
        } // end scan around neighborhood of fTriggerDecisionTick
      
      } // endif waveform is big enough
    
    } // end loop over channels on the board
    
  } // end loop over caen fragments

  return triggerBits.to_ulong();
}  

//===============================================================----------
void FragmentToDigitAlg::makeTPCRawDigits(std::vector<CAENFragment> const& caenFrags,
                                          std::vector<raw::RawDigit>     & tpcDigits)
{
  raw::ChannelID_t tpcChan = 0;
  size_t maxChan = 64;
  size_t boardId = 0;
  float  ped     = 0.;

  // make a list of the starting wire number for each board channel 0
  size_t startWireInd[8] = {239, 175, 111, 47,   0,   0,   0, 0 };
  size_t startWireCol[8] = {0,   0,   0,   239, 223, 159, 95, 31};
 
  for(auto const& frag : caenFrags){
    
    // the TPC mapping has the readout going to boards 0-7 of
    // the CAEN 1751, channels 0-63 of the boards 0-6, channels 0-31 of board 7
    // To make things hard, we decided to count the wires down instead of up
    // Board 0 channel 0  --> wire 239 of the induction plane
    // Board 3 channel 48 --> wire 0   of the induction plane
    // Board 3 channel 49 --> wire 239 of the collection plane
    // Board 7 channel 32 --> wire 0   of the collection plane
    boardId = frag.header.boardId;
    if(boardId > 7) continue;
    else{
      if(boardId < 7) maxChan = 64;
      else maxChan = 32;
      for(size_t chan = 0; chan < maxChan; ++chan){ 
        if(chan > frag.waveForms.size() )
        throw cet::exception("FragmentToDigitAlg") << "attempting to access channel "
						  << chan << " from 1740 fragment with only "
						  << frag.waveForms.size() << " channels";
        
        // get TPC channel for the induction plane
        if( boardId < 3 || (boardId == 3 && chan < 48) )
        tpcChan = startWireInd[boardId] - chan;
        // get TPC Channel for the collection plane
        else if( boardId > 3)
        tpcChan = 240 + startWireCol[boardId] - chan;
        else if(boardId == 3 && chan > 47)
        tpcChan = 240 + startWireCol[boardId] - chan + 48;
        
        // as of v04_13_00 of LArSoft, the event display no longer takes the
        // pedestal value from the RawDigit and uses an interface to a database instead
        // that doesn't really work for LArIAT, so pre-pedestal subtract the data
        // and keep the pedestal value for reference in the RawDigit
        std::vector<short> const padc(frag.waveForms[chan].data.begin(), frag.waveForms[chan].data.end());
        ped = this->findPedestal(padc);
        //REL	fRawDigitPedestals->Fill(ped);
        std::vector<short> adc(padc.size());
        for(size_t a = 0; a < adc.size(); ++a){
          adc[a] = padc[a] - (short)ped;
          //REL	  fRawDigitADC->Fill(adc[a]);
        }

        raw::RawDigit rd(tpcChan, adc.size(), adc);
        rd.SetPedestal(ped);
        tpcDigits.push_back(rd);
      } // end loop to fill channels from this board
    }// end if it is a TPC board      
  }// end loop over caen fragments

  return;
  
}

//===============================================================----------
float FragmentToDigitAlg::findPedestal(const std::vector<short> & adcVec)
{

  // do nothing if there are no values in the vector
  if(adcVec.size() < 1) return 0.;
  
  // first find overall mean and RMS
  float mean = 0.;
  float rms = 0.;
  float n = 1.*adcVec.size();
  for(auto const& adc : adcVec) mean += adc/n;
  for(auto const& adc : adcVec) rms += pow(1.*adc-mean,2);
  rms = sqrt( rms / n );
  
  // recalculate mean, but limiting to 
  // values within 1 sigma of overall mean
  float pedestal = 0.;
  int nn = 0;
  for(auto const& adc : adcVec){
    if( fabs(adc - mean) < rms ) {
      pedestal += 1.*adc;
      nn++;
    }
  }
  pedestal /= 1.*nn;
  
  return pedestal;

}

//===============================================================----------
void FragmentToDigitAlg::makeOpDetPulses(std::vector<CAENFragment>    const& caenFrags,
                                         std::vector<raw::OpDetPulse>      & opDetPulse)
{
  // loop over the caenFrags
  uint32_t boardId        = 0;
  uint32_t boardIdquery   = 0;
  uint32_t chanOff        = 0;

  // keywords for the light detection system (in retrospect we should have 
  // used some unique designator word like "lightsys_" for these channels...)
  const char *keywords[] = {"CRYO","SIPM","ARAPUCA"};
  std::vector< std::string > values(keywords,std::end(keywords));

  std::string board("board_");
  std::string channel("_channel_");
  size_t boardLoc;
  size_t channelLoc;
  size_t temp = -1; // initilize to max value
  uint32_t triggerTimeTag = 0;
  int      firstSample    = 0;

  for(auto const& caenFrag : caenFrags)
  {
    boardId        = caenFrag.header.boardId;
    triggerTimeTag = caenFrag.header.triggerTimeTag;

    for(auto hardwareIter : fHardwareConnections)
    {
      // search for one of the light-system keywords 
      // in this channel's name (from the hardware DB)
      for(size_t i=0; i<values.size(); i++){
        temp = hardwareIter.second.find(values[i]);
        if( temp != std::string::npos) break;
      }

      if(temp != std::string::npos)
      {
        boardLoc   = hardwareIter.first.find(board) + board.size();				//location in string where board ID is
        channelLoc = hardwareIter.first.find(channel) + channel.size();		//location in string where channel ID is
    
        if(boardLoc !=std::string::npos)
        {
          std::string strboardId (hardwareIter.first, boardLoc, channelLoc - channel.size());
          std::string strchannelId (hardwareIter.first, channelLoc, hardwareIter.first.size());
          
          boardIdquery = std::stoi(strboardId); 					//convert string boardID to uint32
          chanOff      = std::stoi(strchannelId); 			  //convert string channel to uint32

        }//end setting board/channel ID


        if(boardId == boardIdquery)
        {
          // loop over the channels on this board connected to optical detectors
          // check that the current channel, ch, is a valid one for grabbing a waveform
          if(chanOff > caenFrag.waveForms.size() )
          throw cet::exception("FragmentToDigitAlg")
            << "requested channel, " << chanOff
            << " from board "        << boardId
            << " is beyond the scope of the waveform vector";
        
          std::vector<short> waveForm(caenFrag.waveForms[chanOff].data.begin(), caenFrag.waveForms[chanOff].data.end());
        
          // calculate first sample
          firstSample = (int)((100.-fV1751PostPercent) * 0.01 * waveForm.size());
        
          LOG_VERBATIM("FragmentToDigitAlg")
          << " Found optical detector "
          << " boardID: " << boardId 
          << " channel: " << chanOff 
          << " wvform size: " << waveForm.size()
          << " Column: "<< hardwareIter.first<<" Value: "<<hardwareIter.second;
          //<< " BoardId: "<< boardId <<" ChannelId: "<< chanOff
          //<< " firstSample: "<< firstSample << " triggerTimeTag: " << triggerTimeTag;
        
          opDetPulse.push_back(raw::OpDetPulse(static_cast <unsigned short> (chanOff),
                                               waveForm,
                                               static_cast <unsigned int> (triggerTimeTag),
                                               static_cast <unsigned int> (firstSample)
                                               )
                               );
        
        } // end if this board has optical channels on it
      } // end if boardLoc found
    } // end loop over HardwareConnections
  }// end loop over fragments
  return;
}

//===============================================================----------
// boardId is the ID of the board we want to grab the digits from
// boardChans holds the channels on that board that we care about for this
// set of digits we want to make
// chanOffset is the value we subtract from the board channel so that the 
// digits have the right channel range for the desired auxiliary detector (ie 
// channel 0 of the muon range stack is not necessarily on channel 0 of the 
// caen board)
// detName is the name of the detector
void FragmentToDigitAlg::caenFragmentToAuxDetDigits(std::vector<CAENFragment>     const& caenFrags,
                                                    std::vector<raw::AuxDetDigit>      & auxDetDigits,
                                                    uint32_t                      const& boardId,
                                                    std::set<uint32_t>            const& boardChans,
                                                    uint32_t                      const& chanOffset,
                                                    std::string                   const& detName)
{
  // loop over the fragments and grab the one corresponding to this board ID
  for(auto const& frag : caenFrags){

    if(frag.header.boardId != boardId) continue;

    // loop over the channels in the set
    for( auto const& ch : boardChans){
      
      // check that ch is larger than chanOffset
      if(ch < chanOffset)
        throw cet::exception("FragmentToDigitAlg")
        << "requested channel, " << ch
        << " is smaller than the requested offest "
        << chanOffset;

      // check that there is a waveform for the chosen channel
      if(ch > frag.waveForms.size() )
        throw cet::exception("FragmentToDigitAlg") << "requested channel, " << ch
        << " from board "        << boardId
        << " is beyond the scope of the waveform vector";
      
      std::vector<short> waveForm(frag.waveForms[ch].data.begin(), frag.waveForms[ch].data.end());
	
      // place the AuxDetDigit in the vector
      auxDetDigits.push_back(raw::AuxDetDigit(static_cast<unsigned short> (ch - chanOffset),
                                              waveForm,
                                              detName,
                                              static_cast<unsigned long long>(frag.header.triggerTimeTag))
			     );

    } // end loop over channels on the board
  } // end loop over fragments
  
  return;
}

//===============================================================----------
void FragmentToDigitAlg::makeMuonRangeDigits(std::vector<CAENFragment>     const& caenFrags,
                                             std::vector<raw::AuxDetDigit>      & mrAuxDigits)
{
  uint32_t boardId = 0;
  uint32_t chanOff = 0;
  std::set<uint32_t> boardChans;
  std::set<std::string> MURSNames;
  MURSNames.insert("MURS1");
  MURSNames.insert("MURS2");
  MURSNames.insert("MURS3");
  MURSNames.insert("MURS4");
  MURSNames.insert("MURS5");
  MURSNames.insert("MURS6");
  MURSNames.insert("MURS7");
  MURSNames.insert("MURS8");
  MURSNames.insert("MURS9");
  MURSNames.insert("MURS10");
  MURSNames.insert("MURS11");
  MURSNames.insert("MURS12");
  MURSNames.insert("MURS13");
  MURSNames.insert("MURS14");
  MURSNames.insert("MURS15");
  MURSNames.insert("MURS16");
  MURSNames.insert("MURS17");
  MURSNames.insert("MURS18");
  MURSNames.insert("MURS19");
  MURSNames.insert("MURS20");
  MURSNames.insert("MURS21");
  MURSNames.insert("MURS22");
  MURSNames.insert("MURS23");
  MURSNames.insert("MURS24");
  MURSNames.insert("MURS25");
  MURSNames.insert("MURS26");
  MURSNames.insert("MURS27");
  MURSNames.insert("MURS28");
  MURSNames.insert("MURS29");
  MURSNames.insert("MURS30");
  MURSNames.insert("MURS31");
  MURSNames.insert("MURS32");
  std::string board("board_");									//looking for board ID number in string
  std::string channel("_channel_");						//looking for channel ID number in string
  size_t boardLoc;
  size_t channelLoc;
  for( auto hardwareIter : fHardwareConnections)
  {
      // there is at most one entry in a set with a given value, so
      // if we find any entries with the name given by the value of
      // hardwareItr, we are looking at the device of interest
    if( MURSNames.count(hardwareIter.second) > 0)
    {
      
      boardLoc   = hardwareIter.first.find(board) + board.size();				//location in string where board ID is
      channelLoc = hardwareIter.first.find(channel) + channel.size();		//location in string where channel ID is
      
      if(boardLoc != std::string::npos)
      {
        std::string strboardId  (hardwareIter.first, boardLoc,   channelLoc - channel.size());
        std::string strchannelId(hardwareIter.first, channelLoc, hardwareIter.first.size());
        
        boardId = std::stoi(strboardId); 								//convert string boardID to uint32
        chanOff = std::stoi(strchannelId); 							//convert string channel to uint32
        LOG_VERBATIM("FragmentToDigitAlg")
        <<" Found MURS"
        <<" Column: "<< hardwareIter.first
        <<" Value: "<<hardwareIter.second
        <<" BoardId: "<< boardId
        <<" ChannelId: "<< chanOff;
      }//end setting board/channel ID
      
      // \todo This value for chanOff may not be correct.  The chanOff value is intended to allow
      // the caenFragmentToAuxDetDigit know how many to subtract from the board channel such that the
      // channels from the auxdet all satisfy the range of 0-N.
      boardChans.insert(chanOff);
    }//end if MURS
  }//end loop over HardwareConnections
  LOG_VERBATIM("FragmentToDigitAlg")
  << "ChannelOffset: " << *(boardChans.begin());
  
  this->caenFragmentToAuxDetDigits(caenFrags, mrAuxDigits, boardId, boardChans,*(boardChans.begin()) , "MuonRangeStack");
  boardChans.clear();

  return;
}

//===============================================================----------
void FragmentToDigitAlg::makeTOFDigits(std::vector<CAENFragment>     const& caenFrags,
                                       std::vector<raw::AuxDetDigit>      & tofAuxDigits)
{
  uint32_t boardId = 0;
  uint32_t chanOff = 0;
  std::set<uint32_t> boardChansUS;
  std::set<uint32_t> boardChansDS;
  std::set<std::string> TOFNames;
  TOFNames.insert("USTOF1");
  TOFNames.insert("USTOF2");
  TOFNames.insert("USTOF3");
  TOFNames.insert("USTOF4");
  TOFNames.insert("DSTOF1");
  TOFNames.insert("DSTOF2");
  TOFNames.insert("DSTOF3");
  TOFNames.insert("DSTOF4");
  std::string board("board_");
  std::string channel("_channel_");
  size_t boardLoc;
  size_t channelLoc;
  for( auto hardwareIter : fHardwareConnections) 
  {
    if( TOFNames.count(hardwareIter.second) > 0)
    {
      boardLoc   = hardwareIter.first.find(board) + board.size();				//location in string where board ID is
      channelLoc = hardwareIter.first.find(channel) + channel.size();				//location in string where channel ID is
      
      if(boardLoc !=std::string::npos)
      {
        std::string strboardId (hardwareIter.first, boardLoc, channelLoc - channel.size());
        std::string strchannelId (hardwareIter.first, channelLoc, hardwareIter.first.size());
        
        boardId = std::stoi(strboardId);								//convert string boardID to uint32
        chanOff = std::stoi(strchannelId); 							//convert string channel to uint32
        LOG_VERBATIM("FragmentToDigitAlg")
        <<" Found TOF"
        <<" Column: "   << hardwareIter.first
        <<" Value: "    << hardwareIter.second
        <<" BoardId: "  << boardId
        <<" ChannelId: "<< chanOff;
      }//end setting board/channel ID

      // \todo These values for chanOff may not be correct.  The chanOff value is intended to allow
      // the caenFragmentToAuxDetDigit know how many to subtract from the board channel such that the
      // channels from the auxdet all satisfy the range of 0-N.
      if     ( hardwareIter.second == "USTOF1" || hardwareIter.second == "USTOF2" || hardwareIter.second == "USTOF3" || hardwareIter.second == "USTOF4" ) boardChansUS.insert(chanOff);
      else if( hardwareIter.second == "DSTOF1" || hardwareIter.second == "DSTOF2" || hardwareIter.second == "DSTOF3" || hardwareIter.second == "DSTOF4" ) boardChansDS.insert(chanOff);
    }//end find TOFNames
  }//end loop over hardwareDatabase

  LOG_VERBATIM("FragmentToDigitAlg")
  << "ChannelOffset: " << *(boardChansUS.begin())
  << " ChannelOffset: " << *(boardChansDS.begin());
  
  this->caenFragmentToAuxDetDigits(caenFrags, tofAuxDigits, boardId, boardChansUS, *(boardChansUS.begin()), "TOFUS");
  this->caenFragmentToAuxDetDigits(caenFrags, tofAuxDigits, boardId, boardChansDS, *(boardChansDS.begin()), "TOFDS");
  
  boardChansUS.clear();
  boardChansDS.clear();

          
  return;
}

//===============================================================----------
void FragmentToDigitAlg::makeAeroGelDigits(std::vector<CAENFragment>     const& caenFrags,
                                           std::vector<raw::AuxDetDigit>      & agAuxDigits)
{

  uint32_t boardId1p06_1 = 0;
  uint32_t boardId1p06_2 = 0;
  uint32_t boardId1p10_1 = 0;
  uint32_t boardId1p10_2 = 0;

  uint32_t boardId = 0;
  uint32_t chanOff = 0;

  std::set<uint32_t> chanAG1p10_1;
  std::set<uint32_t> chanAG1p10_2;
  std::set<uint32_t> chanAG1p06_1;
  std::set<uint32_t> chanAG1p06_2;

  std::set<std::string> AGNames;
  AGNames.insert("AG1p06_1"); //AGUSE
  AGNames.insert("AG1p06_2"); //AGUSW
  AGNames.insert("AG1p10_1"); //AGDSE
  AGNames.insert("AG1p10_2"); //AGDSW
  std::string board("board_");
  std::string channel("_channel_");
  size_t boardLoc;
  size_t channelLoc;
   for( auto hardwareIter : fHardwareConnections) 
   {
     if( AGNames.count(hardwareIter.second) > 0)
     {
       
       boardLoc   = hardwareIter.first.find(board) + board.size();				//location in string where board ID is
       channelLoc = hardwareIter.first.find(channel) + channel.size();				//location in string where channel ID is
       
       if(boardLoc !=std::string::npos)
       {
         std::string strboardId (hardwareIter.first, boardLoc, channelLoc - channel.size());
         std::string strchannelId (hardwareIter.first, channelLoc, hardwareIter.first.size());
	 boardId = std::stoi(strboardId);

	 // Have to create a unique variable for each board
	 // because it doesn't stay consistent over time
	 // and the two Aerogel counters aren't always on the same boards
	 if     (hardwareIter.second == "AG1p10_1") boardId1p10_1 = boardId;
	 else if(hardwareIter.second == "AG1p10_2") boardId1p10_2 = boardId;
	 else if(hardwareIter.second == "AG1p06_1") boardId1p06_1 = boardId;
	 else if(hardwareIter.second == "AG1p06_2") boardId1p06_2 = boardId;

         chanOff = std::stoi(strchannelId);

         LOG_VERBATIM("FragmentToDigitAlg")<< " Found AeroGel"
         << " Column: "  << hardwareIter.first
         << " Value: "   << hardwareIter.second
         << " BoardId: " << boardId
         <<" ChannelId: "<< chanOff;
       }//end setting board/channel ID
       
        // \todo These values for chanOff may not be correct.  The chanOff value is intended to allow
        // the caenFragmentToAuxDetDigit know how many to subtract from the board channel such that the
        // channels from the auxdet all satisfy the range of 0-N.
       if     (hardwareIter.second == "AG1p10_1") chanAG1p10_1.insert(chanOff);
       else if(hardwareIter.second == "AG1p10_2") chanAG1p10_2.insert(chanOff);
       else if(hardwareIter.second == "AG1p06_1") chanAG1p06_1.insert(chanOff);
       else if(hardwareIter.second == "AG1p06_2") chanAG1p06_2.insert(chanOff);

     }//end find AGNAMES
  }//end loop over hardwareDatabase

  LOG_VERBATIM("FragmentToDigitAlg")
  << "ChannelOffset: " <<  *(chanAG1p10_1.begin())
  << " ChannelOffset: " <<  *(chanAG1p10_2.begin())
  << " ChannelOffset: " <<  *(chanAG1p06_1.begin())
  << " ChannelOffset: " <<  *(chanAG1p06_2.begin());


  this->caenFragmentToAuxDetDigits(caenFrags, agAuxDigits, boardId1p10_1, chanAG1p10_1, *(chanAG1p10_1.begin()), "AG1p10_1");
  this->caenFragmentToAuxDetDigits(caenFrags, agAuxDigits, boardId1p10_2, chanAG1p10_2, *(chanAG1p10_2.begin()), "AG1p10_2");

  this->caenFragmentToAuxDetDigits(caenFrags, agAuxDigits, boardId1p06_1, chanAG1p06_1, *(chanAG1p06_1.begin()), "AG1p06_1");
  this->caenFragmentToAuxDetDigits(caenFrags, agAuxDigits, boardId1p06_2, chanAG1p06_2, *(chanAG1p06_2.begin()), "AG1p06_2");
  
  chanAG1p10_1.clear();
  chanAG1p10_2.clear();
  chanAG1p06_1.clear();
  chanAG1p06_2.clear();

  return;
}

//===============================================================----------
// Halo paddles are currently (Jun 4, 2015) attached to board 9, channels 5 and 6
void FragmentToDigitAlg::makeHaloDigits(std::vector<CAENFragment>     const& caenFrags,
                                        std::vector<raw::AuxDetDigit>      & hAuxDigits)
{
   uint32_t boardId = 0;
   uint32_t chanOff = 0;
   std::set<uint32_t> boardChans;
   std::set<std::string> HALONames;
   HALONames.insert("HALO1");
   HALONames.insert("HALO2");
   std::string board("board_");									//looking for board ID number in string
   std::string channel("_channel_");								//looking for channel ID number in string
   size_t boardLoc;
   size_t channelLoc;
  for( auto hardwareIter : fHardwareConnections) 
  {
    if( HALONames.count(hardwareIter.second) > 0)
    {
      
      boardLoc   = hardwareIter.first.find(board) + board.size();				//location in string where board ID is
      channelLoc = hardwareIter.first.find(channel) + channel.size();				//location in string where channel ID is
      
      if(boardLoc !=std::string::npos)
      {
        std::string strboardId (hardwareIter.first, boardLoc, channelLoc - channel.size());
        std::string strchannelId (hardwareIter.first, channelLoc, hardwareIter.first.size());
        
        boardId = std::stoi(strboardId); 								//convert string boardID to uint32
        chanOff = std::stoi(strchannelId); 							//convert string channel to uint32
        
        LOG_VERBATIM("FragmentToDigitAlg")
        <<" Found HALO"
        <<" Column: "   << hardwareIter.first
        <<" Value: "    << hardwareIter.second
        <<" BoardId: "  << boardId
        <<" ChannelId: "<< chanOff;
        
      }//end setting board/channel ID
      
      // \todo This value for chanOff may not be correct.  The chanOff value is intended to allow
      // the caenFragmentToAuxDetDigit know how many to subtract from the board channel such that the
      // channels from the auxdet all satisfy the range of 0-N.
      boardChans.insert(chanOff);
    }//end find HALONames
  }//end loop over HardwareConnections

  LOG_VERBATIM("FragmentToDigitAlg") << "ChannelOffset: " << *(boardChans.begin());
  
  this->caenFragmentToAuxDetDigits(caenFrags, hAuxDigits, boardId, boardChans, *(boardChans.begin()), "Halo");
  boardChans.clear();
  return;
}

//===============================================================----------
void FragmentToDigitAlg::makeTriggerDigits(std::vector<CAENFragment>     const& caenFrags,
                                           std::vector<raw::AuxDetDigit>      & trAuxDigits)
{
  // For runs <  6155: trigger waveforms all come on board 7, channels 48-63.
  // For runs >= 6155: trigger waveforms come on board 24, channels 48-63 and ...(0-15?)
  uint32_t boardId=0;
  uint32_t chanOff=0;
  std::set<uint32_t> boardChans;
  std::string board("board_");
  std::string channel("_channel_");
  size_t boardLoc;
  size_t channelLoc;
  std::set<std::string> trigNames;
  trigNames.insert("WC1");
  trigNames.insert("WC2");    
  trigNames.insert("WC3");      
  trigNames.insert("WC4");    
  trigNames.insert("BEAMON"); 
  trigNames.insert("COSMICAGCOIN");
  trigNames.insert("USTOF");  
  trigNames.insert("DSTOF");    
  trigNames.insert("PUNCH");  
  trigNames.insert("HALO");   
  trigNames.insert("PULSER"); 
  trigNames.insert("COSMICON"); 
  trigNames.insert("COSMIC"); 
  trigNames.insert("LARRY");
  trigNames.insert("PILEUP"); 
  trigNames.insert("MICHEL"); 
  trigNames.insert("LARSCINT"); 
  trigNames.insert("MURS");
  for( auto hardwareIter : fHardwareConnections) 
  {
    if( trigNames.count(hardwareIter.second) > 0)
    {
      boardLoc   = hardwareIter.first.find(board) + board.size();
      channelLoc = hardwareIter.first.find(channel) + channel.size();
      
      if(boardLoc !=std::string::npos)
      {
        std::string strboardId (hardwareIter.first, boardLoc, channelLoc - channel.size());
        std::string strchannelId (hardwareIter.first, channelLoc, hardwareIter.first.size());
        
        boardId = std::stoi(strboardId);
        chanOff = std::stoi(strchannelId);
        LOG_VERBATIM("FragmentToDigitAlg") <<" Found Triggers"
        <<" Column: "    << hardwareIter.first
        <<" Value: "     << hardwareIter.second
        <<" BoardId: "   << boardId
        <<" ChannelId: " << chanOff
        <<" TrigNames: " << hardwareIter.second;
      }//end setting board/channel ID
      
      // \todo This value for chanOff may not be correct.  The chanOff value is intended to allow
      // the caenFragmentToAuxDetDigit know how many to subtract from the board channel such that the
      // channels from the auxdet all satisfy the range of 0-N.
      boardChans.insert(chanOff);
      
      LOG_VERBATIM("FragmentToDigitAlg") << "ChannelOffset: " << *(boardChans.begin());
      
      this->caenFragmentToAuxDetDigits(caenFrags, trAuxDigits, boardId, boardChans, *(boardChans.begin()), hardwareIter.second);
      boardChans.clear();
    }//end loop over trigNames
  }//end loop over HardwareConnections


  return;
}

//===============================================================----------
// The map below indicates how each TDC maps to each Wire Chamber
// channel            wires
//   0   |-----------| 1
// TDC 3 |           |
//   63  |           |
//       |           | Wire Chamber 1
//   0   |           |
// TDC 4 |           |
//   63  |-----------| 128
//       TDC1     TDC2
//       0  63   0  63 channel
// wires 1         128

// channel            wires
//   0   |-----------| 1
// TDC 7 |           |
//   63  |           |
//       |           | Wire Chamber 2
//   0   |           |
// TDC 8 |           |
//   63  |-----------| 128
//       TDC5     TDC6
//       0  63   0  63 channel
// wires 1         128

// channel            wires
//   0   |-----------| 1
// TDC 11|           |
//   63  |           |
//       |           | Wire Chamber 3
//   0   |           |
// TDC 12|           |
//   63  |-----------| 128
//       TDC9     TDC10
//       0  63   0  63 channel
// wires 1         128

// channel            wires
//   0   |-----------| 1
// TDC 15|           |
//   63  |           |
//       |           | Wire Chamber 4
//   0   |           |
// TDC 16|           |
//   63  |-----------| 128
//       TDC13    TDC14
//       0  63   0  63 channel
// wires 1         128
//
// take the convention that vertical wire numbers start at channel 128
void FragmentToDigitAlg::InitializeMWPCContainers()
{
  this->CleanUpMWPCContainers();

  // make the map of TDC number to detector name and tdc to starting channel
  for(size_t tdc = 1; tdc < 17; ++tdc){
    if(tdc < 5)       fTDCToChamber[tdc] = 0;
    else if(tdc < 9)  fTDCToChamber[tdc] = 1;
    else if(tdc < 13) fTDCToChamber[tdc] = 2;
    else              fTDCToChamber[tdc] = 3;

    if     (tdc == 1 || tdc == 5 || tdc == 9  || tdc == 13) fTDCToStartWire[tdc] = 0;
    else if(tdc == 2 || tdc == 6 || tdc == 10 || tdc == 14) fTDCToStartWire[tdc] = 64;
    else if(tdc == 3 || tdc == 7 || tdc == 11 || tdc == 15) fTDCToStartWire[tdc] = 128;
    else if(tdc == 4 || tdc == 8 || tdc == 12 || tdc == 16) fTDCToStartWire[tdc] = 192;
  }

  // had swapped cables in runs 5546 - 5598 for TDC 7 and 8.
  if(fRunNumber > 5545 && fRunNumber < 5599){
    fTDCToStartWire[7] = 192;
    fTDCToStartWire[8] = 128;
  }

  fMWPCNames.resize(4);
  fMWPCNames[0] = "MWPC1";
  fMWPCNames[1] = "MWPC2";
  fMWPCNames[2] = "MWPC3";
  fMWPCNames[3] = "MWPC4";

  return;
}

//===============================================================----------
void FragmentToDigitAlg::CleanUpMWPCContainers()
{
  fMWPCNames     .clear();
  fTDCToStartWire.clear();   
  fTDCToChamber  .clear();

  return;
}

//===============================================================----------
// set the name of the detector in the AuxDetDigit to be of the form
// MWPCXX where XX is the controller Number
void FragmentToDigitAlg::makeMWPCDigits(std::vector<TDCFragment::TdcEventData> const& tdcEventData,
                                        std::vector<raw::AuxDetDigit>               & mwpcAuxDigits)
{

  size_t channelsPerChamber = TDCFragment::N_CHANNELS * TDCFragment::TDCS_PER_CHAMBER;

  // vector to hold the channels for a single MWPC
  std::vector<std::vector<short> > chamberHits(TDCFragment::MAX_CHAMBERS * channelsPerChamber);

  // vector to hold the timeStamps for each channel in the MWPC
  std::vector<unsigned long long> chamberTimeStamps(TDCFragment::MAX_CHAMBERS * channelsPerChamber, 0);

  // LOG_VERBATIM("FragmentToDigitAlg") << "there are " << tdcEventData.size() << " tdcEventData objects in the vector";

  for(auto const& tdced : tdcEventData){ 

    unsigned int tdcNumber = static_cast <unsigned int> (tdced.tdcEventHeader.tdcNumber);

    mf::LogDebug("FragmentToDigitAlg") << "tdcEventHeader.tdcNumber: " << tdcNumber;

    // determine the chamber and start wire
    auto switr = fTDCToStartWire.find(tdced.tdcEventHeader.tdcNumber);
    auto chitr = fTDCToChamber.find(tdced.tdcEventHeader.tdcNumber);

    // return if TDC number is out of range; stops filling of MWPC digits
    if (tdcNumber < 1 or tdcNumber > 16) return;

    if( chitr == fTDCToChamber.end() || switr == fTDCToStartWire.end() )
      throw cet::exception("FragmentToDigitAlg")
      << "TDC number " << tdcNumber
      << " is not present in map to chamber number or start wire";

    // LOG_VERBATIM("FragmentToDigitAlg") << "there are " << tdced.tdcHits.size() << " tdc hit objects in the vector for chamber "
    // 				    << chitr->second << " on tdc " << chitr->first << " start wire " << switr->second;
    
    for(auto const& hit : tdced.tdcHits){
      if(chitr->second >= TDCFragment::MAX_CHAMBERS || 
         switr->second + (size_t)hit.channel >= channelsPerChamber
         )
        throw cet::exception("FragmentToDigitAlg")
        << "Chamber is " << chitr->second << "/" << TDCFragment::MAX_CHAMBERS
        << " hit channel is " << (size_t)hit.channel
        << " first wire in tdc is " << switr->second << "/"
        << channelsPerChamber;

      chamberHits      [chitr->second * channelsPerChamber + switr->second + size_t (hit.channel)].push_back(hit.timeBin);
      chamberTimeStamps[chitr->second * channelsPerChamber + switr->second + size_t (hit.channel)] = tdced.tdcEventHeader.tdcTimeStamp;

      // LOG_VERBATIM("FragmentToDigitAlg") << chamberHits[chitr->second * channelsPerChamber + switr->second + size_t (hit.channel)].size() << " " 
      // 				      << (size_t)hit.channel << " " << switr->second << " " << (size_t)hit.timeBin << "\t" 
      // 				      << chamberTimeStamps[chitr->second * channelsPerChamber + switr->second + size_t (hit.channel)] << " " 
      // 				      << tdced.tdcEventHeader.tdcTimeStamp;
	
    }
      
  } // end loop over tdcEventData

  // now make the AuxDetDigits for this fragment
  for(size_t cham = 0; cham < TDCFragment::MAX_CHAMBERS; ++cham){
    for(size_t chan = 0; chan < channelsPerChamber; ++chan){

      if(chamberHits[cham*channelsPerChamber + chan].size() < 1) continue;
      
      mwpcAuxDigits.push_back(raw::AuxDetDigit(static_cast <unsigned short> (chan),
                                               chamberHits[cham*channelsPerChamber + chan],
                                               fMWPCNames[cham],
                                               static_cast <unsigned long long> (chamberTimeStamps[cham*channelsPerChamber + chan]))
                              );

    }
    
  } // end loops to create AuxDetDigits

  return;
}

//===============================================================----------
std::vector<raw::Trigger> FragmentToDigitAlg::makeTheTriggers(art::EventNumber_t                                    const& EventNumber,
                                                              std::vector<CAENFragment>                             const& caenFrags,
                                                              std::vector< std::vector<TDCFragment::TdcEventData> > const& tdcDataBlocks)
{
  //Hardcoded for now until I can find out how to
  //extract this from the raw event root file
  float eventTime = 0;

  bool caenDataPresent = false;
  std::vector<raw::Trigger> triggerVector;

  //If there are caen fragments present, set the trigger info
  //based on the caen data block information
  if (caenFrags.size() > 0) {
    triggerVector.push_back(raw::Trigger(EventNumber, caenFrags.front().header.triggerTimeTag, eventTime, this->triggerBits(caenFrags)));
    caenDataPresent = true;
  }
  else
    LOG_WARNING("FragmentToDigit")
    << "There are no CAEN Fragments for event " << EventNumber
    << " that may be OK, so continue";

  //If there are no caen fragments present, then set the trigger info
  //based on the TDCDataBlock information
  if (tdcDataBlocks.size() > 0) {
    if (!caenDataPresent) {
      triggerVector.push_back(raw::Trigger(EventNumber, tdcDataBlocks.front().front().tdcEventHeader.tdcTimeStamp,
                              eventTime, this->triggerBits(caenFrags)));
    }
  }
  else
    LOG_WARNING("FragmentToDigit")
    << "There are no TDC Fragments for event " << EventNumber
    << " that may be OK, so continue";

  return triggerVector;

}

//===============================================================----------
//void FragmentToDigitAlg::InitializeRun(art::RunPrincipal* const& run, art::RunNumber_t runNumber, uint64_t timestamp)
void FragmentToDigitAlg::InitializeRun(art::RunNumber_t runNumber, uint64_t timestamp)
{
  fRunNumber=runNumber;
  fRunDateTime = this->TimestampToString(timestamp);	//jess lines
  this->InitializeMWPCContainers();
  LOG_VERBATIM("FragmentToDigitAlg") << "Initializing Run"
                                     << " RunNumber: "     << fRunNumber
                                     <<  "; Date/Time: "   << fRunDateTime
                                     << "; RunTimestamp: " << fRunTimestamp;

  // Set config parameters to get from the lariat_prd database
  fConfigParams.clear();
  fConfigParams.push_back("v1495_config_v1495_delay_ticks");
  fConfigParams.push_back("v1740_config_caen_postpercent");
  fConfigParams.push_back("v1740_config_caen_recordlength");
  fConfigParams.push_back("v1740b_config_caen_postpercent");
  fConfigParams.push_back("v1740b_config_caen_recordlength");
  fConfigParams.push_back("v1751_config_caen_postpercent");
  fConfigParams.push_back("v1751_config_caen_recordlength");
  fConfigParams.push_back("v1740_config_caen_v1740_samplereduction");
  fConfigParams.push_back("v1740b_config_caen_v1740_samplereduction");
  fConfigParams.push_back("tdc_config_tdc_pipelinedelay");
  fConfigParams.push_back("tdc_config_tdc_gatewidth");
  
  // Get V1751 PostPercent settings from database
  fConfigValues.clear();
  fHardwareConnections.clear();											//jess lines
  fHardwareConnections = fDatabaseUtility->GetHardwareConnections(fRunDateTime);			        //jess lines
  fConfigValues = fDatabaseUtility->GetConfigValues(fConfigParams, static_cast <int> (fRunNumber));		
  fV1751PostPercent = std::atof(fConfigValues["v1751_config_caen_postpercent"].c_str());

}
//-------------------jess lines
//=====================================================================
  std::string FragmentToDigitAlg::TimestampToString(std::time_t const& Timestamp) {
    struct tm * TimeInfo;
    char Buffer[30];
    TimeInfo = std::localtime(&Timestamp);
    std::strftime(Buffer, 30, "%Y-%m-%d %H:%M", TimeInfo);
    return std::string(Buffer);
  }

