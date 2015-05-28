////////////////////////////////////////////////////////////////////////
/// \file  ChannelMapLArIATAlg.cxx
/// \brief Interface to algorithm class for the standar, simplest detector channel mapping
///
/// \version $Id:  $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

// ART Includes
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h" 

// LArSoft Includes
#include "Geometry/GeometryCore.h"
#include "Geometry/AuxDetGeo.h"
#include "Geometry/AuxDetSensitiveGeo.h"
#include "Geometry/CryostatGeo.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "Geometry/WireGeo.h"

// LArIATSoft
#include "Geo/ChannelMapLArIATAlg.h"

// C++ standard library
#include <ostream> // std::endl


namespace geo{

  //----------------------------------------------------------------------------
  ChannelMapLArIATAlg::ChannelMapLArIATAlg(fhicl::ParameterSet const& p)
    : fSorter(geo::GeoObjectSorterLArIAT(p))
  {
  }

  //----------------------------------------------------------------------------
  void ChannelMapLArIATAlg::Initialize( GeometryData_t& geodata )
  {
    // start over:
    Uninitialize();
    
    std::vector<geo::CryostatGeo*>& cgeo = geodata.cryostats;
    std::vector<geo::AuxDetGeo*>  & adgeo = geodata.auxDets;
    
    fNcryostat = cgeo.size();
    
    mf::LogInfo("ChannelMapLArIATAlg") << "Initializing LArIAT ChannelMap...";

    // First sort the LArTPC related geometry objects and get their channel mapping
    fSorter.SortCryostats(cgeo);
    for(size_t c = 0; c < cgeo.size(); ++c) 
      cgeo[c]->SortSubVolumes(fSorter);
    
    fNTPC.resize(fNcryostat);
    fWireCounts.resize(fNcryostat);
    fNPlanes.resize(fNcryostat);
    fFirstWireProj.resize(fNcryostat);
    fOrthVectorsY.resize(fNcryostat);
    fOrthVectorsZ.resize(fNcryostat);
    fPlaneBaselines.resize(fNcryostat);
    fWiresPerPlane.resize(fNcryostat);
    fFirstChannelInNextPlane.resize(fNcryostat);
    fFirstChannelInThisPlane.resize(fNcryostat);
    fViews.clear();
    fPlaneIDs.clear();
    fTopChannel = 0;

    int RunningTotal = 0;

    for(unsigned int cs = 0; cs != fNcryostat; ++cs){
      
      fNTPC[cs] = cgeo[cs]->NTPC();
      
      // Size up all the vectors 
      fWireCounts[cs]             .resize(fNTPC[cs]);
      fFirstWireProj[cs] 	  .resize(fNTPC[cs]);
      fOrthVectorsY[cs]  	  .resize(fNTPC[cs]);
      fOrthVectorsZ[cs]  	  .resize(fNTPC[cs]);
      fPlaneBaselines[cs]	  .resize(fNTPC[cs]);
      fWiresPerPlane[cs] 	  .resize(fNTPC[cs]);
      fNPlanes[cs] 	  	  .resize(fNTPC[cs]);
      fFirstChannelInThisPlane[cs].resize(fNTPC[cs]);
      fFirstChannelInNextPlane[cs].resize(fNTPC[cs]);

      for(unsigned int TPCCount = 0; TPCCount != fNTPC[cs]; ++TPCCount){
	unsigned int PlanesThisTPC = cgeo[cs]->TPC(TPCCount).Nplanes();
	fWireCounts[cs][TPCCount]   .resize(PlanesThisTPC);
	fFirstWireProj[cs][TPCCount].resize(PlanesThisTPC);
	fOrthVectorsY[cs][TPCCount] .resize(PlanesThisTPC);
	fOrthVectorsZ[cs][TPCCount] .resize(PlanesThisTPC);
	fNPlanes[cs][TPCCount]=PlanesThisTPC;
	for(unsigned int PlaneCount = 0; PlaneCount != PlanesThisTPC; ++PlaneCount){

	  fViews.emplace(cgeo[cs]->TPC(TPCCount).Plane(PlaneCount).View());
	  fPlaneIDs.emplace(PlaneID(cs, TPCCount, PlaneCount));
	  double ThisWirePitch = cgeo[cs]->TPC(TPCCount).WirePitch(0, 1, PlaneCount);
	  fWireCounts[cs][TPCCount][PlaneCount] = cgeo[cs]->TPC(TPCCount).Plane(PlaneCount).Nwires();
	  
	  double  WireCentre1[3] = {0.,0.,0.};
	  double  WireCentre2[3] = {0.,0.,0.};
	  
	  const geo::WireGeo& firstWire = cgeo[cs]->TPC(TPCCount).Plane(PlaneCount).Wire(0);
	  const double sth = firstWire.SinThetaZ(), cth = firstWire.CosThetaZ();
	  
	  firstWire.GetCenter(WireCentre1,0);
	  cgeo[cs]->TPC(TPCCount).Plane(PlaneCount).Wire(1).GetCenter(WireCentre2,0);
	  
	  // figure out if we need to flip the orthogonal vector 
	  // (should point from wire n -> n+1)
	  double OrthY = cth, OrthZ = -sth;
	  if(((WireCentre2[1] - WireCentre1[1])*OrthY 
	      + (WireCentre2[2] - WireCentre1[2])*OrthZ) < 0){
	    OrthZ *= -1;
	    OrthY *= -1;
	  }
	  
	  // Overall we are trying to build an expression that looks like
	  //  int NearestWireNumber = round((worldPos.OrthVector - FirstWire.OrthVector)/WirePitch);     
	  // That runs as fast as humanly possible.
	  //
	  // Casting to an int is much faster than c rounding commands like floor().  We have to add 0.5
	  // to account for rounding up as well as down.  floor(A) ~ (int)(A+0.5).  We account for the
	  // 0.5 in the first wire constant to avoid adding it every time.
	  //
	  // We can also predivide everything by the wire pitch so we don't do this in the loop
	  //
	  // Putting this together into the useful constants we will use later per plane and tpc:
	  fOrthVectorsY[cs][TPCCount][PlaneCount] = OrthY / ThisWirePitch;
	  fOrthVectorsZ[cs][TPCCount][PlaneCount] = OrthZ / ThisWirePitch;
	  
	  fFirstWireProj[cs][TPCCount][PlaneCount]  = WireCentre1[1]*OrthY + WireCentre1[2]*OrthZ;
	  fFirstWireProj[cs][TPCCount][PlaneCount] /= ThisWirePitch;
	  
	  // now to count up wires in each plane and get first channel in each plane
	  int WiresThisPlane = cgeo[cs]->TPC(TPCCount).Plane(PlaneCount).Nwires();	
	  fWiresPerPlane[cs] .at(TPCCount).push_back(WiresThisPlane);
	  fPlaneBaselines[cs].at(TPCCount).push_back(RunningTotal);
	  
	  RunningTotal += WiresThisPlane;

	  fFirstChannelInThisPlane[cs].at(TPCCount).push_back(fTopChannel);
	  fTopChannel += WiresThisPlane;
	  fFirstChannelInNextPlane[cs].at(TPCCount).push_back(fTopChannel);

	}// end loop over planes
      }// end loop over TPCs
    }// end loop over cryostats

    // calculate the total number of channels in the detector
    fNchannels = fTopChannel;

    LOG_DEBUG("ChannelMapLArIAT") << "# of channels is " << fNchannels;

    // now sort the AuxDetGeo objects and map them to names of the detectors
    // Each raw::AuxDetDigit knows the name of the detector it came from and its
    // channel - sometimes channel maps to sensitive volume and some times it doesn't
    fSorter.SortAuxDets(adgeo);
    for(auto a : adgeo) a->SortSubVolumes(fSorter);

    // map the AuxDetGeo names to their position in the sorted vector
    // Each TOF  detector has 2   channels, and 1  sensitive volume
    // Each AG   detector has 2   channels, and 1  sensitive volume
    // The  MuRS detector has 16  channels, and 16 sensitive volumes
    // Each MWPC detector has 128 channels, and 1  sensitive volume
    fADNameToGeo.clear();
    fADChannelToSensitiveGeo.clear();

    for(size_t a = 0; a < adgeo.size(); ++a){
      std::string volName(adgeo[a]->TotalVolume()->GetName());

      if(volName.find("TOFUS") != std::string::npos){
	fADNameToGeo["TOFUS"] = a;
	for(size_t c = 0; c < 2; ++c) fADChannelToSensitiveGeo[a].push_back(0);
      }
      else if(volName.find("TOFDS") != std::string::npos){
	fADNameToGeo["TOFDS"] = a;
	for(size_t c = 0; c < 2; ++c) fADChannelToSensitiveGeo[a].push_back(0);
      }
      else if(volName.find("MWPC1") != std::string::npos){
	fADNameToGeo["MWPC1"]= a;
	for(size_t c = 0; c < 128; ++c) fADChannelToSensitiveGeo[a].push_back(0);
      }
      else if(volName.find("MWPC2") != std::string::npos){
	fADNameToGeo["MWPC2"]= a;
	for(size_t c = 0; c < 128; ++c) fADChannelToSensitiveGeo[a].push_back(0);
      }
      else if(volName.find("MWPC3") != std::string::npos){
	fADNameToGeo["MWPC3"]= a;
	for(size_t c = 0; c < 128; ++c) fADChannelToSensitiveGeo[a].push_back(0);
      }
      else if(volName.find("MWPC4") != std::string::npos){
	fADNameToGeo["MWPC4"]= a;
	for(size_t c = 0; c < 128; ++c) fADChannelToSensitiveGeo[a].push_back(0);
      }
      else if(volName.find("AeroGelUS") != std::string::npos){
	fADNameToGeo["AeroGelUS"] = a;
	for(size_t c = 0; c < 2; ++c) fADChannelToSensitiveGeo[a].push_back(0);
      }	
      else if(volName.find("AeroGelDS") != std::string::npos){
	fADNameToGeo["AeroGelDS"] = a;
	for(size_t c = 0; c < 2; ++c) fADChannelToSensitiveGeo[a].push_back(0);
      }	
      else if(volName.find("MuonRangeStack") != std::string::npos){
	fADNameToGeo["MuonRangeStack"] = a;

	// check that there are the expected 16 sensitive volumes
	if(adgeo[a]->NSensitiveVolume() != 16)
	  throw cet::exception("ChannelMapLArIATAlg") << "Expected 16 sensitive volumes for MuRS AuxDet, "
						      << volName << " instead have only " 
						      << adgeo[a]->NSensitiveVolume();

	// the sorting should have sorted them by plane going from bottom to top in a plane
	
	
      }

    }

    for(auto nitr : fADNameToGeo){

      

    return;
  }
   
  //----------------------------------------------------------------------------
  void ChannelMapLArIATAlg::Uninitialize()
  {
  }

  //----------------------------------------------------------------------------
  std::vector<geo::WireID> ChannelMapLArIATAlg::ChannelToWire(raw::ChannelID_t channel)  const
  {
    std::vector< geo::WireID > AllSegments; 
    unsigned int wire  = 0;
    
    // first check if this channel ID is legal
    if(channel > fTopChannel)
      throw cet::exception("Geometry") << "ILLEGAL CHANNEL ID for channel " << channel << "\n";
    
    // then go find which plane, tpc and cryostat it is in from the information we stored earlier
    for(auto const& id : fPlaneIDs ){
      
      if(channel < fFirstChannelInNextPlane[id.Cryostat][id.TPC][id.Plane]){
    	wire  = channel - fFirstChannelInThisPlane[id.Cryostat][id.TPC][id.Plane];
	AllSegments.push_back(geo::WireID(id.Cryostat, id.TPC, id.Plane, wire));
    	break;
      }	    

    }// end loop over PlaneIDs


    return AllSegments;
  }

  //----------------------------------------------------------------------------
  unsigned int ChannelMapLArIATAlg::Nchannels() const
  {
    return fNchannels;
  }

  //----------------------------------------------------------------------------
  double ChannelMapLArIATAlg::WireCoordinate
    (double YPos, double ZPos, geo::PlaneID const& planeID) const
  {
    // Returns the wire number corresponding to a (Y,Z) position in PlaneNo 
    // with float precision.
    // B. Baller August 2014
    return YPos*AccessElement(fOrthVectorsY, planeID) 
         + ZPos*AccessElement(fOrthVectorsZ, planeID)
         - AccessElement(fFirstWireProj, planeID);
  }

  
  //----------------------------------------------------------------------------
  WireID ChannelMapLArIATAlg::NearestWireID
    (const TVector3& worldPos, geo::PlaneID const& planeID) const
  {

    // This part is the actual calculation of the nearest wire number, where we assume
    //  uniform wire pitch and angle within a wireplane
    
    // add 0.5 to have the correct rounding
    int NearestWireNumber = int
      (0.5 + WireCoordinate(worldPos.Y(), worldPos.Z(), planeID));
    
    // If we are outside of the wireplane range, throw an exception
    // (this response maintains consistency with the previous
    // implementation based on geometry lookup)
    if(NearestWireNumber < 0 ||
       NearestWireNumber >= AccessElement(fWireCounts, planeID))
    {
      int wireNumber = NearestWireNumber; // save for the output
      
      if(NearestWireNumber < 0 ) NearestWireNumber = 0;
      else                       NearestWireNumber = AccessElement(fWireCounts, planeID) - 1;
      
      throw InvalidWireIDError("Geometry", wireNumber, NearestWireNumber)
        << "Can't Find Nearest Wire for position (" 
        << worldPos[0] << "," << worldPos[1] << "," << worldPos[2] << ")"
        << " approx wire number # " << wireNumber
        << " (capped from " << NearestWireNumber << ")\n";
    }

    return geo::WireID(planeID, (geo::WireID::WireID_t) NearestWireNumber);

  }
  
  //----------------------------------------------------------------------------
  // This method returns the channel number, assuming the numbering scheme
  // is heirachical - that is, channel numbers run in order, for example:
  //                                             (Ben J Oct 2011)                   
  //                    Wire1     | 0
  //           Plane1 { Wire2     | 1
  //    TPC1 {          Wire3     | 2
  //           Plane2 { Wire1     | 3   increasing channel number
  //                    Wire2     | 4     (with no gaps)
  //    TPC2 { Plane1 { Wire1     | 5
  //           Plane2 { Wire1     | 6
  //                    Wire2     v 7
  //
  raw::ChannelID_t ChannelMapLArIATAlg::PlaneWireToChannel
    (geo::WireID const& wireID) const
  {
    // This is the actual lookup part - first make sure coordinates are legal
    if(wireID.TPC   < fNTPC[wireID.Cryostat] &&
       wireID.Plane < fWiresPerPlane[wireID.Cryostat][wireID.TPC].size() &&
       wireID.Wire  < AccessElement(fWiresPerPlane, wireID)){
      // if the channel has legal coordinates, its ID is given by the wire
      // number above the number of wires in lower planes, tpcs and cryostats;
      // wireID is used here as PlaneID
      return AccessElement(fPlaneBaselines, wireID) + wireID.Wire;
    }
    else{  
      // if the coordinates were bad, throw an exception
      throw cet::exception("ChannelMapLArIATAlg") << "NO CHANNEL FOUND for " << std::string(wireID);
    }
    
    // made it here, that shouldn't happen, return raw::InvalidChannelID
    mf::LogWarning("ChannelMapLArIATAlg") << "should not be at the point in the function, returning "
					    << "invalid channel";
    return raw::InvalidChannelID;

  }


  //----------------------------------------------------------------------------
  SigType_t ChannelMapLArIATAlg::SignalType(raw::ChannelID_t const channel) const
  {

    // still assume one cryostat for now -- faster
    unsigned int nChanPerTPC = fNchannels/fNTPC[0];
    // casting wil trunc towards 0 -- faster than floor
    unsigned int tpc = channel / nChanPerTPC;  
    //need number of planes to know Collection 
    unsigned int PlanesThisTPC = fNPlanes[0][tpc];
    
      
    
    SigType_t sigt = geo::kMysteryType;
    if(      (channel >= fFirstChannelInThisPlane[0][tpc][0]) &&
             (channel <  fFirstChannelInNextPlane[0][tpc][PlanesThisTPC-2])    ){ sigt = geo::kInduction; }
    else if( (channel >= fFirstChannelInThisPlane[0][tpc][PlanesThisTPC-1]) &&
             (channel <  fFirstChannelInNextPlane[0][tpc][PlanesThisTPC-1])    ){ sigt = geo::kCollection; }
    else
      mf::LogWarning("BadChannelSignalType") << "Channel " << channel
					     << " not given signal type." << std::endl;
            

    return sigt;
  }


  //----------------------------------------------------------------------------
  View_t ChannelMapLArIATAlg::View(raw::ChannelID_t const channel) const
  {

    // still assume one cryostat for now -- faster
    unsigned int nChanPerTPC = fNchannels/fNTPC[0];
    // casting wil trunc towards 0 -- faster than floor
    unsigned int tpc = channel / nChanPerTPC; 


    View_t view = geo::kUnknown; 

    if(      (channel >= fFirstChannelInThisPlane[0][tpc][0]) &&
             (channel <  fFirstChannelInNextPlane[0][tpc][0])    ){ view = geo::kU; }
    else if( (channel >= fFirstChannelInThisPlane[0][tpc][1]) &&
             (channel <  fFirstChannelInNextPlane[0][tpc][1])    ){ view = geo::kV; }
    else if( (channel >= fFirstChannelInThisPlane[0][tpc][2]) &&
             (channel <  fFirstChannelInNextPlane[0][tpc][2])    ){ view = geo::kZ; }
    else
      mf::LogWarning("BadChannelSignalType") << "Channel " << channel
					     << " not given view type.";

    return view;
  }  

  //----------------------------------------------------------------------------
  std::set<View_t> const& ChannelMapLArIATAlg::Views() const
  {
    return fViews;
  }

  //----------------------------------------------------------------------------
  std::set<PlaneID> const& ChannelMapLArIATAlg::PlaneIDs() const
  {
    return fPlaneIDs;
  }

} // namespace
