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
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"

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
  void ChannelMapLArIATAlg::Initialize( GeometryData_t const& geodata )
  {
    // start over:
    Uninitialize();
    
    std::vector<geo::CryostatGeo*> const& cgeo = geodata.cryostats;
    
    fNcryostat = cgeo.size();
    
    LOG_DEBUG("ChannelMapLArIATAlg")
    << "Initializing LArIAT ChannelMap..."
    << "\t There are "
    << fNcryostat
    << " cryostats";

    fNTPC.resize(fNcryostat);
    fWireCounts             .resize(fNcryostat);
    fNPlanes                .resize(fNcryostat);
    fFirstWireProj          .resize(fNcryostat);
    fOrthVectorsY           .resize(fNcryostat);
    fOrthVectorsZ           .resize(fNcryostat);
    fPlaneBaselines         .resize(fNcryostat);
    fWiresPerPlane          .resize(fNcryostat);
    fFirstChannelInNextPlane.resize(fNcryostat);
    fFirstChannelInThisPlane.resize(fNcryostat);
    fViews.clear();
    fPlaneIDs.clear();
    fTopChannel = 0;

    int RunningTotal = 0;

    for(unsigned int cs = 0; cs != fNcryostat; ++cs){
      
      fNTPC[cs] = cgeo[cs]->NTPC();
      
      LOG_DEBUG("ChannelMapLArIATAlg")
      << "\t Cryostat "
      << cs
      << " has "
      << fNTPC[cs]
      << " TPCs";

      // Size up all the vectors 
      fWireCounts[cs]             .resize(fNTPC[cs]);
      fFirstWireProj[cs] 	        .resize(fNTPC[cs]);
      fOrthVectorsY[cs]  	        .resize(fNTPC[cs]);
      fOrthVectorsZ[cs]  	        .resize(fNTPC[cs]);
      fPlaneBaselines[cs]	        .resize(fNTPC[cs]);
      fWiresPerPlane[cs] 	        .resize(fNTPC[cs]);
      fNPlanes[cs] 	  	          .resize(fNTPC[cs]);
      fFirstChannelInThisPlane[cs].resize(fNTPC[cs]);
      fFirstChannelInNextPlane[cs].resize(fNTPC[cs]);

      for(unsigned int TPCCount = 0; TPCCount != fNTPC[cs]; ++TPCCount){
        
        unsigned int PlanesThisTPC = cgeo[cs]->TPC(TPCCount).Nplanes();
        
        LOG_DEBUG("ChannelMapLArIATAlg")
        << "\t TPC "
        << TPCCount
        << " has "
        << PlanesThisTPC
        << " planes";
        
        fWireCounts[cs][TPCCount]   .resize(PlanesThisTPC);
        fFirstWireProj[cs][TPCCount].resize(PlanesThisTPC);
        fOrthVectorsY[cs][TPCCount] .resize(PlanesThisTPC);
        fOrthVectorsZ[cs][TPCCount] .resize(PlanesThisTPC);
        fNPlanes[cs][TPCCount]=PlanesThisTPC;
        for(unsigned int PlaneCount = 0; PlaneCount != PlanesThisTPC; ++PlaneCount){

          fViews.emplace(cgeo[cs]->TPC(TPCCount).Plane(PlaneCount).View());
          fPlaneIDs.emplace(PlaneID(cs, TPCCount, PlaneCount));

          // the following commented out code is for debugging purposes
//          size_t nwires = cgeo[cs]->TPC(TPCCount).Plane(PlaneCount).Nwires();
//          for(size_t w = 0; w < nwires; ++w){
//            auto const& wire1 = cgeo[cs]->TPC(TPCCount).Plane(PlaneCount).Wire(w);
//            auto const& wire2 = (w == nwires - 1) ? cgeo[cs]->TPC(TPCCount).Plane(PlaneCount).Wire(w - 1) : cgeo[cs]->TPC(TPCCount).Plane(PlaneCount).Wire(w + 1);
//            
//            double pitch = geo::WireGeo::WirePitch(wire1, wire2);
//            double xyz[3] = {0.};
//            wire1.GetCenter(xyz);
//            
//            if(std::abs(pitch - 0.5) > 1.e-3)
//              LOG_VERBATIM("ChannelMapLArIATAlg")
//              << "\t\t wire: "
//              << w
//              << " pitch: "
//              << pitch
//              << " center ("
//              << xyz[0]
//              << ", "
//              << xyz[1]
//              << ", "
//              << xyz[2]
//              << ")";
//          }
          
          double ThisWirePitch = cgeo[cs]->TPC(TPCCount).WirePitch(0, 1, PlaneCount);
          
          fWireCounts[cs][TPCCount][PlaneCount] = cgeo[cs]->TPC(TPCCount).Plane(PlaneCount).Nwires();

          LOG_DEBUG("ChannelMapLArIATAlg")
          << "\t Plane "
          << PlaneCount
          << " has "
          << fWireCounts[cs][TPCCount][PlaneCount]
          << " wires";
          
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
  unsigned int ChannelMapLArIATAlg::Nchannels
    (readout::ROPID const& ropid) const
  {
    if (!HasROP(ropid)) return 0;
    // The number of channels matches the number of wires. Life is easy.
    return WireCount(FirstWirePlaneInROP(ropid));
  } // ChannelMapLArIATAlg::Nchannels(ROPID)
  
  
  //----------------------------------------------------------------------------
  double ChannelMapLArIATAlg::WireCoordinate(double              YPos, 
                                             double              ZPos,
                                             geo::PlaneID const& planeID) const
  {
    // Returns the wire number corresponding to a (Y,Z) position in PlaneNo 
    // with float precision.
    // B. Baller August 2014
    return YPos*AccessElement(fOrthVectorsY, planeID) 
         + ZPos*AccessElement(fOrthVectorsZ, planeID)
         - AccessElement(fFirstWireProj, planeID);
  }

  
  //----------------------------------------------------------------------------
  WireID ChannelMapLArIATAlg::NearestWireID(TVector3     const& worldPos, 
                                            geo::PlaneID const& planeID) const
  {

    // This part is the actual calculation of the nearest wire number, where we assume
    //  uniform wire pitch and angle within a wireplane
    
    // add 0.5 to have the correct rounding
    int NearestWireNumber = int(0.5 + WireCoordinate(worldPos.Y(), worldPos.Z(), planeID));
    
    // If we are outside of the wireplane range, throw an exception
    // (this response maintains consistency with the previous
    // implementation based on geometry lookup)
    if(NearestWireNumber < 0 || (unsigned int)NearestWireNumber >= this->WireCount(planeID) )
    {
      int wireNumber = NearestWireNumber; // save for the output
      
      if(NearestWireNumber < 0 ) NearestWireNumber = 0;
      else                       NearestWireNumber = this->WireCount(planeID) - 1;
      
      throw InvalidWireIDError("Geometry", wireNumber, NearestWireNumber)<< "Can't Find Nearest Wire for position (" 
									 << worldPos[0] << "," 
									 << worldPos[1] << "," 
									 << worldPos[2] << ")"
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
  raw::ChannelID_t ChannelMapLArIATAlg::PlaneWireToChannel(geo::WireID const& wireID) const
  {
    unsigned int const* pBaseLine = GetElementPtr(fPlaneBaselines, wireID);
    // This is the actual lookup part - first make sure coordinates are legal
    if (pBaseLine) {
      // if the channel has legal coordinates, its ID is given by the wire
      // number above the number of wires in lower planes, tpcs and cryostats
      return *pBaseLine + wireID.Wire;
    }
    else{  
      // if the coordinates were bad, throw an exception
      throw cet::exception("ChannelMapLArIATAlg") << "NO CHANNEL FOUND for " << std::string(wireID);
    }
    
    // made it here, that shouldn't happen, return raw::InvalidChannelID
    LOG_WARNING("ChannelMapLArIATAlg")
    << "should not be at the point in the function, returning "
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
      LOG_WARNING("BadChannelSignalType")
      << "Channel " << channel
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
      LOG_WARNING("BadChannelSignalType")
      << "Channel " << channel
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


  //----------------------------------------------------------------------------
  unsigned int ChannelMapLArIATAlg::NTPCsets
    (readout::CryostatID const& cryoid) const
  {
    // return the same number as the number of TPCs
    return (cryoid.isValid && cryoid.Cryostat < fNTPC.size())?
      fNTPC[cryoid.Cryostat]: 0;
  } // ChannelMapLArIATAlg::NTPCsets()
  
  
  //----------------------------------------------------------------------------
  unsigned int ChannelMapLArIATAlg::MaxTPCsets() const
  {
    return MaxTPCs();
  } // ChannelMapLArIATAlg::MaxTPCsets()
  
  
  //----------------------------------------------------------------------------
  bool ChannelMapLArIATAlg::HasTPCset(readout::TPCsetID const& tpcsetid) const
  {
    return tpcsetid.TPCset < NTPCsets(tpcsetid);
  } // ChannelMapLArIATAlg::HasTPCset()
  
  
  //----------------------------------------------------------------------------
  readout::TPCsetID ChannelMapLArIATAlg::TPCtoTPCset
    (geo::TPCID const& tpcid) const
  {
    return ConvertTPCtoTPCset(tpcid);
  } // ChannelMapLArIATAlg::TPCtoTPCset()
  
  
  //----------------------------------------------------------------------------
  std::vector<geo::TPCID> ChannelMapLArIATAlg::TPCsetToTPCs
    (readout::TPCsetID const& tpcsetid) const
  {
    std::vector<geo::TPCID> IDs;
    if (tpcsetid.isValid) IDs.emplace_back(ConvertTPCsetToTPC(tpcsetid));
    return IDs;
  } // ChannelMapLArIATAlg::TPCsetToTPCs()
  
  
  //----------------------------------------------------------------------------
  geo::TPCID ChannelMapLArIATAlg::FirstTPCinTPCset
    (readout::TPCsetID const& tpcsetid) const
  {
    return ConvertTPCsetToTPC(tpcsetid);
  } // ChannelMapLArIATAlg::FirstTPCinTPCset()
  
  
  //----------------------------------------------------------------------------
  unsigned int ChannelMapLArIATAlg::MaxTPCs() const
  {
    unsigned int max = 0;
    for (unsigned int nTPCs: fNTPC) if (nTPCs > max) max = nTPCs;
    return max;
  } // ChannelMapLArIATAlg::MaxTPCs()
  
  
  //----------------------------------------------------------------------------
  unsigned int ChannelMapLArIATAlg::NROPs
      (readout::TPCsetID const& tpcsetid) const
  {
    if (!HasTPCset(tpcsetid)) return 0;
    return AccessElement(fNPlanes, FirstTPCinTPCset(tpcsetid));
  } // ChannelMapLArIATAlg::NROPs()
  
  
  //----------------------------------------------------------------------------
  unsigned int ChannelMapLArIATAlg::MaxROPs() const {
    unsigned int max = 0;
    for (auto const& cryo_tpc: fNPlanes)
      for (unsigned int nPlanes: cryo_tpc)
        if (nPlanes > max) max = nPlanes;
    return max;
  } // ChannelMapLArIATAlg::MaxROPs()
  
  
  //----------------------------------------------------------------------------
  bool ChannelMapLArIATAlg::HasROP(readout::ROPID const& ropid) const {
    return ropid.ROP < NROPs(ropid);
  } // ChannelMapLArIATAlg::HasROP()
  
  
  //----------------------------------------------------------------------------
  readout::ROPID ChannelMapLArIATAlg::WirePlaneToROP
    (geo::PlaneID const& planeid) const
  {
    return ConvertWirePlaneToROP(planeid);
  } // ChannelMapLArIATAlg::WirePlaneToROP()
  
  
  //----------------------------------------------------------------------------
  std::vector<geo::PlaneID> ChannelMapLArIATAlg::ROPtoWirePlanes
    (readout::ROPID const& ropid) const
  {
    std::vector<geo::PlaneID> IDs;
    if (ropid.isValid) IDs.emplace_back(FirstWirePlaneInROP(ropid));
    return IDs;
  } // ChannelMapLArIATAlg::ROPtoWirePlanes()
  
  
  //----------------------------------------------------------------------------
  std::vector<geo::TPCID> ChannelMapLArIATAlg::ROPtoTPCs
    (readout::ROPID const& ropid) const
  {
    std::vector<geo::TPCID> IDs;
    // we take the TPC set of the ROP and convert it straight into a TPC ID
    if (ropid.isValid) IDs.emplace_back(ConvertTPCsetToTPC(ropid.asTPCsetID()));
    return IDs;
  } // ChannelMapLArIATAlg::ROPtoTPCs()
  
  
  //----------------------------------------------------------------------------
  readout::ROPID ChannelMapLArIATAlg::ChannelToROP
    (raw::ChannelID_t channel) const
  {
    if (!raw::isValidChannelID(channel)) return {}; // invalid ROP returned
    
    // which wires does the channel cover?
    std::vector<geo::WireID> wires = ChannelToWire(channel);
    
    // - none:
    if (wires.empty()) return {}; // default-constructed ID, invalid
    
    // - one: maps its plane ID into a ROP ID
    return WirePlaneToROP(wires[0]);
  } // ChannelMapLArIATAlg::ROPtoTPCs()
  
  
  //----------------------------------------------------------------------------
  raw::ChannelID_t ChannelMapLArIATAlg::FirstChannelInROP
    (readout::ROPID const& ropid) const
  {
    if (!ropid.isValid) return raw::InvalidChannelID;
    return (raw::ChannelID_t)
      AccessElement(fPlaneBaselines, ConvertROPtoWirePlane(ropid));
  } // ChannelMapLArIATAlg::FirstChannelInROP()
  
  
  //----------------------------------------------------------------------------
  geo::PlaneID ChannelMapLArIATAlg::FirstWirePlaneInROP
    (readout::ROPID const& ropid) const
  {
    return ConvertROPtoWirePlane(ropid);
  } // ChannelMapLArIATAlg::FirstWirePlaneInROP()
  
  
  //----------------------------------------------------------------------------
  readout::TPCsetID ChannelMapLArIATAlg::ConvertTPCtoTPCset
    (geo::TPCID const& tpcid)
  {
    if (!tpcid.isValid) return {}; // invalid ID, default-constructed
    return {
      (readout::CryostatID::CryostatID_t) tpcid.Cryostat,
      (readout::TPCsetID::TPCsetID_t) tpcid.TPC
      };
  } // ChannelMapLArIATAlg::ConvertTPCtoTPCset()
  
  
  //----------------------------------------------------------------------------
  geo::TPCID ChannelMapLArIATAlg::ConvertTPCsetToTPC
    (readout::TPCsetID const& tpcsetid)
  {
    if (!tpcsetid.isValid) return {};
    return {
      (geo::CryostatID::CryostatID_t) tpcsetid.Cryostat,
      (geo::TPCID::TPCID_t) tpcsetid.TPCset
      };
  } // ChannelMapLArIATAlg::ConvertTPCsetToTPC()
  
  
  //----------------------------------------------------------------------------
  readout::ROPID ChannelMapLArIATAlg::ConvertWirePlaneToROP
    (geo::PlaneID const& planeid)
  {
    if (!planeid.isValid) return {}; // invalid ID, default-constructed
    return {
      (readout::CryostatID::CryostatID_t) planeid.Cryostat,
      (readout::TPCsetID::TPCsetID_t) planeid.TPC,
      (readout::ROPID::ROPID_t) planeid.Plane
      };
    
  } // ChannelMapLArIATAlg::ConvertWirePlaneToROP()
  
  
  //----------------------------------------------------------------------------
  geo::PlaneID ChannelMapLArIATAlg::ConvertROPtoWirePlane
    (readout::ROPID const& ropid)
  {
    if (!ropid.isValid) return {};
    return {
      (geo::CryostatID::CryostatID_t) ropid.Cryostat,
      (geo::TPCID::TPCID_t) ropid.TPCset,
      (geo::PlaneID::PlaneID_t) ropid.ROP
      };
  } // ChannelMapLArIATAlg::ConvertROPtoWirePlane()
  
  
  //----------------------------------------------------------------------------
  
} // namespace
