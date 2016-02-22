//
// Build a dictionary.
//
// $Id: classes.h,v 1.8 2010/04/12 18:12:28  Exp $
// $Author:  $
// $Date: 2010/04/12 18:12:28 $
// 
// Original author Rob Kutschke, modified by klg
//
// Notes:
// 1) The system is not able to deal with
//    art::Wrapper<std::vector<std::string> >;
//    The problem is somewhere inside root's reflex mechanism
//    and Philippe Canal says that it is ( as of March 2010) a
//    known problem.  He also says that they do not have any
//    plans to fix it soon.  We can always work around it 
//    by putting the string inside another object.

#include "art/Persistency/Common/Wrapper.h"
#include "art/Persistency/Common/Assns.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"



#include "lardata/RawData/RawDigit.h"
#include "lardata/RawData/OpDetPulse.h"
#include "lardata/RawData/AuxDetDigit.h"
#include "lardata/RawData/BeamInfo.h"
#include "lardata/RawData/ExternalTrigger.h"
#include "lardata/RawData/TriggerData.h"
#include "lardata/RawData/OpDetWaveform.h"
#include "LArIATDataProducts/AGCounter.h"
#include "LArIATDataProducts/AuxDetParticleID.h"
#include "LArIATDataProducts/ConditionsSummary.h"
#include "LArIATDataProducts/MuonRangeStackHits.h"
#include "LArIATDataProducts/TOF.h"
#include "LArIATDataProducts/WCTrack.h"
#include "LArIATDataProducts/TOF.h"
#include "LArIATDataProducts/WCTrack.h"
#include "LArIATDataProducts/Edge.h"


//
// Only include objects that we would like to be able to put into the event.
// Do not include the objects they contain internally.
//

//Stuff regarding AuxDetParticleID
template class std::vector<ldp::AuxDetParticleID>;
template class art::Wrapper<ldp::AuxDetParticleID>;
template class art::Wrapper<std::vector<ldp::AuxDetParticleID> >;

template class std::vector<ldp::TOF>;

template class art::Wrapper< ldp::TOF                  >;

template class art::Wrapper< std::vector<ldp::TOF >    >;

template class std::vector<ldp::AGCounter>;

template class art::Wrapper< ldp::AGCounter                  >;

template class art::Wrapper< std::vector<ldp::AGCounter >    >;

template class  std::vector<AGCHits>;

template class art::Wrapper< AGCHits                  >;

template class art::Wrapper< std::vector<AGCHits >    >;

template class std::pair< art::Ptr<raw::Trigger>,     art::Ptr<ldp::TOF> >    ;
template class std::pair< art::Ptr<raw::AuxDetDigit>, art::Ptr<ldp::TOF> >    ;

template class std::pair< art::Ptr<ldp::TOF>,     art::Ptr<raw::Trigger> >    ;
template class std::pair< art::Ptr<ldp::TOF>,     art::Ptr<raw::AuxDetDigit> >;

template class art::Assns<ldp::TOF,     raw::Trigger,      void>; 
template class art::Assns<raw::Trigger,     ldp::TOF,      void>; 
template class art::Assns<raw::AuxDetDigit, ldp::TOF,      void>;

template class art::Wrapper<art::Assns<raw::Trigger,     ldp::TOF,     void> >;
template class art::Wrapper<art::Assns<raw::AuxDetDigit, ldp::TOF,     void> >;
template class art::Wrapper<art::Assns<ldp::TOF,     raw::Trigger,     void> >;

template class std::vector<ldp::WCTrack>;

template class art::Wrapper< ldp::WCTrack                  >;

template class art::Wrapper< std::vector<ldp::WCTrack >    >;

template class std::pair< art::Ptr<raw::Trigger>,     art::Ptr<ldp::WCTrack> >    ;
template class std::pair< art::Ptr<raw::AuxDetDigit>, art::Ptr<ldp::WCTrack> >    ;

template class std::pair< art::Ptr<ldp::WCTrack>,     art::Ptr<raw::Trigger> >    ;
template class std::pair< art::Ptr<ldp::WCTrack>,     art::Ptr<raw::AuxDetDigit> >;

template class art::Assns<ldp::WCTrack,     raw::Trigger,      void>; 
template class art::Assns<raw::Trigger,     ldp::WCTrack,      void>; 
template class art::Assns<raw::AuxDetDigit, ldp::WCTrack,      void>;

template class art::Wrapper<art::Assns<raw::Trigger,     ldp::WCTrack,     void> >;
template class art::Wrapper<art::Assns<raw::AuxDetDigit, ldp::WCTrack,     void> >;
template class art::Wrapper<art::Assns<ldp::WCTrack,     raw::Trigger,     void> >;



//Stuff regarding/belonging to Edge Rec
template class std::vector<ldp::Edge>;
template class art::Wrapper< ldp::Edge                  >;
template class art::Wrapper< std::vector<ldp::Edge >    >;

template class std::pair< art::Ptr<raw::Trigger>,     art::Ptr<ldp::Edge> >    ;
template class std::pair< art::Ptr<raw::AuxDetDigit>, art::Ptr<ldp::Edge> >    ;

template class std::pair< art::Ptr<ldp::Edge>,     art::Ptr<raw::Trigger> >    ;
template class std::pair< art::Ptr<ldp::Edge>,     art::Ptr<raw::AuxDetDigit> >;

template class art::Assns<ldp::Edge,     raw::Trigger,      void>; 
template class art::Assns<raw::Trigger,     ldp::Edge,      void>; 
template class art::Assns<raw::AuxDetDigit, ldp::Edge,      void>;

template class art::Wrapper<art::Assns<raw::Trigger,     ldp::Edge,     void> >;
template class art::Wrapper<art::Assns<raw::AuxDetDigit, ldp::Edge,     void> >;
template class art::Wrapper<art::Assns<ldp::Edge,     raw::Trigger,     void> >;



//Stuff regarding/belonging to Muon Range Stack Reco
template class std::vector<ldp::MuonRangeStackHits >;
template class art::Wrapper< ldp::MuonRangeStackHits               >;
template class art::Wrapper< std::vector<ldp::MuonRangeStackHits >    >;
template class std::map<int,std::vector<int> >;
template class art::Wrapper<std::map<int,std::vector<int> > >;
template class art::Wrapper<MuRSTrack>;
template class std::vector<MuRSTrack>;
template class art::Wrapper<std::vector<MuRSTrack> >;
template class std::pair<int,int>;
template class art::Wrapper<std::pair<int,int> >;
template class std::vector<std::pair<int,int> >;
template class art::Wrapper<std::vector<std::pair<int,int> > >;


template class art::Wrapper<ldp::ConditionsSummary>;
