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



#include "RawData/RawDigit.h"
#include "RawData/OpDetPulse.h"
#include "RawData/AuxDetDigit.h"
#include "RawData/BeamInfo.h"
#include "RawData/ExternalTrigger.h"
#include "RawData/TriggerData.h"
#include "RawData/OpDetWaveform.h"

//
// Only include objects that we would like to be able to put into the event.
// Do not include the objects they contain internally.
//
template class std::vector<ldp::WCTrack>;

template class art::Wrapper< ldp::WCTrack                  >;

template class art::Wrapper< std::vector<ldp::WCTrack >    >;

template class std::pair< art::Ptr<raw::Trigger>,     art::Ptr<ldp::WCTrack> >    ;
template class std::pair< art::Ptr<raw::AuxDetDigit>, art::Ptr<ldp::WCTrack> >    ;

template class std::pair< art::Ptr<ldp::WCTrack>,     art::Ptr<raw::Trigger> >    ;
template class std::pair< art::Ptr<ldp::WCTrack>,     art::Ptr<raw::AuxDetDigit> >;

template class art::Assns<raw::Trigger,     ldp::WCTrack,      void> >; 
template class art::Assns<raw::AuxDetDigit, ldp::WCTrack,      void> >;

template class art::Wrapper<art::Assns<raw::Trigger,     ldp::WCTrack,     void> >;
template class art::Wrapper<art::Assns<raw::AuxDetDigit, ldp::WCTrack,     void> >;
