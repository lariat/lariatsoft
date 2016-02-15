#include "art/Persistency/Common/Wrapper.h"
#include "art/Persistency/Common/Assns.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"

#include "LArIATFragments/V1495Fragment.h"

#include "RawData/RawDigit.h"
#include "RawData/AuxDetDigit.h"
#include "RawData/OpDetPulse.h"
#include "RawData/TriggerData.h"
#include "RecoBase/Wire.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Vertex.h"
#include "RecoBase/EndPoint2D.h"
#include "RecoBase/SpacePoint.h"
#include "RecoBase/Track.h"
#include "LArIATDataProducts/WCTrack.h"
#include "LArIATDataProducts/Stopping.h"

#include "RawData/RawDigit.h"
#include "RawData/AuxDetDigit.h"
#include "RawData/OpDetPulse.h"
#include "RawData/TriggerData.h"
#include "SummaryData/RunData.h"

/* template class std::vector<V1495Fragment>; */
/* template class std::vector<V1495ChannelData>; */
/* template class std::vector<V1495TriggerPatternData>; */
/* template class art::Wrapper< std::vector<V1495Fragment>           >; */
/* template class art::Wrapper< std::vector<V1495ChannelData>	  >; */
/* template class art::Wrapper< std::vector<V1495TriggerPatternData> >; */

template class std::pair<art::Ptr<raw::Trigger>, art::Ptr<raw::RawDigit>     >;
template class std::pair<art::Ptr<raw::Trigger>, art::Ptr<recob::Wire>       >;
template class std::pair<art::Ptr<raw::Trigger>, art::Ptr<recob::Hit>        >;
template class std::pair<art::Ptr<raw::Trigger>, art::Ptr<recob::Cluster>    >;
template class std::pair<art::Ptr<raw::Trigger>, art::Ptr<recob::EndPoint2D> >;
template class std::pair<art::Ptr<raw::Trigger>, art::Ptr<recob::SpacePoint> >;
template class std::pair<art::Ptr<raw::Trigger>, art::Ptr<recob::Track>      >;
template class std::pair<art::Ptr<raw::Trigger>, art::Ptr<recob::Vertex>     >;
template class std::pair<art::Ptr<raw::Trigger>, art::Ptr<raw::AuxDetDigit>  >;
template class std::pair<art::Ptr<raw::Trigger>, art::Ptr<raw::OpDetPulse>   >;
//template class art::pair<art::Ptr<ldp::WCTrack>, art::Ptr<recob::Track>      >;


template class std::pair<art::Ptr<raw::RawDigit>,     art::Ptr<raw::Trigger> >;
template class std::pair<art::Ptr<recob::Wire>,       art::Ptr<raw::Trigger> >;
template class std::pair<art::Ptr<recob::Hit>,        art::Ptr<raw::Trigger> >;
template class std::pair<art::Ptr<recob::Cluster>,    art::Ptr<raw::Trigger> >;
template class std::pair<art::Ptr<recob::EndPoint2D>, art::Ptr<raw::Trigger> >;
template class std::pair<art::Ptr<recob::SpacePoint>, art::Ptr<raw::Trigger> >;
template class std::pair<art::Ptr<recob::Track>,      art::Ptr<raw::Trigger> >;
template class std::pair<art::Ptr<recob::Vertex>,     art::Ptr<raw::Trigger> >;
template class std::pair<art::Ptr<raw::AuxDetDigit>,  art::Ptr<raw::Trigger> >;
template class std::pair<art::Ptr<raw::OpDetPulse>,   art::Ptr<raw::Trigger> >;
//template class art::pair<art::Ptr<recob::Track>,      art::Ptr<ldp::WCTrack> >;


template class art::Assns<raw::Trigger, raw::RawDigit,     void>;
template class art::Assns<raw::Trigger, recob::Wire,       void>;
template class art::Assns<raw::Trigger, recob::Hit,        void>;
template class art::Assns<raw::Trigger, recob::Cluster,    void>;
template class art::Assns<raw::Trigger, recob::Vertex,     void>;
template class art::Assns<raw::Trigger, recob::EndPoint2D, void>;
template class art::Assns<raw::Trigger, recob::SpacePoint, void>;
template class art::Assns<raw::Trigger, recob::Track,      void>;
template class art::Assns<raw::Trigger, raw::AuxDetDigit,  void>;
template class art::Assns<raw::Trigger, raw::OpDetPulse,   void>;
template class art::Assns<ldp::WCTrack, recob::Track,      void>;
template class art::Assns<recob::Track, recob::Stopping,   void>; 

template class art::Wrapper<art::Assns<raw::Trigger, raw::RawDigit,     void> >;
template class art::Wrapper<art::Assns<raw::Trigger, recob::Wire,       void> >;
template class art::Wrapper<art::Assns<raw::Trigger, recob::Hit,        void> >;
template class art::Wrapper<art::Assns<raw::Trigger, recob::Cluster,    void> >;
template class art::Wrapper<art::Assns<raw::Trigger, recob::EndPoint2D, void> >;
template class art::Wrapper<art::Assns<raw::Trigger, recob::SpacePoint, void> >;
template class art::Wrapper<art::Assns<raw::Trigger, recob::Track,      void> >;
template class art::Wrapper<art::Assns<raw::Trigger, recob::Vertex,     void> >;
template class art::Wrapper<art::Assns<raw::Trigger, raw::AuxDetDigit,  void> >;
template class art::Wrapper<art::Assns<raw::Trigger, raw::OpDetPulse,   void> >;
template class art::Wrapper<art::Assns<ldp::WCTrack, recob::Track,      void> >;
