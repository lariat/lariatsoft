#include "canvas/Persistency/Common/Wrapper.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"

#include "LArIATFragments/V1495Fragment.h"


#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/AuxDetDigit.h"
#include "lardataobj/RawData/OpDetPulse.h"
#include "lardataobj/RawData/TriggerData.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"


#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/AuxDetDigit.h"
#include "lardataobj/RawData/OpDetPulse.h"
#include "lardataobj/RawData/TriggerData.h"

#include "larcoreobj/SummaryData/RunData.h"

#include "LArIATDataProducts/WCTrack.h"

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
template class std::pair<art::Ptr<ldp::WCTrack>, art::Ptr<recob::Track>      >;
template class std::pair<art::Ptr<recob::Track>, art::Ptr<ldp::WCTrack>      >;


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
//template class std::pair<art::Ptr<recob::Track>,      art::Ptr<ldp::WCTrack> >;
//template class std::pair<art::Ptr<ldp::WCTrack>,      art::Ptr<recob::Track> >;


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
template class art::Assns<recob::Track, ldp::WCTrack,      void>;


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
template class art::Wrapper<art::Assns<recob::Track, ldp::WCTrack,      void> >;

template class art::Wrapper<art::PtrVector<recob::Track> >;

