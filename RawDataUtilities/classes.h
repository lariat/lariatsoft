#include "art/Persistency/Common/Wrapper.h"
#include "art/Persistency/Common/Assns.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"

#include "LArIATFragments/V1495Fragment.h"

#include "lardata/RawData/RawDigit.h"
#include "lardata/RawData/AuxDetDigit.h"
#include "lardata/RawData/OpDetPulse.h"
#include "lardata/RawData/TriggerData.h"
#include "lardata/RecoBase/Wire.h"
#include "lardata/RecoBase/Hit.h"
#include "lardata/RecoBase/Cluster.h"
#include "lardata/RecoBase/Vertex.h"
#include "lardata/RecoBase/EndPoint2D.h"
#include "lardata/RecoBase/SpacePoint.h"
#include "lardata/RecoBase/Track.h"

#include "lardata/RawData/RawDigit.h"
#include "lardata/RawData/AuxDetDigit.h"
#include "lardata/RawData/OpDetPulse.h"
#include "lardata/RawData/TriggerData.h"
#include "larcore/SummaryData/RunData.h"

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
