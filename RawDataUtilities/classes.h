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

/* template class std::vector<V1495Fragment>; */
/* template class std::vector<V1495ChannelData>; */
/* template class std::vector<V1495TriggerPatternData>; */
/* template class art::Wrapper< std::vector<V1495Fragment>           >; */
/* template class art::Wrapper< std::vector<V1495ChannelData>	  >; */
/* template class art::Wrapper< std::vector<V1495TriggerPatternData> >; */

template class std::pair<art::Ptr<raw::Trigger>, art::Ptr<raw::RawDigit>    >;
template class std::pair<art::Ptr<raw::Trigger>, art::Ptr<recob::Wire>      >;
template class std::pair<art::Ptr<raw::Trigger>, art::Ptr<recob::Hit>       >;
template class std::pair<art::Ptr<raw::Trigger>, art::Ptr<recob::Cluster>       >;
template class std::pair<art::Ptr<raw::Trigger>, art::Ptr<std::vector<recob::Cluster> >       >;
template class std::pair<art::Ptr<raw::Trigger>, art::Ptr<raw::AuxDetDigit> >;
template class std::pair<art::Ptr<raw::Trigger>, art::Ptr<raw::OpDetPulse>  >;

template class std::pair<art::Ptr<raw::RawDigit>,    art::Ptr<raw::Trigger> >;
template class std::pair<art::Ptr<recob::Wire>,      art::Ptr<raw::Trigger> >;
template class std::pair<art::Ptr<recob::Hit>,       art::Ptr<raw::Trigger> >;
template class std::pair<art::Ptr<recob::Cluster>,   art::Ptr<raw::Trigger> >;
template class std::pair<art::Ptr<std::vector<recob::Cluster> >,   art::Ptr<raw::Trigger> >;
template class std::pair<art::Ptr<raw::AuxDetDigit>, art::Ptr<raw::Trigger> >;
template class std::pair<art::Ptr<raw::OpDetPulse>,  art::Ptr<raw::Trigger> >;

template class art::Assns<raw::Trigger, raw::RawDigit,    void>;
template class art::Assns<raw::Trigger, recob::Wire, void>;
template class art::Assns<raw::Trigger, recob::Hit, void>;
template class art::Assns<raw::Trigger, recob::Cluster, void>;		
template class art::Assns<raw::Trigger, std::vector<recob::Cluster>, void>;		
template class art::Assns<raw::Trigger, raw::AuxDetDigit, void>;
template class art::Assns<raw::Trigger, raw::OpDetPulse,  void>;

template class art::Wrapper<art::Assns<raw::Trigger, raw::RawDigit,    void> >;
template class art::Wrapper<art::Assns<raw::Trigger, recob::Wire,      void> >;
template class art::Wrapper<art::Assns<raw::Trigger, recob::Hit,       void> >;
template class art::Wrapper<art::Assns<raw::Trigger, recob::Cluster,       void> >;
template class art::Wrapper<art::Assns<raw::Trigger, std::vector<recob::Cluster>,       void> >;
template class art::Wrapper<art::Assns<raw::Trigger, raw::AuxDetDigit, void> >;
template class art::Wrapper<art::Assns<raw::Trigger, raw::OpDetPulse,  void> >;
