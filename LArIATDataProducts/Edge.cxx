/** ****************************************************************************
 * @file Edge.cxx
 * @brief Definition of basic channel signal object.
 * @author brebel@fnal.gov
 * @see  Edge.h
 * 
 * ****************************************************************************/

#include "LArIATDataProducts/Edge.h"

// C/C++ standard libraries
#include <utility> // std::move()

namespace ldp{

  //----------------------------------------------------------------------
  Edge::Edge()
    : fChannel(raw::InvalidChannelID)
    , fView(geo::kUnknown)
    , fSignalROI()
    {}

  //----------------------------------------------------------------------
  Edge::Edge(
    RegionsOfInterest_t const& sigROIlist,
    raw::ChannelID_t channel,
    geo::View_t view
    )
    : fChannel(channel)
    , fView(view)
    , fSignalROI(sigROIlist)
    {}

  //----------------------------------------------------------------------
  Edge::Edge(
    RegionsOfInterest_t&& sigROIlist,
    raw::ChannelID_t channel,
    geo::View_t view
    )
    : fChannel(channel)
    , fView(view)
    , fSignalROI(std::move(sigROIlist))
    {}


  //----------------------------------------------------------------------
  std::vector<float> Edge::Signal() const {
    return { fSignalROI.begin(), fSignalROI.end() };
  } // Edge::Signal()


}
////////////////////////////////////////////////////////////////////////

