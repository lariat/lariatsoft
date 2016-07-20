/** ****************************************************************************
 * @file Edge.h
 * @brief Declaration of basic channel signal object.
 * @author brebel@fnal.gov
 * @see  Edge.cxx
 * 
 * Changes:
 * 20141211 Gianluca Petrillo (petrillo@fnal.gov)
 *   data architecture revision changes:
 *   - fSignalType and SignalType() removed
 *   - fRawDigit and RawDigit() removed
 *   - fChannel and Channel() added
 *   - constructors now take directly a RawDigit, not its art::Ptr
 * 
 * ****************************************************************************/

#ifndef EDGE_H
#define EDGE_H

#include <vector>

#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataobj/Utilities/sparse_vector.h"


// prototype declaration
namespace raw { class RawDigit; }


/// Reconstruction base classes
namespace ldp {

  /**
   * @brief Class holding the deconvoluted signals from a channel
   * @note Despite the name, this class is associated to a readout channel, not
   *       just a wire
   *
   * The channel content is expected to have been filtered from noise and
   * corrected for electronics response.
   * The content is presented as calibrated ADC counts, pedestal removed, as
   * function of time in discrete TDC units (use TimeService to discover the
   * exact extent of each tick).
   * The content is organized as time intervals where some signal is present
   * ("regions of interest", RoI), outside which we assume no signal, i.e,
   * calibrated ADC counts of 0.
   * Strictly speaking, the definition of the regions of interest is a negative
   * one: we first define where we are sure no signal is present; the rest will
   * constitute regions of interest.
   * The identification of such regions is responsibility of the algorithm
   * creating the Wire object. In the simple approach, the whole readout window
   * is stored in a single region of interest, meaning that we don't claim any
   * of the channel signal to be definitely signal free.
   * 
   * Algorithms using the regions of interest can access the channel signal
   * information either ignoring the regions of interest, and being potentially
   * flooded by zeroes from the non-signal regions:
   *     
   *     for (float ADCcount: wire.Signal()) ...
   *     
   * or they can analyze region by region:
   *     
   *     for (auto iROI = wire.begin_range(); wire.end_range(); ++iROI) {
   *       const datarange_t& ROI = *iROI;
   *       const int FirstTick = ROI.begin_index();
   *       const int EndTick = ROI.end_index();
   *       const float FirstADC = ROI[FirstTick]; // index access is by absolute tick number
   *       for (float ADC: ROI) // ... or you can just iterate through
   *       // ...
   *     } // for
   *     
   * Note that the indexed access is always by absolute tick number.
   * More examples of the use of SignalROI() return value are documented in
   * lar::sparse_vector .
   * 
   * Each channel is associated with a raw::RawDigit object. These associations
   * should be stored together with recob::Wire by the producer in a art::Assns
   * data product.
   */
  class Edge {
    public:
      /// a region of interest is a pair (TDC offset, readings)
      typedef lar::sparse_vector<float> RegionsOfInterest_t;
      
      /// Default constructor: a wire with no signal information
      Edge();
      
    private:
      raw::ChannelID_t    fChannel;   ///< ID of the associated channel
      geo::View_t         fView;      ///< View corresponding to the plane of this wire
      RegionsOfInterest_t fSignalROI; ///< Signal on the channel

#ifndef __GCCXML__
    
    public:
      /**
       * @brief Constructor: uses specified signal in regions of interest
       * @param sigROIlist signal organized in regions of interest
       * @param channel the ID of the channel
       * @param view the view the channel belongs to
       *
       * The information used from the raw digit are the channel ID and the
       * length in samples (TDC ticks) of the original readout window.
       * 
       * Signal is copied into the Edge object. If possible, use the other
       * constructor that moves the data instead.
       */
      Edge(
        RegionsOfInterest_t const& sigROIlist,
        raw::ChannelID_t channel,
        geo::View_t view
        );
      
      /**
       * @brief Constructor: uses specified signal in regions of interest
       * @param sigROIlist signal organized in regions of interest
       * @param channel the ID of the channel
       * @param view the view the channel belongs to
       *
       * The information used from the raw digit are the channel ID and the
       * length in samples (TDC ticks) of the original readout window.
       * 
       * Signal information is moved from sigROIlist, that becomes empty.
       */
      Edge(
        RegionsOfInterest_t&& sigROIlist,
        raw::ChannelID_t channel,
        geo::View_t view
        );
      
      ///@{
      ///@name Accessors
      
      /// Return a zero-padded full length vector filled with RoI signal
      std::vector<float>  Signal() const;
      
      /// Returns the list of regions of interest
      const RegionsOfInterest_t& SignalROI()  const;
      
      /// Returns the number of time ticks, or samples, in the channel
      size_t                     NSignal()    const;
      
      /// Returns the view the channel belongs to
      geo::View_t                View()       const;
      
      /// Returns the ID of the channel (or InvalidChannelID)
      raw::ChannelID_t           Channel()    const;
      ///@}
      
      /// Returns whether this channel ID is smaller than the other
      bool operator< (const Edge& than) const;
      
#endif // !__GCCXML__
      
  }; // class Edge
  
} // namespace recob

#ifndef __GCCXML__

inline const ldp::Edge::RegionsOfInterest_t&
                                  ldp::Edge::SignalROI()  const { return fSignalROI;        }
inline size_t                     ldp::Edge::NSignal()    const { return fSignalROI.size(); }
inline geo::View_t                ldp::Edge::View()       const { return fView;             }
inline raw::ChannelID_t           ldp::Edge::Channel()    const { return fChannel;          }
inline bool                       ldp::Edge::operator< (const Edge& than) const
  { return Channel() < than.Channel(); }

#endif // !__GCCXML__

#endif // WIRE_H

////////////////////////////////////////////////////////////////////////
