#ifndef SPILLWRAPPER_H
#define SPILLWRAPPER_H

#include "art/Framework/Principal/Event.h"
#include "artdaq-core/Data/Fragment.hh"

#include <vector>

// JCF, Mar-11-2016

// SpillWrapper serves as a user interface to the raw binary
// representation of a spill which has been broken up into art events
// of one fragment each. All it needs to know is the number of events
// (i.e., fragments) per spill; this value is passed to its
// constructor.

// If we have an instance of SpillWrapper called "sw", to access a
// pointer (const uint8_t*) to the underlying spill, call
// sw->get(). Note, however, that you're not guaranteed to have the
// full spill represented in memory unless sw->ready() == true. This
// is because the internal representation of the spill inside
// SpillWrapper is built up through calls to its add() functions; this
// overloaded function can either directly take an artdaq::Fragment
// object or an art event which contains an artdaq::Fragment. To
// determine the size of SpillWrapper's representation of the spill,
// call sw->size(); this will only be the size of the entire spill if
// sw->ready() == true.

namespace rdu {

  class SpillWrapper {

  public:
    SpillWrapper(std::size_t nfragments) :
      nfragments_(nfragments),
      nfragments_processed_(0),
      size_(0)
    {}

    void add(const art::Event& );

    void add(const artdaq::Fragment& );

    // Consider having get() throw an exception if ready() != true
    const uint8_t* get() const {return &*buffer_.begin(); } 

    std::size_t size() const {return size_; }

    bool ready() const { return (nfragments_processed_ == nfragments_) ? true : false; }

  private:

    const std::size_t nfragments_;
    std::vector<uint8_t> buffer_;

    std::size_t nfragments_processed_;
    std::size_t size_;


  };

}

#endif
