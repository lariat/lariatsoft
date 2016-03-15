
#include "SpillWrapper.h"

#include "artdaq-core/Data/Fragments.hh"

#include "cetlib/exception.h"

#include <iostream>

namespace rdu {

  void SpillWrapper::add(const artdaq::Fragment& frag) {

    if (nfragments_processed_ == 0) {

      // Basically, it's the LAST of the fragments which is the
      // largest if the spill couldn't be evenly split into four
      // fragments, so that's why we add the "+1" padding to the
      // buffer

      size_t bufsize = (frag.dataSize() + 1) * sizeof(artdaq::RawDataType) * nfragments_ ;

      std::cout << "Fragment size in RawDataTypes is " << frag.dataSize() << std::endl;
      std::cout << "buffer size is " << bufsize  << std::endl;

      buffer_.resize( bufsize );
    }

    // JCF, Mar-1-2016: think carefully about what you're doing with
    // this memcpy command...

    std::cout << "About to copy " << frag.dataSizeBytes() << " bytes from " <<
      static_cast<const void*>(frag.dataBeginBytes()) << 
      " to " <<
      static_cast<const void*>(&*buffer_.begin() + size_) << 
      std::endl;

    memcpy( &*buffer_.begin() + size_, frag.dataBeginBytes(), frag.dataSizeBytes() );

    // Bookkeeping
    size_ += frag.dataSizeBytes();
    nfragments_processed_++;
  }

  void SpillWrapper::add(const art::Event& evt) {

    art::ValidHandle<artdaq::Fragments> raw = 
      evt.getValidHandle<artdaq::Fragments>("daq:SPILL");

    if ( raw->size() != 1) {
      throw cet::exception("SpillWrapper") << "Error in SpillWrapper::add(): expected 1 artdaq fragment in event, found " 
					   << raw->size() << "\n";
    }

    add( (*raw)[0] );
  }

}
