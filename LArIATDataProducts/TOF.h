////////////////////////////////////////////////////////////////////////
// $Id: TOF.h,v 1.00 2015/06/03 16:04:20 dsmith Exp $
//
// dansmith@bu.edu
//
////////////////////////////////////////////////////////////////////////

#ifndef LARIATDATAPRODUCTS_TOF_H
#define LARIATDATAPRODUCTS_TOF_H

#include <vector>
#include <iosfwd>
#include <string>

///Raw data description
namespace ldp {
  
  class TOF {

  public:
    TOF();    // Default constructor
    ~TOF(){}; // Default destructor
    
  private:
    std::vector<float> fTOF;
    std::vector<long>  fTimeStamp;
    
#ifndef __GCCXML__

  public:

    //  TOF( std::vector<short> TOF, std::vector<long> TimeStamp ):
    //fTOF(TOF), fTimeStamp(TimeStamp) {}
    
    TOF( std::vector<float> TOF, std::vector<long> TimeStamp ) { 
      fTOF = TOF;
      fTimeStamp = TimeStamp;
    }

    // Get Methods
    float  SingleTOF(size_t iHit) const;
    long   TimeStamp(size_t iHit) const;
    // Set Methods
    void   SetSingleTOF(std::vector<float> value) {fTOF       = value;};
    void   SetTimeStamp(std::vector<long>  value) {fTimeStamp = value;};
    
    size_t NTOF() const;

#endif
  };
}

#ifndef __GCCXML__

inline size_t ldp::TOF::NTOF()    const { return fTOF.size(); }

#endif

#endif // LARIATDATAPRODUCTS_TOF_H

////////////////////////////////////////////////////////////////////////
