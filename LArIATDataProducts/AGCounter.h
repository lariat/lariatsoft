////////////////////////////////////////////////////////////////
//                                                            //
// Data Product for the Aerogel Cherenkov Counter 	      //
//                                                            //
// Authors: UT Austin Karol Lang Group			      //
//							      //
//         *Dung Phan (brianp.dung@gmail.com)		      //
//	    Will Flanagan (will.flanagan@utexas.edu)	      //
//	    Brandon Soubasis (brandon.soubasis@gmail.com      //
//                                                            //
////////////////////////////////////////////////////////////////

#ifndef LARIATDATAPRODUCTS_AGCOUNTER_H
#define LARIATDATAPRODUCTS_AGCOUNTER_H

#include <vector>
#include <string>

namespace ldp {

  struct AGCHits {
    long unsigned int 	TriggerTimeStamp;
    
    long unsigned int 	HitTimeStampUSE;
    long unsigned int 	HitTimeStampUSW;
    long unsigned int 	HitTimeStampDS1;
    long unsigned int 	HitTimeStampDS2;
    
    float			HitPulseAreaUSE;
    float			HitPulseAreaUSW;
    float			HitPulseAreaDS1;
    float			HitPulseAreaDS2;
    
    bool 			HitExistUSE;
    bool 			HitExistUSW;
    bool 			HitExistDS1;
    bool 			HitExistDS2;
  };

  class AGCounter {
  public:
    AGCounter();
    ~AGCounter();

#ifndef __GCCXML__
    AGCounter(std::vector<ldp::AGCHits> const& AGCHits);

    size_t GetNHits() const;
    
    long int		GetTriggerTimeStamp(int);
    
    long unsigned int 	GetHitTimeStampUSE(size_t) const;
    long unsigned int 	GetHitTimeStampUSW(size_t) const;
    long unsigned int 	GetHitTimeStampDS1(size_t) const;
    long unsigned int 	GetHitTimeStampDS2(size_t) const;
    
    float			GetHitPulseAreaUSE(size_t) const;
    float			GetHitPulseAreaUSW(size_t) const;
    float			GetHitPulseAreaDS1(size_t) const;
    float			GetHitPulseAreaDS2(size_t) const;
    
    bool 			GetHitExistUSE(size_t) const;
    bool 			GetHitExistUSW(size_t) const;
    bool 			GetHitExistDS1(size_t) const;
    bool 			GetHitExistDS2(size_t) const;
#endif
    
  private:
  	std::vector<ldp::AGCHits> fAGCHits;
  };
}

#ifndef __GCCXML__

inline size_t	            ldp::AGCounter::GetNHits()                      const { return fAGCHits.size(); }
inline long unsigned int 	ldp::AGCounter::GetHitTimeStampUSE(size_t iHit) const { return fAGCHits[iHit].HitTimeStampUSE; }
inline long unsigned int 	ldp::AGCounter::GetHitTimeStampUSW(size_t iHit) const { return fAGCHits[iHit].HitTimeStampUSW; }
inline long unsigned int 	ldp::AGCounter::GetHitTimeStampDS1(size_t iHit) const { return fAGCHits[iHit].HitTimeStampDS1; }
inline long unsigned int 	ldp::AGCounter::GetHitTimeStampDS2(size_t iHit) const { return fAGCHits[iHit].HitTimeStampDS2; }

inline float			        ldp::AGCounter::GetHitPulseAreaUSE(size_t iHit) const { return fAGCHits[iHit].HitPulseAreaUSE;}
inline float			        ldp::AGCounter::GetHitPulseAreaUSW(size_t iHit) const { return fAGCHits[iHit].HitPulseAreaUSW;}
inline float			        ldp::AGCounter::GetHitPulseAreaDS1(size_t iHit) const { return fAGCHits[iHit].HitPulseAreaDS1;}
inline float			        ldp::AGCounter::GetHitPulseAreaDS2(size_t iHit) const { return fAGCHits[iHit].HitPulseAreaDS2;}

inline bool 			        ldp::AGCounter::GetHitExistUSE(size_t iHit)     const { return fAGCHits[iHit].HitExistUSE; }
inline bool 			        ldp::AGCounter::GetHitExistUSW(size_t iHit)     const { return fAGCHits[iHit].HitExistUSW; }
inline bool 			        ldp::AGCounter::GetHitExistDS1(size_t iHit)     const { return fAGCHits[iHit].HitExistDS1; }
inline bool 			        ldp::AGCounter::GetHitExistDS2(size_t iHit)     const { return fAGCHits[iHit].HitExistDS2; }

#endif // __GCCXML__

#endif // LARIATDATAPRODUCTS_AGCOUNTER_H
