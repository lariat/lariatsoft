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
    
    long unsigned int 	HitTimeStamp1p10_1;
    long unsigned int 	HitTimeStamp1p10_2;
    long unsigned int 	HitTimeStamp1p06_1;
    long unsigned int 	HitTimeStamp1p06_2;
    
    float			HitPulseArea1p10_1;
    float			HitPulseArea1p10_2;
    float			HitPulseArea1p06_1;
    float			HitPulseArea1p06_2;
    
    bool 			HitExist1p10_1;
    bool 			HitExist1p10_2;
    bool 			HitExist1p06_1;
    bool 			HitExist1p06_2;
  };

  class AGCounter {
  public:
    AGCounter();
    ~AGCounter();

#ifndef __GCCXML__
    AGCounter(std::vector<ldp::AGCHits> const& AGCHits);

    size_t GetNHits() const;
    
    long int		GetTriggerTimeStamp(int);
    
    long unsigned int 	GetHitTimeStamp1p10_1(size_t) const;
    long unsigned int 	GetHitTimeStamp1p10_2(size_t) const;
    long unsigned int 	GetHitTimeStamp1p06_1(size_t) const;
    long unsigned int 	GetHitTimeStamp1p06_2(size_t) const;
    
    float			GetHitPulseArea1p10_1(size_t) const;
    float			GetHitPulseArea1p10_2(size_t) const;
    float			GetHitPulseArea1p06_1(size_t) const;
    float			GetHitPulseArea1p06_2(size_t) const;
    
    bool 			GetHitExist1p10_1(size_t) const;
    bool 			GetHitExist1p10_2(size_t) const;
    bool 			GetHitExist1p06_1(size_t) const;
    bool 			GetHitExist1p06_2(size_t) const;
#endif
    
  private:
  	std::vector<ldp::AGCHits> fAGCHits;
  };
}

#ifndef __GCCXML__

inline size_t	            ldp::AGCounter::GetNHits()                      const { return fAGCHits.size(); }
inline long unsigned int 	ldp::AGCounter::GetHitTimeStamp1p10_1(size_t iHit) const { return fAGCHits[iHit].HitTimeStamp1p10_1; }
inline long unsigned int 	ldp::AGCounter::GetHitTimeStamp1p10_2(size_t iHit) const { return fAGCHits[iHit].HitTimeStamp1p10_2; }
inline long unsigned int 	ldp::AGCounter::GetHitTimeStamp1p06_1(size_t iHit) const { return fAGCHits[iHit].HitTimeStamp1p06_1; }
inline long unsigned int 	ldp::AGCounter::GetHitTimeStamp1p06_2(size_t iHit) const { return fAGCHits[iHit].HitTimeStamp1p06_2; }

inline float			        ldp::AGCounter::GetHitPulseArea1p10_1(size_t iHit) const { return fAGCHits[iHit].HitPulseArea1p10_1;}
inline float			        ldp::AGCounter::GetHitPulseArea1p10_2(size_t iHit) const { return fAGCHits[iHit].HitPulseArea1p10_2;}
inline float			        ldp::AGCounter::GetHitPulseArea1p06_1(size_t iHit) const { return fAGCHits[iHit].HitPulseArea1p06_1;}
inline float			        ldp::AGCounter::GetHitPulseArea1p06_2(size_t iHit) const { return fAGCHits[iHit].HitPulseArea1p06_2;}

inline bool 			        ldp::AGCounter::GetHitExist1p10_1(size_t iHit)     const { return fAGCHits[iHit].HitExist1p10_1; }
inline bool 			        ldp::AGCounter::GetHitExist1p10_2(size_t iHit)     const { return fAGCHits[iHit].HitExist1p10_2; }
inline bool 			        ldp::AGCounter::GetHitExist1p06_1(size_t iHit)     const { return fAGCHits[iHit].HitExist1p06_1; }
inline bool 			        ldp::AGCounter::GetHitExist1p06_2(size_t iHit)     const { return fAGCHits[iHit].HitExist1p06_2; }

#endif // __GCCXML__

#endif // LARIATDATAPRODUCTS_AGCOUNTER_H
