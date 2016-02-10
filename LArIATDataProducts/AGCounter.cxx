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

#include "LArIATDataProducts/AGCounter.h"
#include "cetlib/exception.h"

int iHit;
namespace ldp {
  
  AGCounter::AGCounter() {}
  AGCounter::~AGCounter() {}
  
  //	AGCounter::AGCounter(std::vector<std::vector<AGCHits> > AGCHits) { fAGCHits = AGCHits; }
  
  AGCounter::AGCounter(std::vector<AGCHits> AGCHits) { fAGCHits = AGCHits; }
  
  /*	void AGCounter::Linearize() {
	fAGCHitsLinearized.clear();
	for (size_t i = 0; i < fAGCHits.size(); i++) {
	fAGCHitsLinearized.insert(fAGCHitsLinearized.end(), fAGCHits.at(i).begin(), fAGCHits.at(i).end());
	}
	}
	
	
	short unsigned int AGCounter::GetNHits() { return fAGCHitsLinearized.size(); }
	
	long unsigned int 	AGCounter::GetHitTimeStampUSE(int iHit) { return fAGCHitsLinearized.at(iHit).HitTimeStampUSE; }
	long unsigned int 	AGCounter::GetHitTimeStampUSW(int iHit) { return fAGCHitsLinearized.at(iHit).HitTimeStampUSW; }
	long unsigned int 	AGCounter::GetHitTimeStampDS1(int iHit) { return fAGCHitsLinearized.at(iHit).HitTimeStampDS1; }
	long unsigned int 	AGCounter::GetHitTimeStampDS2(int iHit) { return fAGCHitsLinearized.at(iHit).HitTimeStampDS2; }
	
	float			AGCounter::GetHitPulseAreaUSE(int iHit) { return fAGCHitsLinearized.at(iHit).HitPulseAreaUSE;}
	float			AGCounter::GetHitPulseAreaUSW(int iHit) { return fAGCHitsLinearized.at(iHit).HitPulseAreaUSW;}
	float			AGCounter::GetHitPulseAreaDS1(int iHit) { return fAGCHitsLinearized.at(iHit).HitPulseAreaDS1;}
	float			AGCounter::GetHitPulseAreaDS2(int iHit) { return fAGCHitsLinearized.at(iHit).HitPulseAreaDS2;}
	
	bool 			AGCounter::GetHitExistUSE(int iHit) { return fAGCHitsLinearized.at(iHit).HitExistUSE; }
	bool 			AGCounter::GetHitExistUSW(int iHit) { return fAGCHitsLinearized.at(iHit).HitExistUSW; }
	bool 			AGCounter::GetHitExistDS1(int iHit) { return fAGCHitsLinearized.at(iHit).HitExistDS1; }
	bool 			AGCounter::GetHitExistDS2(int iHit) { return fAGCHitsLinearized.at(iHit).HitExistDS2; }
*/

   size_t	AGCounter::GetNHits() const { return fAGCHits.size(); }
//  short unsigned int 	AGCounter::GetNHits() { return fAGCHits.size(); }
  
  long unsigned int 	AGCounter::GetHitTimeStampUSE(int iHit) const { return fAGCHits.at(iHit).HitTimeStampUSE; }
  long unsigned int 	AGCounter::GetHitTimeStampUSW(int iHit) const { return fAGCHits.at(iHit).HitTimeStampUSW; }
  long unsigned int 	AGCounter::GetHitTimeStampDS1(int iHit) const { return fAGCHits.at(iHit).HitTimeStampDS1; }
  long unsigned int 	AGCounter::GetHitTimeStampDS2(int iHit) const { return fAGCHits.at(iHit).HitTimeStampDS2; }
  
  float			AGCounter::GetHitPulseAreaUSE(int iHit) const { return fAGCHits.at(iHit).HitPulseAreaUSE;}
  float			AGCounter::GetHitPulseAreaUSW(int iHit) const { return fAGCHits.at(iHit).HitPulseAreaUSW;}
  float			AGCounter::GetHitPulseAreaDS1(int iHit) const { return fAGCHits.at(iHit).HitPulseAreaDS1;}
  float			AGCounter::GetHitPulseAreaDS2(int iHit) const { return fAGCHits.at(iHit).HitPulseAreaDS2;}
  
  bool 			AGCounter::GetHitExistUSE(int iHit) const { return fAGCHits.at(iHit).HitExistUSE; }
  bool 			AGCounter::GetHitExistUSW(int iHit) const { return fAGCHits.at(iHit).HitExistUSW; }
  bool 			AGCounter::GetHitExistDS1(int iHit) const { return fAGCHits.at(iHit).HitExistDS1; }
  bool 			AGCounter::GetHitExistDS2(int iHit) const { return fAGCHits.at(iHit).HitExistDS2; }
}
