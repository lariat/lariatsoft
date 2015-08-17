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



namespace ldp {
  
  class AGCounter {
  public:
	AGCounter(); 
    	~AGCounter();

//      AGCounter(std::vector<std::vector<AGCHits> >);
       	AGCounter(std::vector<AGCHits>);
	
	short unsigned int 	GetNHits();
	
	long int		GetTriggerTimeStamp(int);
	
	long unsigned int 	GetHitTimeStampUSE(int);
	long unsigned int 	GetHitTimeStampUSW(int);
	long unsigned int 	GetHitTimeStampDS1(int);
	long unsigned int 	GetHitTimeStampDS2(int);
	
	float			GetHitPulseAreaUSE(int);
	float			GetHitPulseAreaUSW(int);
	float			GetHitPulseAreaDS1(int);
	float			GetHitPulseAreaDS2(int);
	
	bool 			GetHitExistUSE(int);
	bool 			GetHitExistUSW(int);
	bool 			GetHitExistDS1(int);
	bool 			GetHitExistDS2(int);

    
  private:
  	std::vector<AGCHits> fAGCHits;
//	std::vector<AGCHits> 		   fAGCHitsLinearized;
//	void Linearize();	 
  };
}


#endif // LARIATDATAPRODUCTS_AGCOUNTER_H
