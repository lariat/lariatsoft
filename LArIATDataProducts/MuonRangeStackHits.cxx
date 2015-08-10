/////////////////////////////////////////
//
//  MuonRangeStackHits Class
//
//  gkpullia@syr.edu
// COMMENT
///////////////////////////////////////////

#include "LArIATDataProducts/MuonRangeStackHits.h"
#include "cetlib/exception.h"
#include <math.h>
#include <vector>
#include <map>
#include <iostream>

namespace ldp{

//#######################################
MuonRangeStackHits::MuonRangeStackHits()
{
std::map<int,std::vector<int> > paddlemap;
std::vector<int> blank1;
blank1.push_back(0);
for(int i=0; i<16;++i){

paddlemap.emplace(i,blank1);
}
}
//##########################################
MuonRangeStackHits::MuonRangeStackHits(std::map<int, std::vector<int> > paddlemap)
{
fPaddleTimeTickMap=paddlemap;

}

//########################################
std::vector<int> MuonRangeStackHits::TimeTick(int iPaddle) const
{
int lastpaddle=fPaddleTimeTickMap.end()->first;
	if(iPaddle > lastpaddle  ){
	   throw cet::exception("MuonRangeStackHits") << "Requested time tick vector for paddle "<<iPaddle<<".  That doesn't exist for this event.  The last paddle number you can reference is ""<<lastpaddle<<" <<"\n";
	   }
return fPaddleTimeTickMap.find(iPaddle)->second;
}


}// end namespace
