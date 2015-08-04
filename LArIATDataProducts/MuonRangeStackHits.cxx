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
std::map<int,std::vector<int>> paddlemap;
std::vector<int> blank1;
blank1.push_back(0);
for(int i=0; i<16;++i){

paddlemap[i] = blank1;
}
}
//##########################################
MuonRangeStackHits::MuonRangeStackHits(std::map<int, std::vector<int>> paddlemap)
{
for(std::map<int,std::vector<int>>::iterator it=paddlemap.begin(); it!=paddlemap.end(); ++it)
{
	std::cout << it->first<< std::endl;
	

}
}

//########################################
std::vector<int> MuonRangeStackHits::TimeTick(std::map<int, std::vector<int>> paddlemap, int iPaddle) const
{
int lastpaddle=paddlemap.end()->first;
	if(iPaddle > lastpaddle  ){
	   throw cet::exception("MuonRangeStackHits") << "Requested time tick vector for paddle "<<iPaddle<<".  That doesn't exist for this event.  The last paddle number you can reference is ""<<lastpaddle<<" <<"\n";
	   }
return paddlemap.find(iPaddle)->second;
}


}// end namespace
