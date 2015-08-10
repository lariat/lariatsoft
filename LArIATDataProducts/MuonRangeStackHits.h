//////////////////////////////////////////////////////////
// Definitions of MuonRangeStackHits object
// 
// gkpullia@syr.edu
// 
///////////////////////////////////////////////////////////

#ifndef LARIATDATAPRODUCTS_MUONRANGESTACKHITS_H
#define LARIATDATAPRODUCTS_MUONRANGESTACKHITS_H

#include <vector>
#include <iosfwd>
#include <string>
#include <map>
#include <iostream>

namespace ldp{

class MuonRangeStackHits{
	
public:
MuonRangeStackHits();

private:
std::map<int, std::vector<int> > fPaddleTimeTickMap;
/* int fPaddle; //	paddle number through which particle passed.
std::vector<int> fTimeTick;   // time tick when paddle was hit */

#ifndef __GCCXML__

public:

	MuonRangeStackHits(std::map<int, std::vector<int>> paddlemap);
	/* int              Paddle()                 const; */ //Find paddle number that was hit, 0 to number of paddles. Currently 16 total paddles at the end of Run I (Shutdown 2015)
	std::vector<int>  TimeTick(int iPaddle)  const; //The vector listing the time ticks when a certain iPaddle was hit.

	
#endif 
};
}

#ifndef __GCCXML__





#endif

#endif //LARIATDATAPRODUCTS_MUONRANGESTACKHITS_H
	
