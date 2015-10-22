/* 
//Framework includes
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "CLHEP/Units/SystemOfUnits.h"
#include "Geometry/AuxDetGeo.h"
#include "art/Framework/Services/Optional/TFileService.h"*/

// LArIAT includes
#include "LArIATRecoAlg/WCTrackBuilderAlgBase.h"


#include <iostream>
#include <cstdlib>
#include <string>
#include <map> 
class WCTrackAlgBase
{
   public:
          virtual void Hello() = 0;
	  //virtual void reconfigure( fhicl::ParameterSet const& pset );
	  //virtual void JumpThroughHoopFactory()=0;
	  	 
};

typedef WCTrackAlgBase* (*CreateWCTrackBuilderAlgFn)(void);

class Picky : public WCTrackAlgBase
{
public:
	void Hello() {std::cout<<"Picky"<<std::endl;}
	static WCTrackAlgBase* Create() {return new Picky();}
	//void JumpThroughHoopFactory() {this->JumpThroughHoop();}

};
class NotPicky : public WCTrackAlgBase
{
public:
	void Hello() {std::cout<<"Not Picky"<<std::endl;}
	static WCTrackAlgBase* Create() {return new NotPicky();}

};


class WCTrackBuilderFactory
{
private:
	WCTrackBuilderFactory();
	WCTrackBuilderFactory(const WCTrackBuilderFactory &) {}
	WCTrackBuilderFactory &operator= (const WCTrackBuilderFactory &) {return *this;}
	typedef std::map < std::string, CreateWCTrackBuilderAlgFn > FactoryMap;
	FactoryMap m_FactoryMap;
public:
	~WCTrackBuilderFactory() {m_FactoryMap.clear();}
	static WCTrackBuilderFactory* Get()
	{
		static WCTrackBuilderFactory instance;
		return &instance;
	}
	
	void Register(const std::string &WCTrackBuilderAlgName, CreateWCTrackBuilderAlgFn pfnCreate);
	WCTrackAlgBase* CreateAlg(const std::string &WCTrackBuilderAlgName);
};


WCTrackBuilderFactory::WCTrackBuilderFactory()
{
	Register("Picky", &Picky::Create);
	Register("NotPicky", &NotPicky::Create);
}

void WCTrackBuilderFactory::Register(const std::string &WCTrackBuilderAlgName, CreateWCTrackBuilderAlgFn pfnCreate)
{
	m_FactoryMap[WCTrackBuilderAlgName] = pfnCreate;
}

WCTrackAlgBase *WCTrackBuilderFactory::CreateAlg(const std::string &WCTrackBuilderAlgName)
{
	FactoryMap::iterator it = m_FactoryMap.find(WCTrackBuilderAlgName);
	if (it != m_FactoryMap.end())
	return it->second();
	return NULL;
}


