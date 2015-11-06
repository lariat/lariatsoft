class Picky : public WCTrackAlgBase
{
public:
	void Hello() {std::cout<<"Picky"<<std::endl;}
	void InitializeGeometry(){std::cout<<"PICKY OVERRIDE"<<std::endl;}
	static WCTrackAlgBase* Create() {return new Picky();}

	
};
