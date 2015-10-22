class Picky : public WCTrackAlgBase
{
public:
	void Hello() {std::cout<<"Picky"<<std::endl;}
	static WCTrackAlgBase* Create() {return new Picky();}
	//void JumpThroughHoopFactory() {this->JumpThroughHoop();}
	
};
