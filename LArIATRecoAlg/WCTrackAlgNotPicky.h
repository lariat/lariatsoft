class NotPicky : public WCTrackAlgBase
{
public:
	void Hello() {std::cout<<"Not Picky"<<std::endl;}
	static WCTrackAlgBase* Create() {return new NotPicky();}
/*         void loadXMLDatabaseTableForBField( int run, int subrun )
{
  fRun = run;
  fSubRun = subrun;
  fB_field_tesla = 0.0035*std::stod(fDatabaseUtility->GetIFBeamValue("mid_f_mc7an",fRun,fSubRun));
  std::cout << "Run: " << fRun << ", Subrun: " << fSubRun << ", B-field: " << fB_field_tesla << std::endl;
} */
};
