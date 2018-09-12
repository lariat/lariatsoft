// =================================================
// === cleaned up and sparse version of ana tree ===
// =================================================

// ##########################
// ### Framework includes ###
// ##########################
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "canvas/Persistency/Common/FindOneP.h" 
#include "canvas/Persistency/Common/FindManyP.h"
#include "messagefacility/MessageLogger/MessageLogger.h" 
//#include "cetlib/maybe_ref.h"

// ########################
// ### LArSoft includes ###
// ########################
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "larcore/Geometry/Geometry.h"
//#include "larcore/Geometry/CryostatGeo.h"
//#include "larcore/Geometry/TPCGeo.h"
//#include "larcore/Geometry/PlaneGeo.h"
//#include "larcore/Geometry/WireGeo.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/SpacePoint.h"
//#include "lardata/RecoBaseArt/TrackUtils.h" // lar::util::TrackPitchInView()
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"

//#include "RawData/ExternalTrigger.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larevt/Filters/ChannelFilter.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "LArIATDataProducts/WCTrack.h"
#include "LArIATDataProducts/TOF.h"
#include "LArIATDataProducts/AGCounter.h"
#include "RawDataUtilities/TriggerDigitUtility.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCStep.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"

// #####################
// ### ROOT includes ###
// #####################
#include "TComplex.h"
#include "TFile.h"
#include "TH2D.h"
#include "TF1.h"
#include "TTree.h"
#include "TTimeStamp.h"

const int kMaxTOF        = 100;      // :/
const int kMaxWCTracks   = 1000;
const int kMaxTrack      = 1000;     //maximum number of tracks
const int kMaxHits       = 20000;    //maximum number of hits
const int kMaxTrackHits  = 1000;     //maximum number of space points
const int kMaxHitIDs     = 100;      //maximum number of space points ids
const int kMaxTrajHits   = 1000;     //maximum number of trajectory points
const int kMaxCluster    = 1000;     //maximum number of clusters
const int kMaxPrimaryPart = 10;	     //maximum number of true primary particles
const int kMaxDaughterPart = 100;     //maximum number of true daughter particles
const int kMaxPrimaries  = 20000;    //maximum number of true particles tracked
const int kMaxTruePrimaryPts = 5000; //maximum number of points in the true primary trajectory 
const int kMaxTrueDaughterPts = 5000; //maximum number of points in the true daughter trajectory 
const int kMaxIDE = 5000; //maximum number of points in the true primary trajectory 

namespace lariat 
{
  class AnaTreeT1034UC;
}

class lariat::AnaTreeT1034UC : public art::EDAnalyzer 
{
public:
  explicit AnaTreeT1034UC(fhicl::ParameterSet const & p);
  virtual ~AnaTreeT1034UC();

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob();
  void reconfigure(fhicl::ParameterSet const & p);
  //void beginJob() override;
  //void beginRun(art::Run const & r) override;
  //void beginSubRun(art::SubRun const & sr) override;
  //void endJob() override;
  //void endRun(art::Run const & r) override;
  //void endSubRun(art::SubRun const & sr) override;
  //void reconfigure(fhicl::ParameterSet const & p) override;
  //void initializeVariables();

private:

  // === Function used to reset all the variables  ===
  void ResetVars();
  
  // === Storing information into TTree ====
  TTree* fTree;
  // === Histos! ===
  TH2D* E_vs_e;
   
  //=== Storing Run Information ===
  int run;			//<---Run Number
  int subrun;			//<---SubRun Number
  int event;			//<---Event Number
  double evttime;		//<---Event Time Stamp
  double efield[3];		//<---Electric Field 
  int t0;
   
  // === Storing Track Information ===
  int ntracks_reco;		//<---Number of reconstructed tracks
  std::vector<int>    track_WC2TPC_match;
  std::vector<double> track_start_x;
  std::vector<double> track_start_y;
  std::vector<double> track_start_z;
  std::vector<double> track_end_x;
  std::vector<double> track_end_y;
  std::vector<double> track_end_z;
  
  std::vector<double> track_length;
  std::vector<int> track_primary;	//<---Is this particle primary (primary = 1, non-primary = 0)

  // === Storing the tracks SpacePoints (Individual 3D points)
  std::vector<int> ntrack_hits;
  std::vector< std::vector<double> > track_xpos;
  std::vector< std::vector<double> > track_ypos;
  std::vector< std::vector<double> > track_zpos;
  std::vector< std::vector<int> > nhit_ids;
  double track_spt_idarr[kMaxTrack][kMaxTrackHits][kMaxHitIDs];
  double track_spt_earr[kMaxTrack][kMaxTrackHits][kMaxHitIDs];

  std::vector< std::vector<std::map<int, double>> > track_spt_g4map;
  std::vector< std::vector<std::vector<int>> > track_spt_id;
  std::vector< std::vector<std::vector<double>> > track_spt_e;



  // === Storing the tracks Calorimetry Information
  std::vector<int> ind_track_hits;
  std::vector<double> ind_track_ke;
  std::vector< std::vector<double> > ind_track_wire;
  std::vector< std::vector<double> > ind_track_dedx;
  std::vector< std::vector<double> > ind_track_dqdx;
  std::vector< std::vector<double> > ind_track_rr;
  std::vector< std::vector<double> > ind_track_pitch_hit;
  std::vector<int> col_track_hits;
  std::vector<double> col_track_ke;
  std::vector< std::vector<double> > col_track_x;
  std::vector< std::vector<double> > col_track_y;
  std::vector< std::vector<double> > col_track_z;
  std::vector< std::vector<double> > col_track_wire;
  std::vector< std::vector<double> > col_track_dedx;
  std::vector< std::vector<double> > col_track_dqdx;
  std::vector< std::vector<double> > col_track_rr;
  std::vector< std::vector<double> > col_track_pitch_hit;


  // === hit info ===
  int nhits;
  std::vector<double> hit_time;
  std::vector<double> hit_wire;
  std::vector<double> hit_view;
  std::vector<double> hit_amp;

  std::vector<int>    InteractionPoint;         //<---Geant 4 Primary Trj Point Corresponding to the Interaction
  std::vector<int>    InteractionPointType;     //<---Geant 4 Primary Interaction Type


  // === Storing Geant4 MC Truth Information ===
  int no_primaries;				//<---Number of primary Geant4 particles in the event
  int geant_list_size;			//<---Number of Geant4 particles tracked
  double primary_p;				//<---Primary particle momentum
  std::vector<int> PDG;
  std::vector<double> StartPointx; 
  std::vector<double> StartPointy;
  std::vector<double> StartPointz;
  std::vector<double> StartEnergy;
  std::vector<double> StartKE;
  std::vector<double> LastKE;
  std::vector<double> StartPx;
  std::vector<double> StartPy;
  std::vector<double> StartPz;

  std::vector<double> EndPointx; 
  std::vector<double> EndPointy;
  std::vector<double> EndPointz;
  std::vector<double> EndEnergy;
  std::vector<double> EndPx;
  std::vector<double> EndPy;
  std::vector<double> EndPz;

  std::vector<int> Process;
  // ### Recording the process as a integer ###
	  // 0 = primary
	  // 3 = hadElastic
	  // 8 = CoulombScat
	  //10 = ProtonInelastic



  std::vector<int> NumberDaughters;
  std::vector<int> TrackId;
  std::vector<int> Mother;
  std::vector<int> process_primary;	//<---Is this particle primary (primary = 1, non-primary = 0)
  std::vector<std::string> G4Process;         //<---The process which created this particle
  std::vector<std::string> G4FinalProcess;    //<---The last process which this particle went under

  // === Storing additional Geant4 MC Truth Information for the primary track only ===	   
  std::vector<int> NTrTrajPts;
  std::vector< std::vector<double> > MidPosX;
  std::vector< std::vector<double> > MidPosY;
  std::vector< std::vector<double> > MidPosZ;
  std::vector< std::vector<double> > MidPx;
  std::vector< std::vector<double> > MidPy;
  std::vector< std::vector<double> > MidPz;
  std::vector<double> SlabX;
  std::vector<double> SlabY;
  std::vector<double> SlabZ;
  std::vector<double> SlabE;
  std::vector<int>    SlabN;
  std::vector<double> SlapX;
  std::vector<double> SlapY;
  std::vector<double> SlapZ;
  int LastSlabInt;
  //// === changing how storing geant4 mc truth info for primary track only ===
  //std::vector<double> PriMidPosX;
  //std::vector<double> PriMidPosY;
  //std::vector<double> PriMidPosZ;
  //std::vector<double> PriMidPx;
  //std::vector<double> PriMidPy;
  //std::vector<double> PriMidPz;


  std::vector<double> NDTrTrajPts;
  int NProtonDaughters = 0;
  int NNeutronDaughters = 0;
  std::vector<int> DTrackId;
  std::vector<int> DPdgCode;
  std::vector<double> DStartKE;
  std::vector<double> DStartEnergy;
  std::vector<double> DStartP;
  std::vector< std::vector<double> > DMidPosX;
  std::vector< std::vector<double> > DMidPosY;
  std::vector< std::vector<double> > DMidPosZ;
  // === Storing additional Geant4 MC Truth Information for the daughter tracks ===  
   

 // === beamline info ===
 int num_tof_objects;
 double tofObject[kMaxTOF];

 int num_wctracks;
 double wctrk_momentum[kMaxWCTracks];
 double wctrk_XFace[kMaxWCTracks];
 double wctrk_YFace[kMaxWCTracks];
 double wctrk_theta[kMaxWCTracks];
 double wctrk_phi[kMaxWCTracks];
 int wctrk_missed[kMaxWCTracks];
 int wctrk_picky[kMaxWCTracks];

  double electron_lifetime;
  
   
  std::string fTreeName;
  std::string fClusterModuleLabel;
  std::string fHitsModuleLabel;
  std::string fTrackModuleLabel;
  std::string fCalorimetryModuleLabel;
  std::string fParticleIDModuleLabel;
  std::string fG4ModuleLabel;
  std::string fTOFModuleLabel;
  std::string fWCTrackLabel;
  std::string fWC2TPCModuleLabel;

  calo::CalorimetryAlg fCalorimetryAlg;

};


lariat::AnaTreeT1034UC::AnaTreeT1034UC(fhicl::ParameterSet const & pset) 
  : EDAnalyzer(pset)
  , fCalorimetryAlg(pset.get<fhicl::ParameterSet>("CalorimetryAlg"))
{
  this->reconfigure(pset);
}

lariat::AnaTreeT1034UC::~AnaTreeT1034UC()
{
  // Clean up dynamic memory and other resources here.
}

void lariat::AnaTreeT1034UC::reconfigure(fhicl::ParameterSet const & pset)
{
  fTreeName                 = pset.get< std::string >("TreeName", "anatreeuc");

  fHitsModuleLabel          = pset.get< std::string >("HitsModuleLabel");
  fTrackModuleLabel         = pset.get< std::string >("TrackModuleLabel");
  fCalorimetryModuleLabel   = pset.get< std::string >("CalorimetryModuleLabel");
  fParticleIDModuleLabel    = pset.get< std::string >("ParticleIDModuleLabel");
  fClusterModuleLabel       = pset.get< std::string >("ClusterModuleLabel");
  fG4ModuleLabel            = pset.get< std::string >("G4ModuleLabel");

  fTOFModuleLabel           = pset.get< std::string >("TOFModuleLabel");
  fWCTrackLabel             = pset.get< std::string >("WCTrackLabel");
  fWC2TPCModuleLabel        = pset.get< std::string >("WC2TPCModuleLabel", "WC2TPCtrk");

  return;
}

void lariat::AnaTreeT1034UC::analyze(art::Event const & evt)
{
  // #############################################
  // ### Reset variables before we get started ###
  // #############################################
  ResetVars();
  //std::cout<<"Check1"<<std::endl;

  // #######################################
  // ### Get potentially useful services ###
  // #######################################
  // === Geometry Service ===
  art::ServiceHandle<geo::Geometry> geom;
  // === Liquid Argon Properties Services ===
  //auto const* larprop = lar::providerFrom<detinfo::LArPropertiesService>();
  // === Detector properties service ===
  auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
  // === BackTracker service ===
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  const sim::ParticleList& plist = pi_serv->ParticleList();
   
   
  // === Run Number ===
  run = evt.run();
  // === Sub-Run Number ===
  subrun = evt.subRun();
  // === Event Number ===
  event = evt.id().event();
   
  std::cout<<std::endl;
  std::cout<<"========================================="<<std::endl;
  std::cout << "Run = "         << run 
            << ", SubRun = "    << subrun 
            << ", Evt = "       << event    << std::endl;
  std::cout<<"========================================="<<std::endl;
  std::cout<<std::endl;
   
  // === Event Time ===
  art::Timestamp ts = evt.time();
  TTimeStamp tts(ts.timeHigh(), ts.timeLow());
  evttime = tts.AsDouble();

  // === Electron lifetime ===
  electron_lifetime = detprop->ElectronLifetime();

  // === Electric Field ===
  // Note: LArProperties::Efield() has moved to DetectorProperties/DetectorPropertiesService
  efield[0] = detprop->Efield(0);
  efield[1] = detprop->Efield(1);
  efield[2] = detprop->Efield(2);
   
  // === Trigger Offset ====
  t0 = detprop->TriggerOffset();
   
  // #####################################
  // ### Getting the Track Information ###
  // #####################################
  art::Handle< std::vector<recob::Track> > trackListHandle; //<---Define trackListHandle as a vector of recob::Track objects
  std::vector<art::Ptr<recob::Track> > tracklist; //<---Define tracklist as a pointer to recob::tracks
   
  // === Filling the tracklist from the tracklistHandle ===
  if (evt.getByLabel(fTrackModuleLabel,trackListHandle))
    {art::fill_ptr_vector(tracklist, trackListHandle);}

  // ###################################
  // ### Getting the Hit Information ###
  // ###################################
  art::Handle< std::vector<recob::Hit> > hitListHandle; //<---Define hitListHandle as a vector of recob::hits objects
  std::vector<art::Ptr<recob::Hit> > hitlist; //<---Define hitklist as a pointer to recob::hit
   
  // === Filling the hitlist from the hitlistHandle ===
  if (evt.getByLabel(fHitsModuleLabel,hitListHandle))
    {art::fill_ptr_vector(hitlist, hitListHandle);}


  // #############################
  // ### beam line information ###
  // #############################
  art::Handle< std::vector<ldp::TOF> > TOFColHandle;
  std::vector< art::Ptr<ldp::TOF> > tof;
  if(evt.getByLabel(fTOFModuleLabel, TOFColHandle))
    {art::fill_ptr_vector(tof, TOFColHandle);}

  art::Handle< std::vector<ldp::WCTrack> > wctrackHandle;
  std::vector<art::Ptr<ldp::WCTrack> > wctrack;
  if(evt.getByLabel(fWCTrackLabel, wctrackHandle))
    {art::fill_ptr_vector(wctrack, wctrackHandle);}


   
  // ##########################################################
  // ### Grabbing associations for use later in the AnaTool ###
  // ##########################################################
  //std::cout<<"Check2"<<std::endl;
  // === Associations between hits and raw digits ===
  art::FindOne<raw::RawDigit>       ford(hitListHandle,   evt, fHitsModuleLabel);
  // === Association between SpacePoints and Tracks ===
  art::FindManyP<recob::SpacePoint> fmsp(trackListHandle, evt, fTrackModuleLabel);
  // === Association between Tracks and 2d Hits ===
  art::FindManyP<recob::Track>       fmtk(hitListHandle,   evt, fTrackModuleLabel);
  // === Association between Calorimetry objects and Tracks ===
  art::FindManyP<anab::Calorimetry> fmcal(trackListHandle, evt, fCalorimetryModuleLabel);
  // === Association between Particle ID objects (PID) and Tracks ===
  art::FindManyP<anab::ParticleID>  fmpid(trackListHandle, evt, fParticleIDModuleLabel);
  // ==== Association between Tracks and Hits
  art::FindManyP<recob::Hit> fmth(trackListHandle, evt, fTrackModuleLabel);
  art::FindManyP<recob::Hit, recob::TrackHitMeta> fmthm(trackListHandle, evt, fTrackModuleLabel);
  
   

    // === not my code... dunno what this is ===
//  // ### Something to do with SimChannels...need to come back to ###
//  std::vector<const sim::SimChannel*> fSimChannels;
//  try
//    {evt.getView("largeant", fSimChannels);}
//  catch (art::Exception const&e){ }
    


  // ##################################################################
  // ### Getting SimChannel Info ??? Want to follow reconstruction! ###
  // ##################################################################
//  auto simChannelHandle
//    = event.getValidHandle<std::vector<sim::SimChannel>>
//    (fSimulationProducerLabel);



   
  // ###################################################################
  // ### Setting a boolian to only output MC info if this is MC-info ###
  // ###################################################################
  bool isdata = false;
  if (evt.isRealData())
    {isdata = true;}
	
  else isdata = false;
   

   
  // ----------------------------------------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------------------------------------
  //							FILLING THE MCTruth Geant4 INFORMATION
  // ----------------------------------------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------------------------------------

  if(!isdata)
    {
    // ### Making a vector of MCParticles ###
    std::vector<const simb::MCParticle* > geant_part;
    for(size_t p = 0; p < plist.size(); ++p) 
	  {geant_part.push_back(plist.Particle(p));}
    
    // ### Setting a string for primary ###
    std::string pri("primary");
    std::string hadElastic("hadElastic");
    std::string CoulombScat("CoulombScat");
    std::string ProtonInelastic("protonInelastic");
    
    int primary=0;
    int geant_particle=0;
    int numTrue_inTPC = 0;
    double slab_size_z = 0.5;
    
    // ### Determine the number of primary particles from geant ###
    for( unsigned int i = 0; i < geant_part.size(); ++i )
	  {
	  geant_particle++;
	  if(geant_part[i]->Process()==pri)
	    {primary++;}
	  }//<---End i loop
	
    // ### Saving the number of primary particles ###
    no_primaries=primary;
    // ### Saving the number of Geant4 particles ###
    geant_list_size=geant_particle;
     
    // ### Looping over all the Geant4 particles ###
    int iPrim = 0;
    //int iDaught = 0;
    std::cout<<"geant par size: "<<geant_part.size()<<std::endl;
    for(unsigned int i = 0; i < geant_part.size(); ++i)
      {
	  if(geant_part[i]->Process()==pri){process_primary.push_back(1);}
	  else                             {process_primary.push_back(10);}
           
	  // ### Recording the process as a integer ###
	  // 0 = primary
	  // 3 = hadElastic
	  // 8 = CoulombScat
	  //10 = ProtonInelastic
	   
	  if(geant_part[i]->Process() == pri)            {Process.push_back(0);}
	  if(geant_part[i]->Process() == hadElastic)     {Process.push_back(3);}
	  if(geant_part[i]->Process() == CoulombScat)	 {Process.push_back(8);}
	  if(geant_part[i]->Process() == ProtonInelastic){Process.push_back(10);}

	  // ### Saving other info ###
      PDG.push_back(geant_part[i]->PdgCode());
      Mother.push_back(geant_part[i]->Mother());
      TrackId.push_back(geant_part[i]->TrackId());
      StartEnergy.push_back(geant_part[i]->E());
	  EndEnergy.push_back(geant_part[i]->EndE());
      double startp = sqrt( pow(geant_part[i]->Px(),2) + pow(geant_part[i]->Py(),2) + pow(geant_part[i]->Pz(),2));
      double mass = geant_part[i]->Mass(); 
      double startke = sqrt( pow(startp,2) + pow(mass,2) ) - mass;
      StartKE.push_back(startke);
      double lastke = 0;
	  simb::MCTrajectory truetraj = geant_part[i]->Trajectory();
      int counter_ii = 0;
	  for(auto itTraj = truetraj.begin(); itTraj != truetraj.end(); ++itTraj) {
        if(truetraj.Z(counter_ii) > 0   && truetraj.Z(counter_ii) < 90   &&
          truetraj.X(counter_ii) > 0   && truetraj.X(counter_ii) < 47.5 &&
          truetraj.Y(counter_ii) > -20 && truetraj.Y(counter_ii) < 20) {
          double thisp;
          if(counter_ii) {
            thisp = sqrt( pow(truetraj.Px(counter_ii-1),2) + 
                          pow(truetraj.Py(counter_ii-1),2) + 
                          pow(truetraj.Pz(counter_ii-1),2));
          }
          else {
            thisp = sqrt( pow(truetraj.Px(counter_ii),2) + 
                          pow(truetraj.Py(counter_ii),2) + 
                          pow(truetraj.Pz(counter_ii),2));
          }
          lastke = sqrt( pow(thisp,2) + pow(mass,2) ) - mass;
        }//<--End if in tpc
        counter_ii++;
	  }//<--End loop on true trajectory points
      LastKE.push_back(lastke);


	   
	  // ### Saving the start and end Px, Py, Pz info ###
	  StartPx.push_back(geant_part[i]->Px());
	  StartPy.push_back(geant_part[i]->Py());
	  StartPz.push_back(geant_part[i]->Pz());
	  EndPx.push_back(geant_part[i]->EndPx());
	  EndPy.push_back(geant_part[i]->EndPy());
	  EndPz.push_back(geant_part[i]->EndPz());
	   
      //std::cout << "    p(x,y,z): (" << geant_part[i]->Px() << ", " << geant_part[i]->Py() << ", " << geant_part[i]->Px() << ")     " 
      //          << "    start E: " << geant_part[i]->E() << std::endl;
	  // ### Saving the Start and End Point for this particle ###
	  StartPointx.push_back(geant_part[i]->Vx());
	  StartPointy.push_back(geant_part[i]->Vx());
	  StartPointz.push_back(geant_part[i]->Vx());
	  EndPointx.push_back(geant_part[i]->EndPosition()[0]);
	  EndPointy.push_back(geant_part[i]->EndPosition()[1]);
	  EndPointz.push_back(geant_part[i]->EndPosition()[2]);

	  // ### Saving the processes for this particle ###
	  G4Process.push_back( geant_part[i]->Process() );
	  G4FinalProcess.push_back( geant_part[i]->EndProcess() );
 	   
	  // ### Saving the number of Daughters for this particle ###
	  NumberDaughters.push_back(geant_part[i]->NumberDaughters());

	  // ### Save intermediary information for the primary track
	  if(geant_part[i]->Process() == pri)
        {
        //std::cout << "    p(x,y,z): (" << geant_part[i]->Px() << ", " << geant_part[i]->Py() << ", " << geant_part[i]->Px() << ")" << std::endl; 
        //std::cout << "    sqrt(p2):  " << sqrt( pow(geant_part[i]->Px(), 2) + pow(geant_part[i]->Py(), 2) + pow(geant_part[i]->Pz(), 2) ) << std::endl;
        //std::cout << "    p       :  " << geant_part[i]->P() << std::endl;
        primary_p = geant_part[i]->P();
	    NTrTrajPts.push_back(geant_part[i]->NumberTrajectoryPoints());
	    //simb::MCTrajectory truetraj = geant_part[i]->Trajectory();
        std::vector<double> midx;  std::vector<double> midy;  std::vector<double> midz;
        std::vector<double> midpx; std::vector<double> midpy; std::vector<double> midpz;
        std::vector<double> slabx; std::vector<double> slaby; std::vector<double> slabz;
        std::vector<double> slabe; std::vector<int> slabn;
        std::vector<double> slapx; std::vector<double> slapy; std::vector<double> slapz; 
	    int iPrimPt = 0;	
        double true_total_distance = 0;
        double previous_mid_xpt = geant_part[i]->Vx();    
        double previous_mid_ypt = geant_part[i]->Vy();    
        double previous_mid_zpt = geant_part[i]->Vz();    
	    for(auto itTraj = truetraj.begin(); itTraj != truetraj.end(); ++itTraj)
	      {
          // ### true pt vars ###
          double mid_xpt = truetraj.X(iPrimPt);
          double mid_ypt = truetraj.Y(iPrimPt);
          double mid_zpt = truetraj.Z(iPrimPt);
          double true_distance = sqrt( pow(mid_xpt - previous_mid_xpt, 2) 
                                     + pow(mid_ypt - previous_mid_ypt, 2) 
                                     + pow(mid_zpt - previous_mid_zpt, 2));
          if(mid_zpt > 0   && mid_zpt < 90   &&
             mid_xpt > 0   && mid_xpt < 47.5 &&
             mid_ypt > -20 && mid_ypt < 20){numTrue_inTPC++;}
          true_total_distance += true_distance; 
          std::cout<<"\t\ttrue total distance (highres): "<<true_total_distance<<std::endl;
          //if(fmod(true_total_distance,1) < 0.03)
          if(fmod(true_total_distance,slab_size_z) < 0.03)
            {
            if(mid_zpt > 0   && mid_zpt < 90   &&
               mid_xpt > 0   && mid_xpt < 47.5 &&
               mid_ypt > -20 && mid_ypt < 20)
              {
              std::cout<<"totalD: "<<true_total_distance<<std::endl;
              std::cout<<"dist_mod_slab: "<<fmod(true_total_distance,slab_size_z)<<std::endl;
              std::cout<<"\tx,y,z: "<<mid_xpt<<", "<<mid_ypt<<", "<<mid_zpt<<std::endl;
              slabx.push_back(mid_xpt);
              slaby.push_back(mid_ypt);
              slabz.push_back(mid_zpt);
              slabn.push_back((int)(true_total_distance/0.5));
              slabe.push_back(truetraj.E(iPrimPt));
              slapx.push_back(truetraj.Px(iPrimPt));
              slapy.push_back(truetraj.Py(iPrimPt));
              slapz.push_back(truetraj.Pz(iPrimPt));
              }
            }
                                 
                                 
          previous_mid_xpt = mid_xpt;
          previous_mid_ypt = mid_ypt;
          previous_mid_zpt = mid_zpt;
          // ### pushing back vars ###
          midx.push_back(truetraj.X(iPrimPt));
          midy.push_back(truetraj.Y(iPrimPt));
          midz.push_back(truetraj.Z(iPrimPt));
          midpx.push_back(truetraj.Px(iPrimPt));
          midpy.push_back(truetraj.Py(iPrimPt));
          midpz.push_back(truetraj.Pz(iPrimPt));
	      iPrimPt++;
          // ## debug ##
          //std::cout << "        point: "            << iPrimPt << std::endl;
          //std::cout << "        mid position px: "  << truetraj.Px(iPrimPt) << std::endl;
          //std::cout << "        mid position py: "  << truetraj.Py(iPrimPt) << std::endl;
          //std::cout << "        mid position pz: "  << truetraj.Pz(iPrimPt) << std::endl;
          //std::cout << "        mid position P : "  << sqrt(   pow(truetraj.Px(iPrimPt), 2) 
          //                                                  + pow(truetraj.Py(iPrimPt), 2) 
          //                                                  + pow(truetraj.Pz(iPrimPt), 2) ) << std::endl;
          //std::cout << "        mid position e : "  << truetraj.E(iPrimPt)  << std::endl;
          //std::cout << "        (x,y,z): ("         << truetraj.X(iPrimPt) << ", " 
          //                                          << truetraj.Y(iPrimPt) << ", " 
          //                                          << truetraj.Z(iPrimPt) << ") " << std::endl;
          //std::cout << "        distance: "         << true_distance << std::endl;
          //std::cout << "        total distance: "   << true_total_distance << std::endl;
          //std::cout << std::endl;
	      }//<--End loop on true trajectory points
        std::cout << "      number of trajectory points: " << geant_part[i]->NumberTrajectoryPoints() << std::endl;
        std::cout << "      prim pt                  : "   << iPrimPt << std::endl;
        std::cout << "      total distance: "              << true_total_distance << std::endl;
        std::cout << "\t\tnumslabs: " << slabn.size() << std::endl;
        for(unsigned int i = 2; i < slabn.size(); i++)
          {
          if(i == 2)
            {
            std::cout<<"\t\t\tslabz: "<<slabz[i]<<std::endl;
            SlabN.push_back(slabn[i]);
            SlabE.push_back(slabe[i]);
            SlabX.push_back(slabx[i]);
            SlabY.push_back(slaby[i]);
            SlabZ.push_back(slabz[i]);
            SlapX.push_back(slapx[i]);
            SlapY.push_back(slapy[i]);
            SlapZ.push_back(slapz[i]);
            }
          else
            {
            if(slabn[i] == slabn[i-1])
              {continue;}
            else
              {
              std::cout<<"\t\t\tslabz: "<<slabz[i]<<std::endl;
              SlabN.push_back(slabn[i]);
              SlabE.push_back(slabe[i]);
              SlabX.push_back(slabx[i]);
              SlabY.push_back(slaby[i]);
              SlabZ.push_back(slabz[i]);
              SlapX.push_back(slapx[i]);
              SlapY.push_back(slapy[i]);
              SlapZ.push_back(slapz[i]);
              }
            }
          }
        std::cout<<"SlabN.size() :"<<SlabN.size()<<std::endl;
        MidPosX.push_back(midx); MidPosY.push_back(midy); MidPosZ.push_back(midz);
        MidPx.push_back(midpx);  MidPy.push_back(midpy);  MidPz.push_back(midpz);
        midx.clear(); midy.clear(); midz.clear();
        midpx.clear(); midpy.clear(); midpz.clear();
        slabx.clear(); slaby.clear(); slabz.clear();
        slabn.clear(); slabe.clear();
	    
	    // ### Recording the process as a integer ###
	    // 0 = NoInteractionNodaughters, thought going
	    // 3 = hadElastic
	    // 8 = CoulombScat
	    //10 = ProtonInelastic
	    //13 = protonInelastic
	    auto thisTracjectoryProcessMap =  truetraj.TrajectoryProcesses();
	    // Ok, if the size of the map is 0, all the action might happen at the end of the track
	    // So we check the daugthers:
	    //    - Case 1. There are daugthers:
	    //               * The interesting point is the last one
	    //               * The interaction type is the one that created the first daugther (this might be improved)
	    //    - Case 2. There are NO daugthers:
	    //              * We assign the interaction type to be 0: nothing happens, thought going particle
	    //              * The interesting point is the last one (might not be in the TPC)
	    if(!thisTracjectoryProcessMap.size())
	      {
	      int interestingPoint = (int) (NTrTrajPts[i] - 1);
	      InteractionPoint.push_back(interestingPoint);
	      if(NumberDaughters[i])
	        {		 
	        auto thePrimaryDaughterID = geant_part[i]-> Daughter(0); 
	        for(unsigned int iD = 0; iD < geant_part.size(); ++iD )
	          {
	          if(geant_part[iD]->TrackId() == thePrimaryDaughterID) 
	     	    {		      
	     	    if(geant_part[iD]->Process() == hadElastic)	    {InteractionPointType.push_back(3);}
	     	    if(geant_part[iD]->Process() == CoulombScat)	{InteractionPointType.push_back(8);}
	     	    if(geant_part[iD]->Process() == ProtonInelastic){InteractionPointType.push_back(10);}
	     	    }//<--- End if is daughter(0)
	          }//<--- End particle loop
	        }//<--- End if there are daughters
          else
	        {InteractionPointType.push_back(0);}
	      }
          else
	        {
	        // The map is not zero: somthing interesting might happen in the middle of the track!!
	        for(auto const& couple: thisTracjectoryProcessMap) 
	          {
	          int interestingPoint = (int) couple.first;
	          InteractionPoint.push_back(interestingPoint);         	   
	          if ((truetraj.KeyToProcess(couple.second)).find("hadElastic")      != std::string::npos) InteractionPointType.push_back(3);           
	          if ((truetraj.KeyToProcess(couple.second)).find("pi-Inelastic")    != std::string::npos) InteractionPointType.push_back(1);           
	          if ((truetraj.KeyToProcess(couple.second)).find("pi+Inelastic")    != std::string::npos) InteractionPointType.push_back(14);           
	          if ((truetraj.KeyToProcess(couple.second)).find("kaon-Inelastic")  != std::string::npos) InteractionPointType.push_back(12);           
	          if ((truetraj.KeyToProcess(couple.second)).find("kaon+Inelastic")  != std::string::npos) InteractionPointType.push_back(11);           
	          if ((truetraj.KeyToProcess(couple.second)).find("protonInelastic") != std::string::npos) InteractionPointType.push_back(13);           
	          if ((truetraj.KeyToProcess(couple.second)).find("neutronInelastic")!= std::string::npos) InteractionPointType.push_back(2);           
                
              //std::cout << std::endl;
              //std::cout << "                interaction point pt  : " << (int)couple.first << std::endl;
              //std::cout << "                interaction point P   : " << sqrt(   pow(1000*truetraj.Px((int)couple.first), 2) 
              //                                                                 + pow(1000*truetraj.Py((int)couple.first), 2) 
              //                                                                 + pow(1000*truetraj.Pz((int)couple.first), 2) ) << std::endl;
              //std::cout << "                interaction point P-1 : " << sqrt(   pow(1000*truetraj.Px((int)couple.first-1), 2) 
              //                                                                 + pow(1000*truetraj.Py((int)couple.first-1), 2) 
              //                                                                 + pow(1000*truetraj.Pz((int)couple.first-1), 2) ) << std::endl;
              //
              //std::cout << std::endl;
	          }
	        }	
   
	    iPrim++;
	    }//<--End if primary

	  if(geant_part[i]->Process() != pri)
        {
        if(geant_part[i]->Mother() == 1) 
          {
          //std::cout << "    DPdgCode: " << geant_part[i]->PdgCode() << std::endl;
          //if(geant_part[i]->PdgCode() < 10000){DPdgCode.push_back(geant_part[i]->PdgCode());}
          //DPdgCode.push_back(geant_part[i]->PdgCode());
          //DStartEnergy.push_back(geant_part[i]->E());
	      //NDTrTrajPts.push_back(geant_part[i]->NumberTrajectoryPoints());
          if(geant_part[i]->PdgCode() == 2112)
            {NNeutronDaughters++;}
          if(geant_part[i]->PdgCode() == 2212)
            {NProtonDaughters++;}
          if(geant_part[i]->PdgCode() < 10000)
            {
	        NDTrTrajPts.push_back(geant_part[i]->NumberTrajectoryPoints());
            DTrackId.push_back(geant_part[i]->TrackId());
            DPdgCode.push_back(geant_part[i]->PdgCode());
            DStartEnergy.push_back(geant_part[i]->E());
            DStartP.push_back(geant_part[i]->P());

            std::cout<<"daughter info: "<<std::endl;
            std::cout<<"\tpdg: "<<geant_part[i]->PdgCode()<<std::endl;
            std::cout<<"\tenergy: "<<geant_part[i]->E()<<std::endl;
            std::cout<<"\tmomentum: "<<geant_part[i]->P()<<std::endl;

	        simb::MCTrajectory truetraj = geant_part[i]->Trajectory();
	        int jDaughtPt = 0;	
            std::vector<double> midx; std::vector<double> midy; std::vector<double> midz;
	        for(auto itTraj = truetraj.begin(); itTraj != truetraj.end(); ++itTraj)
	          {
              midx.push_back(truetraj.X(jDaughtPt));
              midy.push_back(truetraj.Y(jDaughtPt));
              midz.push_back(truetraj.Z(jDaughtPt));
	          jDaughtPt++;
	          }//<--End loop on true trajectory points
            DMidPosX.push_back(midx); DMidPosY.push_back(midy); DMidPosZ.push_back(midz);
            midx.clear(); midy.clear(); midz.clear();
            }//<-- End if proton
          }//<-- End if primary daughter
        }//<-- End if not primary
      }//<--End loop on geant particles

      // ## workspace to add *intersting* slab branch? ##
      std::cout << "\ninteractions debug" << std::endl;
      //if(InteractionPointType[InteractionPoint.size()-1] == 13)
      //  {
      //  if(SlabN.size())
      //    {
      //    double int_x = MidPosX[0][InteractionPoint[InteractionPoint.size()-1]];
      //    double int_y = MidPosY[0][InteractionPoint[InteractionPoint.size()-1]];
      //    double int_z = MidPosZ[0][InteractionPoint[InteractionPoint.size()-1]];
      //    double dist_last_slab = sqrt( pow(int_x - SlabX[SlabN.size()-1], 2) +
      //                                  pow(int_y - SlabY[SlabN.size()-1], 2) +
      //                                  pow(int_z - SlabZ[SlabN.size()-1], 2) );
      //    std::cout << "\t\tinelastic int: " << int_x << ", " << int_y << ", " << int_z << std::endl;
      //    std::cout << "\t\t\tdist to last slab: " << dist_last_slab << std::endl;
      //    if(dist_last_slab < 1)
      //      {
      //      std::cout << "\t\t\t\tlast slab goes in int..." << std::endl;
      //      LastSlabInt = 1;
      //      }
      //    else
      //      {
      //      std::cout << "\t\t\t\tlast slab too far from int..." << std::endl;
      //      LastSlabInt = 0;
      //      }
      //    }
      //  }
      //else
      //  {
      //  std::cout << "\tno inelastic int..." << std::endl;
      //  LastSlabInt = 0;
      //  } 
      if(SlabN.size())
        {
        if(InteractionPointType[InteractionPoint.size()-1] == 13)
          {
          double int_x = MidPosX[0][InteractionPoint[InteractionPoint.size()-1]];
          double int_y = MidPosY[0][InteractionPoint[InteractionPoint.size()-1]];
          double int_z = MidPosZ[0][InteractionPoint[InteractionPoint.size()-1]];
          std::cout << "inelastic int: " << int_x << ", " << int_y << ", " << int_z << std::endl;
          std::cout << "slab loop" << std::endl;
          for(unsigned int slab = 0; slab < SlabN.size(); slab++)
            {
            double dist_slab = sqrt( pow(int_x - SlabX[slab], 2) +
                                     pow(int_y - SlabY[slab], 2) +
                                     pow(int_z - SlabZ[slab], 2) );
            std::cout<<"\tslab x,y,z: "<<SlabX[slab]<<", "<<SlabY[slab]<<", "<<SlabZ[slab]<<std::endl;
            std::cout << "\t\tdist to int: " << dist_slab << std::endl; 
            }
          }
        else
          {
          double lastx = MidPosX[0][numTrue_inTPC];
          double lasty = MidPosY[0][numTrue_inTPC];
          double lastz = MidPosZ[0][numTrue_inTPC];
          std::cout << "didn't interact..." << std::endl;
          std::cout << "last pt: " << lastx << ", " << lasty << ", " << lastz << std::endl;
          std::cout << "slab loop" << std::endl;
          for(unsigned int slab = 0; slab < SlabN.size(); slab++)
            {
            double dist_slab = sqrt( pow(lastx - SlabX[slab], 2) +
                                     pow(lasty - SlabY[slab], 2) +
                                     pow(lastz - SlabZ[slab], 2) );
            std::cout<<"\tslab x,y,z: "<<SlabX[slab]<<", "<<SlabY[slab]<<", "<<SlabZ[slab]<<std::endl;
            std::cout << "\t\tdist to last: " << dist_slab << std::endl; 
            }
          
          }
        }
    std::cout<<std::endl;
    }//<---End checking if this is MC 

// ------------------------------------------------------   
//                      tof stuff 
// ------------------------------------------------------
    num_tof_objects = tof.size();
    size_t tof_counter = 0; // book-keeping
    for(size_t i = 0; i < tof.size(); i++){
      size_t number_tof = tof[i]->NTOF();
      for(size_t tof_idx = 0; tof_idx < number_tof; ++tof_idx){
        tofObject[tof_counter] =  tof[i]->SingleTOF(tof_idx);
        ++tof_counter;
      } // loop over TOF
    }//<---End tof_count loop

// --------------------------------------------------------
//                      wc track stuff
// --------------------------------------------------------
    num_wctracks = wctrack.size();
    for(size_t wct_count = 0; wct_count < wctrack.size(); wct_count++){
      wctrk_momentum[wct_count] = wctrack[wct_count]->Momentum();
      wctrk_XFace[wct_count] = wctrack[wct_count]->XYFace(0);
      wctrk_YFace[wct_count] = wctrack[wct_count]->XYFace(1);
      wctrk_theta[wct_count] = wctrack[wct_count]->Theta();
      wctrk_phi[wct_count] = wctrack[wct_count]->Phi();
      wctrk_missed[wct_count] = wctrack[wct_count]->WCMissed();
      wctrk_picky[wct_count] = static_cast< int > (wctrack[wct_count]->IsPicky());
    }


  // ----------------------------------------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------------------------------------
  //							FILLING THE 3-D TRACK INFORMATION
  // ----------------------------------------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------------------------------------
   
  // ### Calling the track momentum calculator ###
  trkf::TrackMomentumCalculator trkm;
  trkm.SetMinLength(10); //change the minimal track length requirement to 10 cm

  // === Saving the number of tracks per event ===
  ntracks_reco=tracklist.size();
  //TVector3 larStart;
  //TVector3 larEnd;
  
  // === Association between WC Tracks and TPC Tracks ===
  int TempTrackMatchedID = -1;
  if(evt.getByLabel(fWCTrackLabel, wctrackHandle))
    {
    std::cout<<"\t\t\tGOT LABEL\n";
    art::FindOneP<recob::Track> fWC2TPC(wctrackHandle, evt, fWC2TPCModuleLabel);
    if(fWC2TPC.isValid())
      {
      std::cout<<"\t\t\t\t\tmodule is valid!\n";
      // === Loop on all the Assn WC-TPC tracks === 
      for(unsigned int indexAssn = 0; indexAssn < fWC2TPC.size(); ++indexAssn )
        {
        // =========================                                                                                       
        // === Get the TPC track ===
        // =========================                                                                      
        cet::maybe_ref<recob::Track const> trackWC2TPC(*fWC2TPC.at(indexAssn));

        if(!trackWC2TPC) continue;
        recob::Track const& aTrack(trackWC2TPC.ref());
        TempTrackMatchedID = aTrack.ID();
        std::cout<<"\t\t\t\t\t\t\tTEMPTRACKID: "<<TempTrackMatchedID<<std::endl;

        }//<----End indexAssn loop          
      }//<---End checking that the WC2TPC   
    }



  // ### Looping over tracks ###
  //double maxtrackenergy = -1;



  //bool is_primary = false;
  //std::cout << "Reco info! trying to track down calo info..." << std::endl;
  for(int iTrack = 0; iTrack < ntracks_reco; iTrack++)
    {
    bool is_primary = false;

    // ### Storing an integer for the match of a WC to TPC track ###
    if(TempTrackMatchedID == tracklist[iTrack]->ID() )
      {track_WC2TPC_match.push_back(1);std::cout<<"\n\n\t\tMATCHED\n";}//<---End match
    else{track_WC2TPC_match.push_back(0);std::cout<<"\t\tNOTMATCHED\n";}


    
    // ### Setting the track information into memory ###
    auto trackStartEnd = tracklist[iTrack]->Extent();
    //larStart = tracklist[iTrack]->VertexDirection();
    //larEnd = tracklist[iTrack]->EndDirection();
    
    // ### track start and end ###
    track_start_x.push_back(trackStartEnd.first.X()); 
    track_start_y.push_back(trackStartEnd.first.Y()); 
    track_start_z.push_back(trackStartEnd.first.Z()); 
    track_end_x.push_back(trackStartEnd.second.X()); 
    track_end_y.push_back(trackStartEnd.second.Y()); 
    track_end_z.push_back(trackStartEnd.second.Z()); 

    
    // ### Recording the track length as calculated by the tracking module ###
    track_length.push_back(tracklist[iTrack]->Length());
    std::cout << "\ttrack length: " << tracklist[iTrack]->Length() << std::endl;
    
    // ### Grabbing the SpacePoints associated with this track ###
    std::vector< art::Ptr<recob::SpacePoint > > spts = fmsp.at(iTrack);
    ntrack_hits.push_back(fmsp.at(iTrack).size());
    std::cout << "\tntrack hits: " << spts.size() << std::endl;
    if(spts[0]->XYZ()[2] < 1){is_primary = true;}
    if(is_primary == true)
      {
      std::cout<<"\tprimary candidtate"<<std::endl;
      std::cout<<"\t\tfirst z: "<<spts[0]->XYZ()[2]<<std::endl;
      }
	if(is_primary==true){track_primary.push_back(1);}
	else                {track_primary.push_back(10);}

    // ### Looping over all the SpacePoints ###
    std::vector<double> x_spts;
    std::vector<double> y_spts;
    std::vector<double> z_spts;
    for(size_t jSpt = 0; jSpt < spts.size(); jSpt++)
      {
      x_spts.push_back(spts[jSpt]->XYZ()[0]);
      y_spts.push_back(spts[jSpt]->XYZ()[1]);
      z_spts.push_back(spts[jSpt]->XYZ()[2]);
      //if(is_primary == true)
      //  {
      //  std::cout << "\t\treco pt (x,y,z): " << spts[jSpt]->XYZ()[0] << ", "
      //                                       << spts[jSpt]->XYZ()[1] << ", "
      //                                       << spts[jSpt]->XYZ()[2] << ") " << std::endl;
      //  }
      }//<----End SpacePoint loop (j)
    track_xpos.push_back(x_spts);
    track_ypos.push_back(y_spts);
    track_zpos.push_back(z_spts);
    x_spts.clear(); y_spts.clear(); z_spts.clear();


    // ----------------------------------------------------------------------------------------------------------------------------
    // ----------------------------------------------------------------------------------------------------------------------------
    //							CALORIMERTY FROM THIS TRACK INFORMATION
    // ----------------------------------------------------------------------------------------------------------------------------
    // ----------------------------------------------------------------------------------------------------------------------------      
    
    // ########################################################## 
    // ### Looping over Calorimetry information for the track ###
    // ########################################################## 
    if(fmcal.isValid())
      {
      std::vector<art::Ptr<anab::Calorimetry> > calos = fmcal.at(iTrack);
      for (size_t j = 0; j<calos.size(); ++j)
        {
        if(!calos[j]->PlaneID().isValid){continue;}
        // ### Grabbing this calorimetry points plane number (0 == induction, 1 == collection) ###
        int pl = calos[j]->PlaneID().Plane;
        if(pl<0||pl>1){continue;}
    
        // ### Recording the number of calorimetry points for this track in this plane ####
        if(pl==0)
          {
          //std::cout << "        number of induction calorimetry points: " << calos[j]->dEdx().size() << std::endl;
          ind_track_hits.push_back(calos[j]->dEdx().size());
          ind_track_ke.push_back(calos[j]->KineticEnergy());
          //std::vector<double> ind_wire;
          std::vector<double> ind_dedx; std::vector<double> ind_dqdx;
          std::vector<double> ind_rr;   std::vector<double> ind_pitch_hit;
          for(size_t k = 0; k<calos[j]->dEdx().size(); ++k)
            {
            if(k>=1000){continue;}
            //ind_wire.push_back(calos[j]->dEdx()[k]);
            ind_dedx.push_back(calos[j]->dEdx()[k]);
            ind_dqdx.push_back(calos[j]->dQdx()[k]);
            ind_rr.push_back(calos[j]->ResidualRange()[k]);
            ind_pitch_hit.push_back(calos[j]->TrkPitchVec()[k]);
            }//<---End calo points (k)
          //ind_track_wire.push_back(ind_wire);
          ind_track_dedx.push_back(ind_dedx); ind_track_dqdx.push_back(ind_dqdx); 
          ind_track_rr.push_back(ind_dqdx);   ind_track_pitch_hit.push_back(ind_pitch_hit);
          ind_dedx.clear(); ind_dqdx.clear(); ind_rr.clear(); ind_pitch_hit.clear(); 
          }
        if(pl==1)
          {
          std::cout << "\t\t\tnumber of collection calorimetry points: " << calos[j]->dEdx().size() << std::endl;
          col_track_hits.push_back(calos[j]->dEdx().size());
          col_track_ke.push_back(calos[j]->KineticEnergy());
          //std::vector<double> col_wire;
          //std::vector<double> col_x;std::vector<double> col_y;std::vector<double> col_z;
          std::vector<double> col_dedx; std::vector<double> col_dqdx;
          std::vector<double> col_rr;   std::vector<double> col_pitch_hit;
          std::vector<double> col_x; std::vector<double> col_y; std::vector<double> col_z;
          for(size_t k = 0; k<calos[j]->dEdx().size(); ++k)
            {
            if(k>=1000){continue;}
            //std::cout << "\t\t\t\tpoint xyz: "<<calos[j]->XYZ()[k][0]<<", "<<calos[j]->XYZ()[k][1]<<", "<<calos[j]->XYZ()[k][0]<<std::endl;
            //std::cout << "\t\t\t\tpoint xyz: " << calos[j]->xyz().X() << ", " << calos[j]->xyz().Y() << ", " << calos[j]->xyz().Z() << std::endl;
            //Double_t col_x = calos[j]->XYZ()[0];
            //auto col_x = calos[j]->XYZ()(0);
            //auto col_x = calos[j]->XYZ().x();
            //auto col_x = calos[j]->XYZ()[k][0];
            //std::cout << "\t\t\t\tcolx: " << col_x << std::endl;
            //std::cout << "              calo pt wire: " << calos[j]->wire()[k] << std::endl;
            //col_wire.push_back(calos[j]->wire()[k]);
            //if(is_primary == true)
            //  {
            //  //std::cout << "\t\t\t\tpoint xyz: "<<calos[j]->XYZ()[k][0]<<", "<<calos[j]->XYZ()[k][1]<<", "<<calos[j]->XYZ()[k][2]<<std::endl;
            //  //double min_dist_to_true = 1.;
            //  //int closest_g4pt = 999;
            //  //for(int ng4 = 0; ng4 < geant_list_size; ng4++)
            //  //  {
            //  //  if(process_primary[ng4] == 1)
            //  //    {
            //  //    for(int trpt = 0; trpt < NTrTrajPts[ng4]; trpt++)
            //  //      {
            //  //      double dist_to_true = sqrt(pow(calos[j]->XYZ()[k][0] - MidPosX[ng4][trpt],2) 
            //  //                               + pow(calos[j]->XYZ()[k][1] - MidPosY[ng4][trpt],2) 
            //  //                               + pow(calos[j]->XYZ()[k][2] - MidPosZ[ng4][trpt],2)); 
            //  //      std::cout << "dist to true: " << dist_to_true << std::endl;
            //  //      if(dist_to_true < min_dist_to_true)
            //  //        {
            //  //        std::cout << "\t\t\t\t\tdistance: " << dist_to_true << std::endl;
            //  //        min_dist_to_true = dist_to_true;
            //  //        closest_g4pt = trpt;
            //  //        }//<--- End if close
            //  //      }//<--- End spt loop
            //  //    }//<--- End if primary
            //  //  }//<--- End g4 loop
            //  //std::cout << "\t\t\t\t\t\tclosest g4 pt: (pt, d): (" << closest_g4pt << ", " << min_dist_to_true << ")" << std::endl;
            //  }//<---End if primary
            col_x.push_back(calos[j]->XYZ()[k][0]);
            col_y.push_back(calos[j]->XYZ()[k][1]);
            col_z.push_back(calos[j]->XYZ()[k][2]);
            col_dedx.push_back(calos[j]->dEdx()[k]);
            col_dqdx.push_back(calos[j]->dQdx()[k]);
            col_rr.push_back(calos[j]->ResidualRange()[k]);
            col_pitch_hit.push_back(calos[j]->TrkPitchVec()[k]);
            }//<---End calo points (k)
          //col_track_wire.push_back(col_wire);
          col_track_dedx.push_back(col_dedx); col_track_dqdx.push_back(col_dqdx); 
          col_track_rr.push_back(col_dqdx);   col_track_pitch_hit.push_back(col_pitch_hit);
          col_track_x.push_back(col_x); col_track_y.push_back(col_y); col_track_z.push_back(col_z);
          col_dedx.clear(); col_dqdx.clear(); col_rr.clear(); col_pitch_hit.clear(); 
          col_x.clear(); col_y.clear(); col_z.clear();

          }

        }//<---End looping over calo points (j)
      }//<---End checking Calo info is valid  


      // ## Still in Track loop, want to try and use back tracker! ##
      // ## going to grab all the hits (really spts) and then use btracker ##
      // ## to get all the trackIDEs associated w/ this hit ##
      // ## should also do something like this but for clusters -- later ##


      if(!isdata&&fmth.isValid()) {
	    //int TrackID = 0;
	    std::vector< art::Ptr<recob::Hit> > allHits = fmth.at(iTrack);
        std::cout<<"\n\n\t\tBACKTRACKER STUDY\n";
        std::cout<<"allHits.size(): "<<allHits.size()<<std::endl;
          
        std::vector<int> n_ids;
        std::vector<std::map<int,double>> spt_g4map;
        std::vector<std::vector<int>> spt_id;
        std::vector<std::vector<double>> spt_e;
	    for(size_t h = 0; h < allHits.size(); ++h){
	      art::Ptr<recob::Hit> hit = allHits[h];
          std::vector<sim::TrackIDE> TrackIDs = bt_serv->HitToTrackIDEs(hit);
          std::cout<<"\th: "<<h<<std::endl;
          std::cout<<"\thit: "<<hit<<std::endl;
          std::cout<<"\tTrackIDs.size(): "<<TrackIDs.size()<<std::endl;

          n_ids.push_back(TrackIDs.size());
	      std::map<int,double> trkide;
	      std::vector<int> trkid;
	      std::vector<double> trke;
	      for(size_t e = 0; e < TrackIDs.size(); ++e){
            std::cout<<"\t\te: "<<e<<std::endl;
            std::cout<<"\t\t(trkid, e): "<<TrackIDs[e].trackID<<", "<<TrackIDs[e].energy<<std::endl;
	        trkide[TrackIDs[e].trackID] = TrackIDs[e].energy;
            trkid.push_back(TrackIDs[e].trackID);
            trke.push_back(TrackIDs[e].energy);
            track_spt_idarr[iTrack][h][e] = TrackIDs[e].trackID;
            track_spt_earr[iTrack][h][e] = TrackIDs[e].energy;
	      }	    
          // # push back a map to a vec where each index is a hit
        
          spt_g4map.push_back(trkide);
          spt_id.push_back(trkid);
          spt_e.push_back(trke);

          trkide.clear();
          trkid.clear();
          trke.clear();
	    }
        nhit_ids.push_back(n_ids);
        n_ids.clear();
        // # push back a vector of maps to a track vector
        track_spt_g4map.push_back(spt_g4map);
        track_spt_id.push_back(spt_id);
        track_spt_e.push_back(spt_e);
        spt_g4map.clear();
        spt_id.clear();
        spt_e.clear();
      }//<-- End if mc + fmth    



    }//<---End track loop (i)

    //// ### testing trackg4map ###
    //std::cout<<"\n\n\nTESTING MAP STUFF\n";
    //for(int ii = 0; ii < ntracks_reco; ++ii) {
    //  for(int jj = 0; jj < ntrack_hits[ii]; ++jj) {
    //    std::cout<<"\tspt#: "<<jj<<std::endl;
    //    std::cout<<"\t\t(x,y,z): "<<track_xpos[ii][jj]<<", "<<track_ypos[ii][jj]<<", "<<track_zpos[ii][jj]<<std::endl;
    //    std::cout<<"\t\tnum ids: "<<nhit_ids[ii][jj]<<std::endl;
    //    for(int kk = 0; kk < nhit_ids[ii][jj]; ++kk) {
    //      std::cout<<"\t\t\t(trackid, e): "<<track_spt_idarr[ii][jj][kk]<<", "<<track_spt_earr[ii][jj][kk]<<std::endl;
    //    }
    //    for(unsigned int ll = 0; ll < track_spt_id[ii][jj].size(); ++ll) {
    //      std::cout<<"\t\t\t(trackid, e): "<<track_spt_id[ii][jj][ll]<<", "<<track_spt_e[ii][jj][ll]<<std::endl;
    //    }
    //  }
    //} 


    // ### trying to get all the raw hits ? ###
    nhits = hitlist.size();
    std::cout<<"\n\nGoing to try and get Hit information for an event panel.\n";
    std::cout<<"hitlist.size(): "<<hitlist.size()<<std::endl;;
    for(int iHit = 0; iHit < nhits; iHit++) {
      std::cout<<"\thit #: "<<iHit<<std::endl;
      std::cout<<"\t\tt: "<<hitlist[iHit]->PeakTime()<<std::endl;
      std::cout<<"\t\ta: "<<hitlist[iHit]->PeakAmplitude()<<std::endl;
      std::cout<<"\t\tw: "<<hitlist[iHit]->Channel()<<std::endl;
      std::cout<<"\t\tv: "<<hitlist[iHit]->View()<<std::endl;
      hit_time.push_back(hitlist[iHit]->PeakTime());
      hit_wire.push_back(hitlist[iHit]->Channel());
      hit_view.push_back(hitlist[iHit]->View());
      hit_amp.push_back(hitlist[iHit]->PeakAmplitude());
    }
    


  fTree->Fill();

}//<---End analyze()


void lariat::AnaTreeT1034UC::beginJob()
{
   
  //std::cout<<"Check-1"<<std::endl;
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>(fTreeName.c_str(),  fTreeName.c_str());
  fTree->Branch("run",                         &run,"run/I");
  fTree->Branch("subrun",                      &subrun,"subrun/I");
  fTree->Branch("event",                       &event,"event/I");
  fTree->Branch("evttime",                     &evttime,"evttime/D");
  fTree->Branch("efield",                      efield,"efield[3]/D");
  fTree->Branch("t0",                          &t0,"t0/I");
  fTree->Branch("ntracks_reco",                &ntracks_reco,"ntracks_reco/I");
  fTree->Branch("track_WC2TPC_match",          &track_WC2TPC_match);
  fTree->Branch("track_primary",               &track_primary);
  fTree->Branch("track_start_x",               &track_start_x);
  fTree->Branch("track_start_y",               &track_start_y);
  fTree->Branch("track_start_z",               &track_start_z);
  fTree->Branch("track_end_x",                 &track_end_x);
  fTree->Branch("track_end_y",                 &track_end_y);
  fTree->Branch("track_end_z",                 &track_end_z);
  fTree->Branch("track_length",                &track_length);
  fTree->Branch("ntrack_hits",                 &ntrack_hits);
  fTree->Branch("track_xpos",                  &track_xpos);
  fTree->Branch("track_ypos",                  &track_ypos);
  fTree->Branch("track_zpos",                  &track_zpos);
  fTree->Branch("nhit_ids",                    &nhit_ids);

  //fTree->Branch("track_spt_g4map",             &track_spt_g4map);
  fTree->Branch("track_spt_idarr",             track_spt_idarr, "track_spt_idarr[ntracks_reco][1000][100]/D");
  fTree->Branch("track_spt_earr",              track_spt_earr, "track_spt_earr[ntracks_reco][1000][100]/D");
  
  fTree->Branch("ind_track_hits",              &ind_track_hits);
  fTree->Branch("ind_track_ke",                &ind_track_ke);
  fTree->Branch("ind_track_wire",              &ind_track_wire);
  fTree->Branch("ind_track_dedx",              &ind_track_dedx);
  fTree->Branch("ind_track_dqdx",              &ind_track_dqdx);
  fTree->Branch("ind_track_rr",                &ind_track_rr);
  fTree->Branch("ind_track_pitch_hit",         &ind_track_pitch_hit);
  fTree->Branch("col_track_hits",              &col_track_hits);
  fTree->Branch("col_track_ke",                &col_track_ke);
  fTree->Branch("col_track_x",                 &col_track_x);
  fTree->Branch("col_track_y",                 &col_track_y);
  fTree->Branch("col_track_z",                 &col_track_z);
  fTree->Branch("col_track_wire",              &col_track_wire);
  fTree->Branch("col_track_dedx",              &col_track_dedx);
  fTree->Branch("col_track_dqdx",              &col_track_dqdx);
  fTree->Branch("col_track_rr",                &col_track_rr);
  fTree->Branch("col_track_pitch_hit",         &col_track_pitch_hit);
  //fTree->Branch("num_hit",                     &num_hit,"num_hits/I");
  //fTree->Branch("hit_channel",                 &hit_channel);
  //fTree->Branch("hit_integral",                &hit_integral);

  fTree->Branch("nhits",                       &nhits,"nhits/I");
  fTree->Branch("hit_time",                    &hit_time);
  fTree->Branch("hit_amp",                     &hit_amp);
  fTree->Branch("hit_wire",                    &hit_wire);
  fTree->Branch("hit_view",                    &hit_view);
  
  fTree->Branch("no_primaries",                &no_primaries,"no_primaries/I");
  fTree->Branch("geant_list_size",             &geant_list_size,"geant_list_size/I");
  fTree->Branch("primary_p",                   &primary_p,"primary/D");
  fTree->Branch("PDG",                         &PDG);
  fTree->Branch("StartKE",                     &StartKE);
  fTree->Branch("LastKE",                      &LastKE);
  fTree->Branch("StartEnergy",                 &StartEnergy);
  fTree->Branch("StartPx",                     &StartPx);
  fTree->Branch("StartPy",                     &StartPy);
  fTree->Branch("StartPz",                     &StartPz);
  fTree->Branch("EndEnergy",                   &EndEnergy);
  fTree->Branch("EndPx",                       &EndPx);
  fTree->Branch("EndPy",                       &EndPy);
  fTree->Branch("EndPz",                       &EndPz);
  fTree->Branch("StartPointx",                 &StartPointx);
  fTree->Branch("StartPointy",                 &StartPointy);
  fTree->Branch("StartPointz",                 &StartPointz);
  fTree->Branch("EndPointx",                   &EndPointx);
  fTree->Branch("EndPointy",                   &EndPointy);
  fTree->Branch("EndPointz",                   &EndPointz);
  fTree->Branch("Process",                     &Process);
  fTree->Branch("NumberDaughters",             &NumberDaughters);
  //fTree->Branch("primary_simChannel_num_voxel",&primary_simChannel_num_voxel);
  //fTree->Branch("primary_simChannel_voxel_dr", &primary_simChannel_voxel_dr);
  //fTree->Branch("primary_simChannel_voxel_E",  &primary_simChannel_voxel_E);
  //fTree->Branch("primary_simChannel_voxel_e",  &primary_simChannel_voxel_e);
  //fTree->Branch("primary_simChannel_voxel_x",  &primary_simChannel_voxel_x);
  //fTree->Branch("primary_simChannel_voxel_y",  &primary_simChannel_voxel_y);
  //fTree->Branch("primary_simChannel_voxel_z",  &primary_simChannel_voxel_z);
  //fTree->Branch("primary_num_simChannel",      &primary_num_simChannel,"primary_num_simChannel/I");
  //fTree->Branch("primary_simChannel",          &primary_simChannel);
  //fTree->Branch("primary_simChannel_dr",       &primary_simChannel_dr);
  //fTree->Branch("primary_simChannel_E",        &primary_simChannel_E);
  //fTree->Branch("primary_simChannel_e",        &primary_simChannel_e);
  fTree->Branch("Mother",                      &Mother);
  fTree->Branch("TrackId",                     &TrackId);
  fTree->Branch("process_primary",             &process_primary);
  fTree->Branch("G4Process",                   &G4Process);
  fTree->Branch("G4FinalProcess",              &G4FinalProcess);  
  fTree->Branch("NTrTrajPts",                  &NTrTrajPts);
  fTree->Branch("NProtonDaughters",            &NProtonDaughters,"NProtonDaughters/I");
  fTree->Branch("NNeutronDaughters",           &NNeutronDaughters,"NNeutronDaughters/I");
  fTree->Branch("NDTrTrajPts",                 &NDTrTrajPts);
  fTree->Branch("DTrackId",                    &DTrackId);
  fTree->Branch("DPdgCode",                    &DPdgCode);
  fTree->Branch("DStartEnergy",                &DStartEnergy);
  fTree->Branch("DStartP",                     &DStartP);
  fTree->Branch("MidPosX",                     &MidPosX);
  fTree->Branch("MidPosY",                     &MidPosY);
  fTree->Branch("MidPosZ",                     &MidPosZ);
  fTree->Branch("MidPx",                       &MidPx);
  fTree->Branch("MidPy",                       &MidPy);
  fTree->Branch("MidPz",                       &MidPz);
  fTree->Branch("SlabX",                       &SlabX);
  fTree->Branch("SlabY",                       &SlabY);
  fTree->Branch("SlabZ",                       &SlabZ);
  fTree->Branch("SlapX",                       &SlapX);
  fTree->Branch("SlapY",                       &SlapY);
  fTree->Branch("SlapZ",                       &SlapZ);
  fTree->Branch("SlabN",                       &SlabN);
  fTree->Branch("SlabE",                       &SlabE);
  fTree->Branch("LastSlabInt",                 &LastSlabInt,"LastSlabInt/I");
  //fTree->Branch("PriMidPosX",                  &PriMidPosX);
  //fTree->Branch("PriMidPosY",                  &PriMidPosY);
  //fTree->Branch("PriMidPosZ",                  &PriMidPosZ);
  //fTree->Branch("PriMidPx",                    &PriMidPx);
  //fTree->Branch("PriMidPy",                    &PriMidPy);
  //fTree->Branch("PriMidPz",                    &PriMidPz);
  fTree->Branch("DMidPosX",                   &DMidPosX);
  fTree->Branch("DMidPosY",                   &DMidPosY);
  fTree->Branch("DMidPosZ",                   &DMidPosZ);
  fTree->Branch("InteractionPoint",            &InteractionPoint);
  fTree->Branch("InteractionPointType",        &InteractionPointType);

  fTree->Branch("num_tof_objects",             &num_tof_objects,"num_tof_objects/I");
  fTree->Branch("tofObject",                   tofObject,"tojObject[num_tof_objects]/D");
  fTree->Branch("num_wctracks",                &num_wctracks,"num_wctracks/I");
  fTree->Branch("wctrk_momentum",              wctrk_momentum,"wctrk_momentum[num_wctracks]/D");
  fTree->Branch("wctrk_XFace",                 wctrk_XFace,"wctrk_XFace[num_wctracks]/D");
  fTree->Branch("wctrk_YFace",                 wctrk_YFace,"wctrk_YFace[num_wctracks]/D");
  fTree->Branch("wctrk_theta",                 wctrk_theta,"wctrk_theta[num_wctracks]/D");
  fTree->Branch("wctrk_phi",                   wctrk_phi,"wctrk_phi[num_wctracks]/D");

  fTree->Branch("electron_lifetime",           &electron_lifetime, "electron_lifetime/D");

  // ### subdir for truth ###
  art::TFileDirectory truthDir = tfs->mkdir("truth");
  // ### histos! ###
  //E_vs_e   = truthDir.make<TH2D>("E_vs_e","IDE E(MeV) vs e(#)", 100, 0., 1., 600, 0., 6000.);

}

void lariat::AnaTreeT1034UC::ResetVars()
{
  G4Process.clear();
  G4FinalProcess.clear();
  InteractionPoint.clear();
  InteractionPointType.clear();
  PDG.clear();
  Mother.clear();
  TrackId.clear();
  NumberDaughters.clear();
  //primary_num_simChannel = -999;
  //primary_simChannel_num_voxel.clear();
  //primary_simChannel_voxel_dr.clear();
  //primary_simChannel_voxel_E.clear();
  //primary_simChannel_voxel_e.clear();
  //primary_simChannel_voxel_x.clear();
  //primary_simChannel_voxel_y.clear();
  //primary_simChannel_voxel_z.clear();
  //primary_simChannel.clear();
  //primary_simChannel_dr.clear();
  //primary_simChannel_E.clear();
  //primary_simChannel_e.clear();
  process_primary.clear();
  StartPointx.clear();  
  StartPointy.clear();
  StartPointz.clear();
  StartKE.clear();
  LastKE.clear();
  StartEnergy.clear();
  StartPx.clear();
  StartPy.clear();
  StartPz.clear(); 
  EndPointx.clear();
  EndPointy.clear(); 
  EndPointz.clear();  
  EndEnergy.clear();  
  EndPx.clear();  
  EndPy.clear();  
  EndPz.clear();  
  NDTrTrajPts.clear();
  DTrackId.clear();
  DPdgCode.clear();
  NTrTrajPts.clear();
  MidPosX.clear();
  MidPosY.clear();
  MidPosZ.clear();
  MidPx.clear();
  MidPy.clear();
  MidPz.clear();
  DMidPosX.clear();
  DMidPosY.clear();
  DMidPosZ.clear();
  DStartEnergy.clear();
  DStartP.clear();
  SlabX.clear();
  SlabY.clear();
  SlabZ.clear();
  SlapX.clear();
  SlapY.clear();
  SlapZ.clear();
  SlabN.clear();
  SlabE.clear();
  LastSlabInt = -99999;

  track_primary.clear();
  track_WC2TPC_match.clear();
  track_start_x.clear();
  track_start_y.clear();
  track_start_z.clear();
  track_end_x.clear();
  track_end_y.clear();
  track_end_z.clear();
  track_length.clear();
  ntrack_hits.clear();
  track_xpos.clear();
  track_ypos.clear();
  track_zpos.clear();
  nhit_ids.clear();
  track_spt_g4map.clear();
  track_spt_id.clear();
  track_spt_e.clear();
  ind_track_ke.clear();
  ind_track_wire.clear();
  ind_track_dedx.clear();
  ind_track_dqdx.clear();
  ind_track_rr.clear();
  ind_track_pitch_hit.clear();
  col_track_hits.clear();
  col_track_ke.clear();
  col_track_x.clear();
  col_track_y.clear();
  col_track_z.clear();
  col_track_wire.clear();
  col_track_dedx.clear();
  col_track_dqdx.clear();
  col_track_rr.clear();
  col_track_pitch_hit.clear();

  nhits = -99999;
  hit_time.clear();
  hit_wire.clear();
  hit_view.clear();
  hit_amp.clear();


  for(int i=0; i<kMaxTrack; ++i) {
    for(int j=0; j<kMaxTrackHits; ++j) {
      for(int k=0; k<kMaxHitIDs; ++k) {
        track_spt_idarr[i][j][k] = -9999;
        track_spt_earr[i][j][k] = -9999;
      }
    }
  }

  run = -99999;
  subrun = -99999;
  event = -99999;
  evttime = -99999;
  for (int i = 0; i<3; ++i){
    efield[i] = -99999;
  }
  t0 = -99999;
  NProtonDaughters = 0;
  NNeutronDaughters = 0;
  ntracks_reco = -99999;

  no_primaries = -99999;
  geant_list_size=-999;
  primary_p=-999;

  num_tof_objects = -99999;
  for(int i = 0; i < kMaxTOF; i++){
    tofObject[i] = -99999;
  }//<---End i loop

  num_wctracks = -9999999;
  for (int i = 0; i < kMaxWCTracks; i++){
    wctrk_momentum[i] = -999999;
    wctrk_XFace[i] = -999999;
    wctrk_YFace[i] = -999999;
    wctrk_theta[i] = -999999;
    wctrk_phi[i] = -999999;
    wctrk_missed[i] = -999999;
    wctrk_picky[i] = -999999;
  }
  
  electron_lifetime = -99999;

}

DEFINE_ART_MODULE(lariat::AnaTreeT1034UC)
