art_make(LIB_LIBRARIES larcorealg_Geometry
         larreco_RecoAlg_ClusterRecoUtil
         lardataobj_RecoBase
         ${ROOT_BASIC_LIB_LIST}
         ${ROOT_SPECTRUM}
         MODULE_LIBRARIES LArIATRecoAlg
         LArIATDataProducts
	       LArIATAnaModule
         larevt_Filters
         larreco_ClusterFinder
         larreco_RecoAlg_ClusterRecoUtil
         Utilities_DatabaseUtilityT1034_service
         larreco_RecoAlg
	       larreco_Calorimetry
         lardataobj_RawData
         larsim_Simulation nug4_ParticleNavigation lardataobj_Simulation
         larcoreobj_SummaryData
         larcorealg_Geometry
         larcore_Geometry_Geometry_service
         RawDataUtilities
         lardata_ArtDataHelper
         lardata_RecoObjects
         lardata_Utilities
         
         
         larsim_MCCheater_BackTrackerService_service
         larsim_MCCheater_ParticleInventoryService_service
         nurandom_RandomUtils_NuRandomService_service
         nusimdata_SimulationBase
         ${ART_FRAMEWORK_CORE}
         ${ART_FRAMEWORK_PRINCIPAL}
         ${ART_FRAMEWORK_SERVICES_REGISTRY}
         ${ART_FRAMEWORK_SERVICES_BASIC}
         ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
         ${ART_ROOT_IO_TFILESERVICE_SERVICE}
         ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR}
         art_Persistency_Common canvas
         art_Persistency_Provenance canvas
         art_Utilities canvas
         ${MF_MESSAGELOGGER}
         
         ${FHICLCPP}
         cetlib cetlib_except
         ${CLHEP}
         ${ROOT_BASIC_LIB_LIST}
	
)

install_headers()
install_fhicl()
install_source()
