art_make(BASENAME_ONLY
         LIB_LIBRARIES larcore_Geometry
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
         DatabaseUtilityT1034_service
         larreco_RecoAlg
	 lardata_AnalysisAlg
         lardataobj_RawData
         larsim_Simulation lardataobj_Simulation
         larcoreobj_SummaryData
         larcore_Geometry
         larcore_Geometry_Geometry_service
         RawDataUtilities
         lardataobj_RecoBase
         lardata_RecoBaseArt
         lardata_RecoObjects
         lardata_Utilities
         
         
         larsim_MCCheater_BackTracker_service
         larsim_RandomUtils_LArSeedService_service
         nusimdata_SimulationBase
         ${ART_FRAMEWORK_CORE}
         ${ART_FRAMEWORK_PRINCIPAL}
         ${ART_FRAMEWORK_SERVICES_REGISTRY}
         ${ART_FRAMEWORK_SERVICES_BASIC}
         ${ART_FRAMEWORK_SERVICES_OPTIONAL}
         ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
         ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR}
         art_Persistency_Common canvas_Persistency_Common
         art_Persistency_Provenance canvas_Persistency_Provenance
         art_Utilities canvas_Utilities
         ${MF_MESSAGELOGGER}
         ${MF_UTILITIES}
         ${FHICLCPP}
         ${CETLIB}
         ${CLHEP}
         ${ROOT_BASIC_LIB_LIST}
	
)

install_headers()
install_fhicl()
install_source()
