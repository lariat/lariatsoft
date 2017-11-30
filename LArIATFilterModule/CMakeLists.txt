art_make(LIB_LIBRARIES larcorealg_Geometry
         larreco_RecoAlg_ClusterRecoUtil
         lardataobj_RecoBase
         MODULE_LIBRARIES LArIATRecoAlg
         larreco_RecoAlg_ImagePatternAlgs_PointIdAlg
         LArIATDataProducts
         larevt_Filters
         larreco_ClusterFinder
         larreco_RecoAlg_ClusterRecoUtil
         Utilities_DatabaseUtilityT1034_service
         larreco_RecoAlg
         lardataobj_RawData
         larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
         larcoreobj_SummaryData
         larcorealg_Geometry
         larcore_Geometry_Geometry_service
         RawDataUtilities
         lardataobj_RecoBase
         lardata_ArtDataHelper
         lardata_RecoObjects larreco_Deprecated
         lardata_Utilities
         larsim_MCCheater_BackTrackerService_service
         larsim_MCCheater_ParticleInventoryService_service
         nutools_RandomUtils_NuRandomService_service
         nusimdata_SimulationBase
         ${ART_FRAMEWORK_CORE}
         ${ART_FRAMEWORK_PRINCIPAL}
         ${ART_FRAMEWORK_SERVICES_REGISTRY}
         ${ART_FRAMEWORK_SERVICES_BASIC}
         ${ART_FRAMEWORK_SERVICES_OPTIONAL}
         ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
         ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR}
         art_Persistency_Common canvas
         art_Persistency_Provenance canvas
         art_Utilities canvas
         ${MF_MESSAGELOGGER}
         ${MF_UTILITIES}
         ${FHICLCPP}
         cetlib cetlib_except
         ${CLHEP}
         ${ROOT_BASIC_LIB_LIST}
)

install_headers()
install_fhicl()
install_source()
