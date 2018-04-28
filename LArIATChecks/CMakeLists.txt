#
art_make(
          EXCLUDE NNTrackShowerEff_module.cc NNMCTrackShowerEff_module.cc NNPlotter_module.cc
  LIB_LIBRARIES    larcorealg_Geometry
  larreco_RecoAlg_ClusterRecoUtil
  lardataobj_RecoBase
  MODULE_LIBRARIES LArIATRecoAlg
  LArIATDataProducts
  larevt_Filters
  larreco_ClusterFinder
  larreco_RecoAlg_ClusterRecoUtil
  larreco_RecoAlg
  lardataobj_RawData
  larcoreobj_SummaryData
  larcorealg_Geometry
  larcore_Geometry_Geometry_service
  RawDataUtilities
  lardataobj_RecoBase
  lardata_ArtDataHelper
  lardata_Utilities
  larsim_MCCheater_BackTrackerService_service
  larsim_MCCheater_ParticleInventoryService_service
  nutools_RandomUtils_NuRandomService_service
  larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
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
  ${TBB}
  )

# tensorflow is now optional
if (DEFINED ENV{TENSORFLOW_DIR})
simple_plugin(NNTrackShowerEff "module"
  LArIATRecoAlg
  larreco_RecoAlg_ImagePatternAlgs_DataProvider
  larreco_RecoAlg_ImagePatternAlgs_Tensorflow_TF
  larreco_RecoAlg_ImagePatternAlgs_Tensorflow_PointIdAlg
  LArIATDataProducts
  larevt_Filters
  larreco_ClusterFinder
  larreco_RecoAlg_ClusterRecoUtil
  larreco_RecoAlg
  lardataobj_RawData
  larcoreobj_SummaryData
  larcorealg_Geometry
  larcore_Geometry_Geometry_service
  RawDataUtilities
  lardataobj_RecoBase
  lardata_ArtDataHelper
  lardata_Utilities
  larsim_MCCheater_BackTrackerService_service
  larsim_MCCheater_ParticleInventoryService_service
  nutools_RandomUtils_NuRandomService_service
  larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
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
  ${TBB}
)

simple_plugin(NNMCTrackShowerEff "module"
  LArIATRecoAlg
  larreco_RecoAlg_ImagePatternAlgs_DataProvider
  larreco_RecoAlg_ImagePatternAlgs_Tensorflow_TF
  larreco_RecoAlg_ImagePatternAlgs_Tensorflow_PointIdAlg
  LArIATDataProducts
  larevt_Filters
  larreco_ClusterFinder
  larreco_RecoAlg_ClusterRecoUtil
  larreco_RecoAlg
  lardataobj_RawData
  larcoreobj_SummaryData
  larcorealg_Geometry
  larcore_Geometry_Geometry_service
  RawDataUtilities
  lardataobj_RecoBase
  lardata_ArtDataHelper
  lardata_Utilities
  larsim_MCCheater_BackTrackerService_service
  larsim_MCCheater_ParticleInventoryService_service
  nutools_RandomUtils_NuRandomService_service
  larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
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
  ${TBB}
)

simple_plugin(NNPlotter "module"
  LArIATRecoAlg
  larreco_RecoAlg_ImagePatternAlgs_DataProvider
  larreco_RecoAlg_ImagePatternAlgs_Tensorflow_TF
  larreco_RecoAlg_ImagePatternAlgs_Tensorflow_PointIdAlg
  LArIATDataProducts
  larevt_Filters
  larreco_ClusterFinder
  larreco_RecoAlg_ClusterRecoUtil
  larreco_RecoAlg
  lardataobj_RawData
  larcoreobj_SummaryData
  larcorealg_Geometry
  larcore_Geometry_Geometry_service
  RawDataUtilities
  lardataobj_RecoBase
  lardata_ArtDataHelper
  lardata_Utilities
  larsim_MCCheater_BackTrackerService_service
  larsim_MCCheater_ParticleInventoryService_service
  nutools_RandomUtils_NuRandomService_service
  larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
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
  ${TBB}
)
endif (DEFINED ENV{TENSORFLOW_DIR})

install_headers()
install_fhicl()
install_source()
