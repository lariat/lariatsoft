include_directories ( $ENV{GEANT4_FQ_DIR}/include )
include_directories ( $ENV{GEANT4_FQ_DIR}/lib64 )

art_make(LIB_LIBRARIES    larcorealg_Geometry
                          larreco_RecoAlg_ClusterRecoUtil
                          lardataobj_MCBase
                          lardataobj_RecoBase
                          LArIATChecks
         MODULE_LIBRARIES LArIATRecoAlg
                          LArIATDataProducts
                          larevt_Filters
                          larreco_ClusterFinder
                          larreco_RecoAlg_ClusterRecoUtil
                          larreco_RecoAlg
                          lardataobj_MCBase
                          lardataobj_RawData
                          lardataobj_RecoBase
                          larcoreobj_SummaryData
                          larcorealg_Geometry
                          larcore_Geometry_Geometry_service
                          RawDataUtilities
                          lardataobj_RecoBase
                          lardata_ArtDataHelper
                          lardata_Utilities
                          larsim_MCCheater_BackTrackerService_service
                          larsim_MCCheater_ParticleInventoryService_service
                          larcore_Geometry_AuxDetGeometry_service
                          nurandom_RandomUtils_NuRandomService_service
                          nusimdata_SimulationBase
                          ${ART_FRAMEWORK_CORE}
                          ${ART_FRAMEWORK_PRINCIPAL}
                          ${ART_FRAMEWORK_SERVICES_REGISTRY}
                          ${ART_FRAMEWORK_SERVICES_BASIC}
                          ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                          ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
	                        ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
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
	                        ${G4_LIB_LIST}
	                        nug4_G4Base
	                        larsim_LArG4
	                        ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
	                        ${ART_ROOT_IO_TFILESERVICE_SERVICE}
)

install_headers()
install_fhicl()
install_source()
