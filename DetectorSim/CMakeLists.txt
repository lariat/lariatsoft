
include_directories( $ENV{LARDATA_INC}    )
include_directories( $ENV{LARCORE_INC}    )
include_directories( $ENV{LARSIM_INC}     )
include_directories( $ENV{POSTGRESQL_INC} )

art_make(LIB_LIBRARIES larcore_Geometry
                       larreco_RecoAlg_ClusterRecoUtil
                       lardataobj_RecoBase
                       LArIATChecks
         MODULE_LIBRARIES 
                          larevt_Filters
                          lardataobj_Simulation
                          larcore_Geometry
                          larcore_Geometry_Geometry_service
                          lardataobj_RecoBase
			                    lardataobj_RawData
                          lardata_ArtDataHelper
                          lardata_Utilities
                          larcore_Geometry_AuxDetGeometry_service
			                    Utilities_SignalShapingServiceT1034_service
                          nutools_RandomUtils_NuRandomService_service
			                    lardata_DetectorInfoServices_DetectorClocksServiceStandard_service
                          nusimdata_SimulationBase
                          ${ART_FRAMEWORK_CORE}
                          ${ART_FRAMEWORK_PRINCIPAL}
                          ${ART_FRAMEWORK_SERVICES_REGISTRY}
                          ${ART_FRAMEWORK_SERVICES_BASIC}
                          ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                          ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                          ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR}
			                    ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
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
