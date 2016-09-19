
include_directories( $ENV{LARDATA_INC}    )
include_directories( $ENV{LARCORE_INC}    )
include_directories( $ENV{LARSIM_INC}     )
include_directories( $ENV{POSTGRESQL_INC} )

art_make(BASENAME_ONLY

         LIB_LIBRARIES larcore_Geometry
                       larreco_RecoAlg_ClusterRecoUtil
                       lardataobj_RecoBase
                       LArIATChecks
         MODULE_LIBRARIES 
                          larevt_Filters
                          larsimobj_Simulation
                          larcore_Geometry
                          larcore_Geometry_Geometry_service
                          lardataobj_RecoBase
			  lardataobj_RawData
                          lardata_RecoBaseArt
                          lardata_Utilities
                          larcore_Geometry_AuxDetGeometry_service
			  SignalShapingServiceT1034_service
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
