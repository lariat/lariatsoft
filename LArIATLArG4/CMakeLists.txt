
include_directories ( $ENV{GENIE_INC}/GENIE )
include_directories ( $ENV{LIBXML2_FQ_DIR}/include/libxml2 )
include_directories ( $ENV{GEANT4_FQ_DIR}/include )
include_directories ( $ENV{XERCES_C_INC} )

art_make( 
          LIBRARY_NAME     LArIATLArG4
          LIB_LIBRARIES    larsim_PhotonPropagation
                           larsim_PhotonPropagation_PhotonVisibilityService_service
                           lardata_Utilities
                           larsim_LArG4
                           larsim_Simulation 
                           nug4_ParticleNavigation
                           lardataobj_Simulation
                           larcorealg_Geometry
                           larcore_Geometry_Geometry_service
                           nusimdata_SimulationBase
                           nug4_G4Base
                           ${ART_FRAMEWORK_CORE}
                           ${ART_FRAMEWORK_PRINCIPAL}
                           ${ART_FRAMEWORK_SERVICES_REGISTRY}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL}
#                           ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                           ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                           ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                           art_Persistency_Common
                           art_Persistency_Provenance
                           art_Utilities
                           canvas
                           ${MF_MESSAGELOGGER}
                           ${MF_UTILITIES}
                           ${FHICLCPP}
                           cetlib cetlib_except
                           ${ROOT_GENVECTOR}
                           ${ROOT_EG}
                           ${ROOT_TREEPLAYER} 
                           ${ROOT_FFTW}
                           ${ROOT_REFLEX}
                           ${ROOT_EGPYTHIA6}
                           ${ROOT_GUI}
                           ${ROOT_BASIC_LIB_LIST}
                           ${ROOT_GEOM}
                           ${CLHEP}
                           ${G4_LIB_LIST}
                           ${Boost_SYSTEM_LIBRARY}
          MODULE_LIBRARIES LArIATLArG4
                           larsim_LArG4
                           larsim_MCCheater_ParticleInventoryService_service
                           larsim_Simulation
                           nug4_ParticleNavigation
                           lardataobj_Simulation
                           larcorealg_Geometry
                           larcore_Geometry_Geometry_service
                           nurandom_RandomUtils_NuRandomService_service
                           nusimdata_SimulationBase
                           ${ART_FRAMEWORK_SERVICES_REGISTRY}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                           #${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                           ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                           ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                           canvas
                           ${MF_MESSAGELOGGER}
                           ${MF_UTILITIES}
                           ${FHICLCPP}
                           cetlib cetlib_except
                           ${ROOT_GENVECTOR}
                           ${ROOT_BASIC_LIB_LIST}
                           ${ROOT_GEOM}
                           nug4_G4Base
                           ${CLHEP}
                           ${G4_LIB_LIST}
)

install_headers()
install_fhicl()
install_source()

