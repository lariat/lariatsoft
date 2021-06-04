include_directories ( $ENV{GEANT4_FQ_DIR}/include )
include_directories ( $ENV{GEANT4_FQ_DIR}/lib64 )

art_make( LIBRARY_NAME  LArIATRecoAlg
          LIB_LIBRARIES lardata_Utilities
          larcorealg_Geometry
          lardataobj_RawData
          Utilities_DatabaseUtilityT1034_service
          LArIATDataProducts
	        nug4_G4Base
	        larsim_LArG4
	        nusimdata_SimulationBase
	        larsim_Simulation nug4_ParticleNavigation lardataobj_Simulation
	        ${G4_LIB_LIST}
          ${ART_FRAMEWORK_CORE}
          ${ART_FRAMEWORK_PRINCIPAL}
          ${ART_FRAMEWORK_SERVICES_REGISTRY}
          ${ART_FRAMEWORK_SERVICES_BASIC}
          ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
          ${ART_ROOT_IO_TFILESERVICE_SERVICE}
          ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
          art_Persistency_Common
          art_Persistency_Provenance
          art_Utilities
          canvas
          ${MF_MESSAGELOGGER}
          
          ${FHICLCPP}
          ${CETLIB}
          cetlib cetlib_except
          ${CLHEP}
          ${ROOT_BASIC_LIB_LIST}
          ${ROOT_GEOM}
          ${ROOT_MINUIT}
          ${ROOT_MINUIT2}
          ${Boost_SYSTEM_LIBRARY}
        )

install_headers()
install_fhicl()
install_source()
