
art_make( LIBRARY_NAME  LArIATRecoAlg
          LIB_LIBRARIES lardata_Utilities
          larcorealg_Geometry
          lardataobj_RawData
          Utilities_DatabaseUtilityT1034_service
          ${ART_FRAMEWORK_CORE}
          ${ART_FRAMEWORK_PRINCIPAL}
          ${ART_FRAMEWORK_SERVICES_REGISTRY}
          ${ART_FRAMEWORK_SERVICES_BASIC}
          ${ART_FRAMEWORK_SERVICES_OPTIONAL}
          ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
          ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
          art_Persistency_Common
          art_Persistency_Provenance
          art_Utilities
          canvas
          ${MF_MESSAGELOGGER}
          ${MF_UTILITIES}
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
