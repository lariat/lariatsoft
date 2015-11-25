art_make(BASENAME_ONLY
         
	 LIB_LIBRARIES Geometry
	               ClusterRecoUtil
		       RecoBase
         MODULE_LIBRARIES LArIATRecoAlg
	 		  LArIATDataProducts
			  Filters
			  ClusterFinder
			  ClusterRecoUtil
			  RecoAlg
	 		  RawData
	                  SummaryData
			  Geometry
			  Geometry_service
			  RawDataUtilities
			  RecoBase
			  RecoBaseArt
			  Utilities
                          BackTracker_service
			  ${SEEDSERVICE_SERVICE} # artextensions
			  ${SIMULATIONBASE}
			  ${ART_FRAMEWORK_CORE}
			  ${ART_FRAMEWORK_PRINCIPAL}
			  ${ART_FRAMEWORK_SERVICES_REGISTRY}
			  ${ART_FRAMEWORK_SERVICES_BASIC}
			  ${ART_FRAMEWORK_SERVICES_OPTIONAL}
			  ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			  ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR}
			  ${ART_PERSISTENCY_COMMON}
			  ${ART_PERSISTENCY_PROVENANCE}
			  ${ART_UTILITIES}
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
