art_make(BASENAME_ONLY
         
	 LIB_LIBRARIES larcore_Geometry
	               larreco_RecoAlg_ClusterRecoUtil
		       lardata_RecoBase
		       LArIATChecks
         MODULE_LIBRARIES LArIATRecoAlg
	 		  LArIATDataProducts
			  larevt_Filters
			  larreco_ClusterFinder
			  larreco_RecoAlg_ClusterRecoUtil
			  larreco_RecoAlg
	 		  lardata_RawData
	                  larcore_SummaryData
			  larcore_Geometry
			  larcore_Geometry_Geometry_service
			  RawDataUtilities
			  lardata_RecoBase
			  lardata_RecoBaseArt
			  lardata_Utilities
			  
                          lardata_Utilities_DetectorProperties_service
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
