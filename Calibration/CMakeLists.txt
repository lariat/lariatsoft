
add_subdirectory(constants)

simple_plugin(CalWireT1034 "module"
			        lardata_Utilities
              larevt_Filters
			        lardataobj_RawData
              RawDataUtilities
              lardataobj_RecoBase
              lardata_ArtDataHelper
			        RawDataUtilities
		          SignalShapingServiceT1034_service
              ${ART_FRAMEWORK_CORE}
			        ${ART_FRAMEWORK_PRINCIPAL}
			        ${ART_FRAMEWORK_SERVICES_REGISTRY}
			        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
			        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			        art_Persistency_Common canvas
			        art_Persistency_Provenance canvas
			        art_Utilities canvas
			        ${MF_MESSAGELOGGER}
			        ${MF_UTILITIES}
              ${FHICLCPP}
			        cetlib cetlib_except
              ${ROOT_GEOM}
              ${ROOT_XMLIO}
              ${ROOT_GDML}
			        ${ROOT_BASIC_LIB_LIST}
              BASENAME_ONLY
)

simple_plugin(CalWireROIT1034 "module"
			        lardata_Utilities
              larevt_Filters
			        lardataobj_RawData
              lardataobj_RecoBase
              lardata_ArtDataHelper
			        RawDataUtilities
		          SignalShapingServiceT1034_service
              ${ART_FRAMEWORK_CORE}
			        ${ART_FRAMEWORK_PRINCIPAL}
              ${ART_FRAMEWORK_SERVICES_REGISTRY}
			        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
			        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			        art_Persistency_Common canvas
			        art_Persistency_Provenance canvas
			        art_Utilities canvas
			        ${MF_MESSAGELOGGER}
			        ${MF_UTILITIES}
              ${FHICLCPP}
			        cetlib cetlib_except
              ${ROOT_GEOM}
              ${ROOT_XMLIO}
              ${ROOT_GDML}
			        ${ROOT_BASIC_LIB_LIST}
              BASENAME_ONLY
)

install_headers()
install_fhicl()
install_source()
