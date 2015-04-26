art_make_library( LIBRARY_NAME Geo
                  SOURCE ChannelMapLArIATAlg.cxx GeoObjectSorterLArIAT.cxx
		  LIBRARIES Geometry
		            ${MF_MESSAGELOGGER}
		            ${MF_UTILITIES}
			    ${FHICLCPP}
			    ${ROOT_BASIC_LIB_LIST}                
			    ${CETLIB}
)

simple_plugin( LArIATGeometryHelper service
	  	Geometry
		Geo
		${ART_FRAMEWORK_CORE}
		${ART_FRAMEWORK_IO_SOURCES}
		${ART_FRAMEWORK_PRINCIPAL}
		${ART_PERSISTENCY_COMMON}
		${ART_PERSISTENCY_PROVENANCE}
		${ART_UTILITIES}
		${ART_FRAMEWORK_SERVICES_REGISTRY}
		${MF_MESSAGELOGGER}
		${MF_UTILITIES}
		${ROOT_BASIC_LIB_LIST}                
		${CETLIB}
		BASENAME_ONLY
)

install_headers()
install_fhicl()
install_source()

# testing
add_subdirectory(gdml)
