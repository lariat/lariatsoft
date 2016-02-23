art_make_library( LIBRARY_NAME Geo
                  SOURCE     AuxDetChannelMapLArIATAlg.cxx 
			     AuxDetGeoObjectSorterLArIAT.cxx  
		             ChannelMapLArIATAlg.cxx 
			     GeoObjectSorterLArIAT.cxx 
		  LIBRARIES larcore_Geometry
		            ${MF_MESSAGELOGGER}
		            ${MF_UTILITIES}
			    ${FHICLCPP}
			    ${ROOT_BASIC_LIB_LIST} 
			    ${ROOT_GEOM}
			    ${CETLIB}
)

simple_plugin( LArIATGeometryHelper service
	  	larcore_Geometry
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

simple_plugin( LArIATAuxDetGeometryHelper service
	  	larcore_Geometry
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
