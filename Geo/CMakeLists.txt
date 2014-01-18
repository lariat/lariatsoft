
art_make( BASENAME_ONLY
          LIBRARY_NAME  Geo
          LIB_LIBRARIES Geometry
	  		${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_IO_SOURCES}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_PERSISTENCY_PROVENANCE}
			${ART_UTILITIES}
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
	                ${FHICLCPP}
			${CETLIB}
			${ROOT_BASIC_LIB_LIST}
                        ${ROOT_GEOM}
                        ${ROOT_XMLIO}
                        ${ROOT_GDML}
                        ${ROOT_EG}
          SERVICE_LIBRARIES Geo
	  		    Geometry
			    ${ART_FRAMEWORK_CORE}
			    ${ART_FRAMEWORK_IO_SOURCES}
			    ${ART_FRAMEWORK_PRINCIPAL}
			    ${ART_PERSISTENCY_COMMON}
			    ${ART_PERSISTENCY_PROVENANCE}
			    ${ART_UTILITIES}
			    ${ART_FRAMEWORK_SERVICES_REGISTRY}
			    ${MF_MESSAGELOGGER}
			    ${MF_UTILITIES}
	                    ${FHICLCPP}
			    ${CETLIB}
          )

install_headers()
install_fhicl()
install_source()

# testing
add_subdirectory(gdml)
