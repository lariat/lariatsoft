
include_directories( $ENV{LARDATA_INC}    )
include_directories( $ENV{LARCORE_INC}    )
include_directories( $ENV{LARSIM_INC}     )
include_directories( $ENV{POSTGRESQL_INC} )

simple_plugin(SimWireLArIAT "module"
			Simulation
	                RawData
	                Geometry
	                Geometry_service
			LArFFT_service
			LArProperties_service
			DetectorProperties_service
			Utilities
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
			${ART_PERSISTENCY_COMMON}
			${ART_PERSISTENCY_PROVENANCE}
			${ART_UTILITIES}
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
			${FHICLCPP}
			${CETLIB}
			${CLHEP}
			${ROOT_BASIC_LIB_LIST}
              BASENAME_ONLY
)

install_headers()
install_fhicl()
install_source()
