
include_directories( $ENV{ARTDAQ_CORE_INC}                     )
include_directories( $ENV{LARIATFRAGMENTS_INC}/LArIATFragments )
include_directories( $ENV{LARDATA_INC}                         )
include_directories( $ENV{LARCORE_INC}          	       )
include_directories( $ENV{LARSIM_INC}           	       )
include_directories( $ENV{POSTGRESQL_INC}       	       )

simple_plugin(LArIATFragmentReader "module"
              RawData
	      LArIATFragments
              ${ART_FRAMEWORK_CORE}
	      ${ART_FRAMEWORK_PRINCIPAL}
	      ${ART_FRAMEWORK_SERVICES_REGISTRY}
	      ${ART_FRAMEWORK_SERVICES_OPTIONAL}
	      ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
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
