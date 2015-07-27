# for classes_def.xml!!
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

include_directories( $ENV{ARTDAQ_CORE_INC}                     )
include_directories( $ENV{LARIATFRAGMENTS_INC}/LArIATFragments )
include_directories( $ENV{LARDATA_INC}                         )
include_directories( $ENV{LARCORE_INC}          	       )
include_directories( $ENV{LARSIM_INC}           	       )
include_directories( $ENV{POSTGRESQL_INC}       	       )

art_make(BASENAME_ONLY
         LIBRARY_NAME     RawDataUtilities
	 LIB_LIBRARIES    RawData
	                  LariatFragments
			  ${ART_FRAMEWORK_CORE}
			  ${ART_FRAMEWORK_PRINCIPAL}
			  ${ART_PERSISTENCY_COMMON}
			  ${ART_PERSISTENCY_PROVENANCE}
			  ${ART_UTILITIES}
			  ${MF_MESSAGELOGGER}
			  ${MF_UTILITIES}
			  ${CETLIB}
			  ${ROOT_BASIC_LIB_LIST}
			  ${Boost_SYSTEM_LIBRARY}	
	                  ${ART_FRAMEWORK_SERVICES_REGISTRY}
			  ${ART_FRAMEWORK_SERVICES_BASIC}
			  ${ART_FRAMEWORK_SERVICES_OPTIONAL}
			  ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                          ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                          ${FHICLCPP}	  
         DICT_LIBRARIES   RawData
	 		  SummaryData
	 SOURCE_LIBRARIES RawData
	                  RawDataUtilities
	 		  LArIATFragments
			  SummaryData
	                  ${ART_FRAMEWORK_SERVICES_REGISTRY}
			  ${ART_FRAMEWORK_SERVICES_BASIC}
			  ${ART_FRAMEWORK_SERVICES_OPTIONAL}
			  ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                          ${ART_FRAMEWORK_IO_ROOT}
			  ${MF_MESSAGELOGGER}
                          ${ROOT_BASIC_LIB_LIST}
         MODULE_LIBRARIES RawData
	                  SummaryData
			  Geometry_service
                          LariatFragments
			  RawDataUtilities
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
)


install_headers()
install_fhicl()
install_source()
