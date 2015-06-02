# for classes_def.xml!!
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

include_directories( $ENV{LARDATA_INC}                         )
include_directories( $ENV{LARCORE_INC}          	       )

art_make(BASENAME_ONLY
         LIBRARY_NAME     LArIATDataProducts
	 LIB_LIBRARIES    RawData
			  ${ART_FRAMEWORK_CORE}
			  ${ART_FRAMEWORK_PRINCIPAL}
			  ${ART_PERSISTENCY_COMMON}
			  ${ART_PERSISTENCY_PROVENANCE}
			  ${ART_UTILITIES}
			  ${MF_MESSAGELOGGER}
			  ${MF_UTILITIES}
			  ${CETLIB}
			  ${ROOT_BASIC_LIB_LIST}
         DICT_LIBRARIES   RawData
)


install_headers()
install_source()
