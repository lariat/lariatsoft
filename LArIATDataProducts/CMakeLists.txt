# for classes_def.xml!!
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

include_directories( $ENV{LARDATA_INC} )
include_directories( $ENV{LARCORE_INC} )

art_make(LIBRARY_NAME     LArIATDataProducts
      	 LIB_LIBRARIES    lardataobj_RawData
			                    ${ART_FRAMEWORK_CORE}
			                    ${ART_FRAMEWORK_PRINCIPAL}
			                    art_Persistency_Common canvas
			                    art_Persistency_Provenance canvas
			                    art_Utilities canvas
                          ${MF_MESSAGELOGGER}
			                    ${MF_UTILITIES}
                          cetlib
                          cetlib_except
			                    ${ROOT_BASIC_LIB_LIST}
         DICT_LIBRARIES   lardataobj_RawData
)


install_headers()
install_source()
