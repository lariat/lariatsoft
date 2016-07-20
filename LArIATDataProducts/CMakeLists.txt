# for classes_def.xml!!
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

include_directories( $ENV{LARDATA_INC}                         )
include_directories( $ENV{LARCORE_INC}          	       )

art_make(BASENAME_ONLY
         LIBRARY_NAME     LArIATDataProducts
	 LIB_LIBRARIES    lardataobj_RawData
			  ${ART_FRAMEWORK_CORE}
			  ${ART_FRAMEWORK_PRINCIPAL}
			  art_Persistency_Common canvas_Persistency_Common
			  art_Persistency_Provenance canvas_Persistency_Provenance
			  art_Utilities canvas_Utilities
			  ${MF_MESSAGELOGGER}
			  ${MF_UTILITIES}
			  ${CETLIB}
			  ${ROOT_BASIC_LIB_LIST}
         DICT_LIBRARIES   lardataobj_RawData
)


install_headers()
install_source()
