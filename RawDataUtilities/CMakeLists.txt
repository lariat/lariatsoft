# for classes_def.xml!!
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

include_directories( $ENV{ARTDAQ_CORE_INC}                     )
include_directories( $ENV{LARIATFRAGMENTS_INC}/LArIATFragments )
include_directories( $ENV{LARDATA_INC}                         )
include_directories( $ENV{LARCORE_INC}                         )
include_directories( $ENV{LARSIM_INC}                          )
include_directories( $ENV{POSTGRESQL_INC}                      )

art_make(LIBRARY_NAME     RawDataUtilities
         LIB_LIBRARIES    lardataobj_RawData
                          LariatFragments
                          Utilities_DatabaseUtilityT1034_service
                          ${ART_FRAMEWORK_CORE}
                          ${ART_FRAMEWORK_PRINCIPAL}
                          art_Persistency_Common canvas
                          art_Persistency_Provenance canvas
                          art_Utilities canvas
                          ${MF_MESSAGELOGGER}
                          ${ARTDAQ-CORE_UTILITIES}
                          cetlib cetlib_except
                          ${ROOT_BASIC_LIB_LIST}
                          ${Boost_SYSTEM_LIBRARY}
                          ${ART_FRAMEWORK_SERVICES_REGISTRY}
                          ${ART_FRAMEWORK_SERVICES_BASIC}
                          ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                          ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                          ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                          ${FHICLCPP}
         DICT_LIBRARIES   lardataobj_RawData
                          larcoreobj_SummaryData
         SOURCE_LIBRARIES larcorealg_Geometry
                          lardataobj_RawData
                          lardata_Utilities
                          LariatFragments
                          LArIATDataProducts
                          larcoreobj_SummaryData
                          RawDataUtilities
                          Utilities_DatabaseUtilityT1034_service
                          ${ART_FRAMEWORK_SERVICES_REGISTRY}
                          ${ART_FRAMEWORK_SERVICES_BASIC}
                          ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                          ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                          ${ART_ROOT_IO}
                          ${MF_MESSAGELOGGER}
                          ${ROOT_BASIC_LIB_LIST}
         MODULE_LIBRARIES lardataobj_RawData
                          larcoreobj_SummaryData
                          larcore_Geometry_Geometry_service
                          LariatFragments
                          RawDataUtilities
                          ${ART_FRAMEWORK_CORE}
                          ${ART_FRAMEWORK_PRINCIPAL}
                          ${ART_FRAMEWORK_SERVICES_REGISTRY}
                          ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                          ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                          art_Persistency_Common canvas
                          art_Persistency_Provenance canvas
                          art_Utilities canvas
                          ${MF_MESSAGELOGGER}
                          
                          ${FHICLCPP}
                          cetlib cetlib_except
                          ${CLHEP}
                          ${ROOT_BASIC_LIB_LIST}
)


install_headers()
install_fhicl()
install_source()
