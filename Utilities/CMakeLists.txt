cet_enable_asserts()

set( lariat_util_lib_list lardata_Utilities_LArFFT_service
     larcorealg_Geometry
     larcore_Geometry_Geometry_service
     lardata_Utilities
     ${ART_FRAMEWORK_CORE}
     ${ART_FRAMEWORK_PRINCIPAL}
     ${ART_FRAMEWORK_SERVICES_REGISTRY}
     art_Persistency_Common canvas
     art_Persistency_Provenance canvas
     ${ART_FRAMEWORK_SERVICES_BASIC}
     ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
     ${ART_ROOT_IO_TFILESERVICE_SERVICE}
     ${ART_FRAMEWORK_SERVICES_SYSTEM_FILECATALOGMETADATA_SERVICE}
     art_Utilities canvas
     ${MF_MESSAGELOGGER}
     
     cetlib cetlib_except
     ${ROOT_BASIC_LIB_LIST}
     ${ROOT_GEOM}
)

simple_plugin( FileCatalogMetadataLArIAT "service"
               ${lariat_util_lib_list}
)

simple_plugin( SignalShapingServiceT1034  "service"
               ${lariat_util_lib_list}
)

simple_plugin( DatabaseUtilityT1034  "service"
               lardata_Utilities_DatabaseUtil_service
               ${MF_MESSAGELOGGER}
               
               ${FHICLCPP}
               cetlib cetlib_except
               ${PQ}
               art_Utilities canvas
               ${ART_FRAMEWORK_CORE}
               ${ART_FRAMEWORK_PRINCIPAL}
               ${ART_ROOT_IO_ROOTDB}
               art_Persistency_Common canvas
               art_Persistency_Provenance canvas
               ${ART_FRAMEWORK_SERVICES_REGISTRY}
)

simple_plugin( DetectorPropertiesServiceLArIAT  "service"
               LArIATDataProviders
               ${lariat_util_lib_list}
               ${SQLITE3}
               ${ART_ROOT_IO_ROOTDB}
)

file(GLOB pwd_files *passwd)
install( FILES ${pwd_files} DESTINATION ${product}/${version}/pwd/lardata_Utilities )

install_headers()
install_fhicl()
install_source()

