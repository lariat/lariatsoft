cet_enable_asserts()

set( lariat_util_lib_list lardata_Utilities_LArFFT_service
     larcore_Geometry
     larcore_Geometry_Geometry_service
     lardata_Utilities
     ${ART_FRAMEWORK_CORE}
     ${ART_FRAMEWORK_PRINCIPAL}
     ${ART_FRAMEWORK_SERVICES_REGISTRY}
     art_Persistency_Common canvas
     art_Persistency_Provenance canvas
     ${ART_FRAMEWORK_SERVICES_BASIC}
     ${ART_FRAMEWORK_SERVICES_OPTIONAL}
     ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
     ${ART_FRAMEWORK_SERVICES_SYSTEM_FILECATALOGMETADATA_SERVICE}
     art_Utilities canvas
     ${MF_MESSAGELOGGER}
     ${MF_UTILITIES}
     cetlib cetlib_except
     ${ROOT_BASIC_LIB_LIST}
     ${ROOT_GEOM}
)

simple_plugin( FileCatalogMetadataLArIAT "service"
               ${lariat_util_lib_list}
               BASENAME_ONLY
)

simple_plugin( SignalShapingServiceT1034  "service"
               ${lariat_util_lib_list}
               BASENAME_ONLY
)

simple_plugin( DatabaseUtilityT1034  "service"
               lardata_Utilities_DatabaseUtil_service
               ${MF_MESSAGELOGGER}
               ${MF_UTILITIES}
               ${FHICLCPP}
               cetlib cetlib_except
               ${PQ}
               art_Utilities canvas
               ${ART_FRAMEWORK_CORE}
               ${ART_FRAMEWORK_PRINCIPAL}
               ${ART_PERSISTENCY_ROOTDB}
               art_Persistency_Common canvas
               art_Persistency_Provenance canvas
               ${ART_FRAMEWORK_SERVICES_REGISTRY}
               BASENAME_ONLY
)

simple_plugin( DetectorPropertiesServiceLArIAT  "service"
               LArIATDataProviders
               ${lariat_util_lib_list}
               ${SQLITE3}
               ${ART_PERSISTENCY_ROOTDB}
               BASENAME_ONLY
)

file(GLOB pwd_files *passwd)
install( FILES ${pwd_files} DESTINATION ${product}/${version}/pwd/lardata_Utilities )

install_headers()
install_fhicl()
install_source()

