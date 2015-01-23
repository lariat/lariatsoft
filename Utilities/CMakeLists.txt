cet_enable_asserts()

set( lariat_util_lib_list LArFFT_service
                        LArProperties_service
                        Geometry
                        HitFinder
                        RecoBase
                        Geometry_service
                        Utilities
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_PERSISTENCY_COMMON}
                        ${ART_PERSISTENCY_PROVENANCE}
                        ${ART_FRAMEWORK_SERVICES_BASIC}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        ${ART_FRAMEWORK_SERVICES_SYSTEM_FILECATALOGMETADATA_SERVICE}
                        ${ART_UTILITIES}
                        ${MF_MESSAGELOGGER}
                        ${MF_UTILITIES}
                        ${CETLIB}
                        ${ROOT_BASIC_LIB_LIST}
    )




simple_plugin( FileCatalogMetadataLArIAT "service"
               ${lariat_util_lib_list}
	       BASENAME_ONLY
	)
             
 
install_headers()
# do not install duplicate standard_reco_uboone.fcl file
# hack until install_fhicl( EXCLUDES standard_reco_uboone.fcl ) does the right thing
install_fhicl( LIST sam_lariat.fcl
                    services_lariat.fcl
                    )
install_source()


