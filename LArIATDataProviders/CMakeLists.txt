cet_enable_asserts()

art_make( NO_PLUGINS
          LIBRARY_NAME  LArIATDataProviders
          LIB_LIBRARIES
  	                wda
	                nuevdb_IFDatabase
                        larcorealg_Geometry
                        ${MF_MESSAGELOGGER}
                        
                        cetlib cetlib_except
                        ${FHICLCPP}
                        ${ROOT_BASIC_LIB_LIST}
                        ${Boost_SYSTEM_LIBRARY}
        )



install_headers()
install_fhicl()
install_source()
