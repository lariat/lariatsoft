cet_enable_asserts()

art_make( NO_PLUGINS
          LIBRARY_NAME  LArIATDataProviders
          LIB_LIBRARIES
  	                wda
	                IFDatabase
                        larcore_Geometry
                        ${MF_MESSAGELOGGER}
                        ${MF_UTILITIES}
                        ${CETLIB}
                        ${FHICLCPP}
                        ${ROOT_BASIC_LIB_LIST}
                        ${Boost_SYSTEM_LIBRARY}
        )



install_headers()
install_fhicl()
install_source()
