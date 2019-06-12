# ======================================================================
#  lariatsoft main build file
#
#  cd .../path/to/build/directory
#  source .../path/to/lariatsoft/ups/setup_for_development <-d|-p>
#  cmake [-DCMAKE_INSTALL_PREFIX=/install/path]
#        -DCMAKE_BUILD_TYPE=$CETPKG_TYPE
#        $CETPKG_SOURCE
#  make
#  make test
#  make install
#  make package (builds distribution tarfile)
# ======================================================================


# use cmake 2.8 or later
cmake_minimum_required (VERSION 3.3)

project(lariatsoft)

# cetbuildtools contains our cmake modules
SET ( CETBUILDTOOLS_VERSION $ENV{CETBUILDTOOLS_VERSION} )
IF (NOT CETBUILDTOOLS_VERSION)
    MESSAGE (FATAL_ERROR "ERROR: setup cetbuildtools to get the cmake modules")
ENDIF()

set(CMAKE_MODULE_PATH $ENV{CANVAS_ROOT_IO_DIR}/Modules
		      $ENV{CETBUILDTOOLS_DIR}/Modules
		      ${CMAKE_MODULE_PATH})

include(CetCMakeEnv)
cet_cmake_env()

set_install_root()

cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs
  EXTRA_CXX_FLAGS
)

cet_report_compiler_flags()

# these are minimum required versions, not the actual product versions

# these are the larsoft products 
find_ups_product( larsoft )
find_ups_product( larcoreobj )
find_ups_product( larcore )
find_ups_product( lardataobj )
find_ups_product( lardataalg )
find_ups_product( lardata )
find_ups_product( larevt )
find_ups_product( larsim )
find_ups_product( larreco )
find_ups_product( larana )
find_ups_product( larpandora )
find_ups_product( lareventdisplay )
find_ups_product( cetbuildtools )


# art related products
find_ups_product( art )
find_ups_product(art_root_io)
find_ups_product( artdaq_core )

# products that LArIAT needs for input data
find_ups_product( LariatFilters )
find_ups_product( lariatfragments )

# non-art library dependences
find_ups_product( postgresql )
find_ups_product( nutools )
find_ups_product( nurandom )
find_ups_root( root )
find_ups_geant4( geant4 )

# macros for artdaq_dictionary and simple_plugin
include(ArtDictionary)
include(ArtMake)
include(BuildPlugins)

# source
add_subdirectory(Calibration)
add_subdirectory(DetectorSim)
add_subdirectory(Geo)
add_subdirectory(LArIATAnaModule)
add_subdirectory(LArIATDataProducts)
add_subdirectory(LArIATDataProviders)
add_subdirectory(LArIATEventGenerator)
add_subdirectory(LArIATFilterModule)
add_subdirectory(LArIATRecoAlg)
add_subdirectory(LArIATRecoModule)
add_subdirectory(LArIATChecks)
add_subdirectory(LongBo)
add_subdirectory(RawDataUtilities)
add_subdirectory(Utilities)

# fcl files
add_subdirectory(JobConfigurations)
add_subdirectory(LArSoftConfigurations)
add_subdirectory(mccconfigs)

# tests
add_subdirectory(test)

# ups - table and config files
add_subdirectory(ups)

# packaging utility
include(UseCPack)
