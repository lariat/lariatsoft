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
cmake_minimum_required (VERSION 2.8)

project(lariatsoft)

# cetbuildtools contains our cmake modules
SET ( CETBUILDTOOLS_VERSION $ENV{CETBUILDTOOLS_VERSION} )
IF (NOT CETBUILDTOOLS_VERSION)
    MESSAGE (FATAL_ERROR "ERROR: setup cetbuildtools to get the cmake modules")
ENDIF()

set(CMAKE_MODULE_PATH $ENV{ART_DIR}/Modules
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
find_ups_product( larsoft         v03_04_05 )
find_ups_product( larcore         v03_03_01 )
find_ups_product( lardata         v03_03_01 )
find_ups_product( larevt          v03_02_04 )
find_ups_product( larsim          v03_02_04 )
find_ups_product( larreco         v03_03_01 )
find_ups_product( larana          v03_03_01 )
find_ups_product( larpandora      v03_03_03 )
find_ups_product( lareventdisplay v03_02_04 )

# art related products
find_ups_product( art             v1_12_02 )
find_ups_product( artdaq_core     v1_04_06 )

# products that LArIAT needs for input data
find_ups_product( LariatFilters   v01_00    )
find_ups_product( lariatfragments v01_00_00 )

# non-art library dependences
find_ups_product( postgresql v9_1_5 )
find_ups_product( nutools v1_00_00 )
find_ups_root( root v5_34_18 )

# macros for artdaq_dictionary and simple_plugin
include(ArtDictionary)
include(ArtMake)
include(BuildPlugins)

# source
add_subdirectory(Calibration)
add_subdirectory(DetectorSim)
add_subdirectory(Geo)
add_subdirectory(JobConfigurations)

add_subdirectory(LArSoftConfigurations)


add_subdirectory(LArIATDataProducts)
add_subdirectory(LArIATEventGenerator)
add_subdirectory(LArIATRecoAlg)
add_subdirectory(LArIATRecoModule)
add_subdirectory(LArIATChecks)
add_subdirectory(LArIATAnaModule)

add_subdirectory(RawDataUtilities)
add_subdirectory(Utilities)
add_subdirectory(mccconfigs)
add_subdirectory(LongBo)


# tests
add_subdirectory(test)

# ups - table and config files
add_subdirectory(ups)

# packaging utility
include(UseCPack)
