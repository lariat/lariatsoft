
# geometry test on LArIAT geometry

# thid test runs a FCL file and fails only if lar exits with a bad exit code;
# it is fast and runs in all test suites.
# FCL files need to be copied to the test area (DATAFILES directive) since they
# are not installed.
cet_test(geometry_lariat HANDBUILT
  DATAFILES test_geometry_lariat.fcl
  TEST_EXEC lar
  TEST_ARGS --rethrow-all --config ./test_geometry_lariat.fcl
)


# unit tests (use the hard-coded configuration for LArIAT geometry)
#
# the following test is equivalent to the module-based one above, but with no framework
#
cet_test(geometry_lariat_test
  SOURCES geometry_lariat_test.cxx
  DATAFILES test_geometry_lariat.fcl
  TEST_ARGS ./test_geometry_lariat.fcl
  LIBRARIES larcore_Geometry
            GeometryTestLib
            Geo # this is LArIAT's geometry
            ${MF_MESSAGELOGGER}
            ${MF_UTILITIES}
            ${FHICLCPP}
            ${CETLIB}
)


# this uses BOOST for the test
cet_test(geometry_iterator_lariat_test
  SOURCES geometry_iterator_lariat_test.cxx
  LIBRARIES larcore_Geometry
            GeometryTestLib
            Geo # this is LArIAT's geometry
            ${MF_MESSAGELOGGER}
            ${MF_UTILITIES}
            ${FHICLCPP}
            ${CETLIB}
  USE_BOOST_UNIT
)


cet_test(geometry_iterator_loop_lariat_test
  SOURCES geometry_iterator_loop_lariat_test.cxx
  LIBRARIES larcore_Geometry
            GeometryTestLib
            Geo # this is LArIAT's geometry
            ${MF_MESSAGELOGGER}
            ${MF_UTILITIES}
            ${FHICLCPP}
            ${CETLIB}
)


# test for channel mapping (BOOST unit test)
cet_test(geometry_lariatchannelmapping_test
  SOURCES geometry_lariatchannelmapping_test.cxx
  LIBRARIES larcore_Geometry
            GeometryTestLib
            Geo # this is LArIAT's geometry
            ${MF_MESSAGELOGGER}
            ${MF_UTILITIES}
            ${FHICLCPP}
            ${CETLIB}
  USE_BOOST_UNIT
  DATAFILES test_geometry_lariat.fcl
  TEST_ARGS ./test_geometry_lariat.fcl
)



