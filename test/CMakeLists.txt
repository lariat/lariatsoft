# various integration tests

include(CetTest)

cet_enable_asserts()

# library unit tests
add_subdirectory(Geo)

# integration tests
add_subdirectory(ci)

