# This CMakeLists.txt should be invoked via dftbp_config_hybrid_dependency()!

# set(BUILD_EXPORTED_TARGETS_ONLY TRUE)

get_target_property(LAPACK_LIBRARIES LAPACK::LAPACK INTERFACE_LINK_LIBRARIES)
add_subdirectory(${OPENMMPOL_SOURCE_DIR} ${OPENMMPOL_BINARY_DIR})

add_library(openmmpol::openmmpol INTERFACE IMPORTED GLOBAL)
target_link_libraries(openmmpol::openmmpol INTERFACE openmmpol)