################################################################################
# MULTI
################################################################################

cmake_minimum_required( VERSION 3.12 )

project( multi LANGUAGES Fortran )

## Ecbuild integration
find_package( ecbuild QUIET )
include( ecbuild_system NO_POLICY_SCOPE )
ecbuild_declare_project()
list( APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake )

## Dependencies
find_package( NetCDF REQUIRED COMPONENTS Fortran )
find_package( LAPACK REQUIRED )
find_package( FFTW REQUIRED )

## Main executable
file (GLOB MULTI_SRC
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.f90" )
ecbuild_add_executable( TARGET  multi
                        SOURCES ${MULTI_SRC})
target_link_libraries( ${PROJECT_NAME} PUBLIC NetCDF::NetCDF_Fortran )
target_link_libraries( ${PROJECT_NAME} PUBLIC ${LAPACK_LIBRARIES} )
target_link_libraries( ${PROJECT_NAME} PUBLIC ${FFTW_LIBRARIES} )
include_directories( ${FFTW_INCLUDE_DIRS} )
