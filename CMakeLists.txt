# Set your Fortran/C/C++ compiler (GNU or Intel):
set (COMPILER GNU)
#
# Set your build type (DEBUG or RELEASE):
set (BUILD_TYPE RELEASE)
#
# Set FFTW3 directories
set (EXTERN_FFTW3_INCLUDE /usr/include)
set (EXTERN_FFTW3_LIBPATH /usr/lib/x86_64-linux-gnu)
#
# Set LAPACK directories
set (EXTERN_LAPACK_INCLUDE /usr/include)
set (EXTERN_LAPACK_LIBPATH /usr/lib/x86_64-linux-gnu)
#
#----------------------------------------------------------------------
# That's all, do not modify what is below
#----------------------------------------------------------------------
#
#----------------------------------------------------------------------
# Check
#----------------------------------------------------------------------
#
# Check compiler
if (NOT ${COMPILER} MATCHES GNU AND NOT ${COMPILER} MATCHES Intel)
   message (FATAL_ERROR  "Wrong compiler, abort!")
endif (NOT ${COMPILER} MATCHES GNU AND NOT ${COMPILER} MATCHES Intel)
#
# Check build type
if (NOT ${BUILD_TYPE} MATCHES DEBUG AND NOT ${BUILD_TYPE} MATCHES RELEASE)
   message (FATAL_ERROR  "Wrong build type, abort!")
endif (NOT ${BUILD_TYPE} MATCHES DEBUG AND NOT ${BUILD_TYPE} MATCHES RELEASE)
#
# Check CMAKE version
cmake_minimum_required (VERSION 2.8)
#
#----------------------------------------------------------------------
# Setup compiler
#----------------------------------------------------------------------
#
# GNU compiler
if (${COMPILER} MATCHES "GNU")
   set (CMAKE_Fortran_COMPILER gfortran)
   set (CMAKE_C_COMPILER gcc)
   set (CMAKE_CXX_COMPILER g++)
endif (${COMPILER} MATCHES "GNU")
#
# Intel compiler
if (${COMPILER} MATCHES "Intel")
   set (CMAKE_Fortran_COMPILER ifort)
   set (CMAKE_C_COMPILER icc)
   set (CMAKE_C++_COMPILER icpc)
endif (${COMPILER} MATCHES "Intel")
#
#----------------------------------------------------------------------
# Setup project
#----------------------------------------------------------------------
#
cmake_minimum_required (VERSION 2.6)
project (multi Fortran)
file (GLOB MULTI_SRC
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.f90"
)
add_executable (multi ${MULTI_SRC})
include_directories (${EXTERN_FFTW3_INCLUDE} ${EXTERN_LAPACK_INCLUDE})
target_link_libraries(multi fftw3 lapack)
#
#----------------------------------------------------------------------
# Define compiler flags
#----------------------------------------------------------------------
#
# GNU compiler
if (${COMPILER} MATCHES "GNU")
   set (CMAKE_Fortran_FLAGS_RELEASE "-g -funroll-all-loops -O3 -ffree-line-length-none -fdefault-real-8 -fdefault-double-8 -fopenmp")
   set (CMAKE_Fortran_FLAGS_DEBUG   "-g -Og -Wextra -Wall -pedantic -fbacktrace -ftrapv -fall-intrinsics -fcheck=all -fimplicit-none -ffree-line-length-none -fdefault-real-8 -fdefault-double-8 -ffpe-trap=invalid,zero,overflow -fopenmp")
endif (${COMPILER} MATCHES "GNU")
#
# Set linker flags
set (LINKER_FLAGS "-L${EXTERN_FFTW3_LIBPATH} -Wl,-rpath,${EXTERN_FFTW3_LIBPATH} -L${EXTERN_LAPACK_LIBPATH} -Wl,-rpath,${EXTERN_LAPACK_LIBPATH} -lfftw3 -llapack")
set_target_properties (multi PROPERTIES LINKER_LANGUAGE Fortran)
set (CMAKE_EXE_LINKER_FLAGS ${LINKER_FLAGS})
#
# Set build type
set (CMAKE_BUILD_TYPE ${BUILD_TYPE} CACHE STRING "Choose the type of build, options are: Debug/Release" FORCE)
#
#----------------------------------------------------------------------
# Define building directories
#----------------------------------------------------------------------
#
set (EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/run)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/run)
