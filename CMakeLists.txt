PROJECT (atus2)

cmake_minimum_required(VERSION 2.8)

# specify compiler
# cmake -DCMAKE_CXX_COMPILER=g++-6 -DCMAKE_C_COMPILER=gcc-6 ..

find_package(MPI)
find_package(VTK NO_MODULE)
find_package(Boost REQUIRED)
find_package(Doxygen)

#SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})

set(CMAKE_CXX_FLAGS_RELEASE "-std=gnu++14 -march=native -O3 -funroll-loops -ftree-vectorize -fopenmp -w -s -Wall")
#set(CMAKE_CXX_FLAGS_DEBUG "-std=gnu++14 -g -Wall -Wextra -fopenmp -fsanitize=thread")
#set(CMAKE_CXX_FLAGS_DEBUG "-std=gnu++14 -g -Wall -Wextra -fopenmp -fsanitize=address")
set(CMAKE_CXX_FLAGS_DEBUG "-std=gnu++14 -g -Wall -Wextra -fopenmp -Wno-unused-parameter -Wno-unused-variable -Wno-unused-but-set-variable")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Release" "Debug" )
endif()

option( BUILD_DOCUMENTATION "Create and install the HTML based API documentation (requires Doxygen)" ${DOXYGEN_FOUND})

set( HOME $ENV{HOME} CACHE STRING INTERNAL )
set( DIR_INC ${PROJECT_SOURCE_DIR}/include CACHE STRING INTERNAL )
set( DIR_MYLIB ${PROJECT_SOURCE_DIR}/source/libs/lib_myutils CACHE STRING INTERNAL )
set( DIR_MYLIB_MPI ${PROJECT_SOURCE_DIR}/source/libs/lib_myutils_mpi_2 CACHE STRING INTERNAL )

if(BUILD_DOCUMENTATION)
    if(NOT DOXYGEN_FOUND)
        message(FATAL_ERROR "Doxygen is needed to build the documentation.")
    endif()

    set(doxyfile ${PROJECT_SOURCE_DIR}/Doxyfile)

    add_custom_target(doc
        COMMAND ${DOXYGEN_EXECUTABLE} ${doxyfile}
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM)
endif()

### GSL library
find_library( GSL_LIBRARY_1 NAMES gsl  PATHS ENV LD_LIBRARY_PATH NO_DEFAULT_PATH )
if( NOT GSL_LIBRARY_1 )
  message( FATAL_ERROR "Could not find libgsl." )
endif()
find_library( GSL_LIBRARY_2 NAMES gslcblas PATHS ENV LD_LIBRARY_PATH NO_DEFAULT_PATH )
if( NOT GSL_LIBRARY_2 )
  message( FATAL_ERROR "Could not find libgslcblas." )
endif()
get_filename_component( TMP ${GSL_LIBRARY_1} PATH )
get_filename_component( TMP ${TMP} PATH )
set( GSL_INC_DIR ${TMP}/include CACHE STRING INTERNAL )

# FFTW library
find_library( FFTW_LIBRARY_1 NAMES fftw3  PATHS ENV LD_LIBRARY_PATH NO_DEFAULT_PATH )
if( NOT FFTW_LIBRARY_1 )
  message( FATAL_ERROR "Could not find libfftw3." )
endif()
find_library( FFTW_LIBRARY_2 NAMES fftw3_omp  PATHS ENV LD_LIBRARY_PATH NO_DEFAULT_PATH )
if( NOT FFTW_LIBRARY_2 )
  message( FATAL_ERROR "Could not find libfftw3_omp." )
endif()
if( MPI_CXX_FOUND )
  find_library( FFTW_LIBRARY_3 NAMES fftw3_mpi  PATHS ENV LD_LIBRARY_PATH NO_DEFAULT_PATH )
  if( NOT FFTW_LIBRARY_3 )
    message( FATAL_ERROR "Could not find libfftw3_mpi." )
  endif()
endif()
get_filename_component( TMP ${FFTW_LIBRARY_1} PATH )
get_filename_component( TMP ${TMP} PATH )
set( FFTW_INC_DIR ${TMP}/include CACHE STRING INTERNAL )

# muParser library
find_library( MUPARSER_LIBRARY NAMES muparser  PATHS ENV LD_LIBRARY_PATH NO_DEFAULT_PATH )
if( NOT MUPARSER_LIBRARY )
  message( FATAL_ERROR "Could not find libmuparser." )
endif()
get_filename_component( TMP ${MUPARSER_LIBRARY} PATH )
get_filename_component( TMP ${TMP} PATH )
set( MUPARSER_INC_DIR ${TMP}/include CACHE STRING INTERNAL )

# lis library http://www.ssisc.org/lis/index.en.html
find_library( LIS_LIBRARY NAMES lis PATHS ENV LD_LIBRARY_PATH NO_DEFAULT_PATH )
if( LIS_LIBRARY )
  get_filename_component( TMP ${LIS_LIBRARY} PATH )
  get_filename_component( TMP ${TMP} PATH )
  set( LIS_INC_DIR ${TMP}/include CACHE STRING INTERNAL )
  include_directories( ${LIS_INC_DIR} )
endif()

set( EXECUTABLE_OUTPUT_PATH  ${HOME}/bin )
set( LIBRARY_OUTPUT_PATH  ${HOME}/bin/lib )
#SET( EXECUTABLE_INSTALL_PATH  ${HOME}/bin )

# Zuerst die Pfade angeben
#=========================
include_directories(SYSTEM ${MPI_CXX_INCLUDE_PATH} ${Boost_INCLUDE_DIRS} ${VTK_INCLUDE_DIRS}  ${GSL_INC_DIR} ${FFTW_INC_DIR} ${MUPARSER_INC_DIR} )
include_directories( ${DIR_INC} ${DIR_MYLIB} ${DIR_MYLIB_MPI}/src_mpi )
#LINK_DIRECTORIES( ${GSL_LIB_DIR} ${FFTW_LIB_DIR} )

# enable profiling
#set( CMAKE_EXE_LINKER_FLAGS -pg )

add_subdirectory( source/libs/lib_myutils )
add_subdirectory( source/SOLVER/src )
add_subdirectory( source/utils/src )
add_subdirectory( source/vis/src )
#add_subdirectory( source/misc )
add_subdirectory( source/interferometry/src )
#add_subdirectory( source_students/tammy/src )

if( LIS_LIBRARY )
  add_subdirectory( source_students/frank/src )
endif()

if( MPI_CXX_FOUND )
  add_subdirectory( ${DIR_MYLIB_MPI}/src_mpi )
  add_subdirectory( source/SOLVER/src_mpi )
  add_subdirectory( source/interferometry/src_mpi )
endif()

# FILE(GLOB bash_sh "${PROJECT_SOURCE_DIR}/Bash/*")
#
# FOREACH( file_i ${bash_sh})
#     MESSAGE(STATUS ${file_i} )
#     INSTALL(FILES ${file_i} PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ DESTINATION "${HOME}/bin" )
# ENDFOREACH( file_i )
