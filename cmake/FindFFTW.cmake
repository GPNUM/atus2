
# Find the FFTW library
# Usage:
#   find_package( FFTW [REQUIRED] [QUIET] )
#

if( FFTW_ROOT )

    find_library(  FFTW_LIBRARY_1
                   NAMES "fftw3" 
                   PATHS ${FFTW_ROOT}
                   PATH_SUFFIXES "lib" "lib64"
                   NO_DEFAULT_PATH
    )

    find_library(  FFTW_LIBRARY_2
                   NAMES "fftw3_omp"
                   PATHS ${FFTW_ROOT}
                   PATH_SUFFIXES "lib" "lib64"
                   NO_DEFAULT_PATH
    )

    find_library(  FFTW_LIBRARY_3
                   NAMES "fftw3_mpi"
                   PATHS ${FFTW_ROOT}
                   PATH_SUFFIXES "lib" "lib64"
                   NO_DEFAULT_PATH
    )    

    find_path(  FFTW_INCLUDE_DIR
                NAMES "fftw3.h"
                PATHS ${FFTW_ROOT}
                PATH_SUFFIXES "include"
                NO_DEFAULT_PATH
    )

else()

    find_library(   FFTW_LIBRARY_1 
                    NAMES "fftw3"  
                    PATHS ENV LD_LIBRARY_PATH NO_DEFAULT_PATH 
    )

    find_library(   FFTW_LIBRARY_2 
                    NAMES "fftw3_omp"  
                    PATHS ENV LD_LIBRARY_PATH NO_DEFAULT_PATH 
    )

    find_library(   FFTW_LIBRARY_3
                    NAMES "fftw3_mpi"  
                    PATHS ENV LD_LIBRARY_PATH NO_DEFAULT_PATH 
    )    

    get_filename_component( TMP ${FFTW_LIBRARY_1} PATH )
    get_filename_component( TMP ${TMP} PATH )
    set( FFTW_INCLUDE_DIR ${TMP}/include CACHE STRING INTERNAL )

endif()

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(FFTW
      REQUIRED_VARS FFTW_INCLUDE_DIR FFTW_LIBRARY_1 FFTW_LIBRARY_2 FFTW_LIBRARY_3
      HANDLE_COMPONENTS
      )

mark_as_advanced(
      FFTW_LIBRARY_1
      FFTW_LIBRARY_2
      FFTW_LIBRARY_3
      FFTW_INCLUDE_DIR
      )
